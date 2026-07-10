! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! Strongly-coupled (implicit) FSI scheme using sub-iteration.
module fluid_pnpn_fsi_subiteration
  use fsi_dynamics, only : fsi_body_t, assemble_structural_inertial_terms, &
       add_fsi_non_linear_matrices
  use fsi_manager, only : fsi_manager_init, linsolve_dense, &
       fsi_prep_checkpoint, fsi_restart_restore
  use fluid_pnpn, only : fluid_pnpn_t
  use field, only : field_t
  use field_math, only : field_copy
  use num_types, only : rp, dp
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : HOST_TO_DEVICE, device_sync
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use projection, only : projection_t
  use projection_vel, only : projection_vel_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  use utils, only : neko_error
  use logger, only : neko_log, LOG_SIZE
  use mesh, only : mesh_t
  use user_intf, only : user_t
  use checkpoint, only : chkp_t
  use mpi_f08, only : MPI_Wtime
  use math, only : rzero, copy
  use device_math, only : device_copy
  use dofmap, only : dofmap_t
  use ale_manager
  implicit none
  private

  ! Coupling-acceleration ("relaxation") methods
  integer, parameter :: ACCEL_CONSTANT = 0   !< fixed under-relaxation
  integer, parameter :: ACCEL_AITKEN = 1   !< per-body dynamic Aitken
  integer, parameter :: ACCEL_IQN = 2   !< IQN-ILS (@ToDo)

  ! Floor for the per-DOF relative error: a DOF whose own velocity/force is
  ! below REL_FLOOR times the peak of its type (translation/rotation, force/
  ! torque) is normalised by that floor instead of its own (near-zero) scale,
  ! so a momentarily-tiny DOF neither blows up nor hides. Purely a
  ! normalisation choice; does not affect the converged solution.
  real(kind=rp), parameter :: REL_FLOOR = 1.0e-3_rp

  ! Predictor for the sub-iteration: none = v_s := v^n (zeroth order)
  !, ext = EXT-k extrapolation of the body velocity history.
  integer, parameter :: PRED_NONE = 0
  integer, parameter :: PRED_EXT  = 1

  type, public, extends(fluid_pnpn_t) :: fluid_pnpn_fsi_subiter_t
     logical :: if_fsi = .false.
     logical :: res_long_print = .false.
     logical :: non_linear_correction_term = .false.
     real(kind=rp) :: gravity_vec(3) = 0.0_rp

     !> Structure integrator: .false. => BDF structure + BDF ALE (default),
     !> .true. => Newmark (constant average acceleration) structure + CN ALE.
     logical :: structure_newmark = .false.
     !> Newmark initial acceleration a^0: .false. => a^0 = 0 (default),
     !> .true. => a^0 from the initial force balance.
     logical :: newmark_a0_balance = .false.
     !> First-step a^0 initialisation done (lazy, needs dt).
     logical :: newmark_a0_done = .false.

     ! N-body storage
     integer :: nbodies_fsi = 0
     type(fsi_body_t), allocatable :: fsi_bodies(:)

     ! Global dense structural system  M * dv = B
     integer, allocatable :: fsi_dof_map(:,:)
     integer :: total_active_dofs = 0
     real(kind=rp), allocatable :: M_global(:,:)
     real(kind=rp), allocatable :: B_global(:)
     real(kind=rp), allocatable :: X_sol(:)
     real(kind=rp), allocatable :: rel_dof_v(:,:)
     real(kind=rp), allocatable :: rel_dof_f(:,:)

     ! Batch arrays for the ALE mesh-velocity override
     integer, allocatable :: batch_ids(:)
     real(kind=rp), allocatable :: batch_trans(:,:)
     real(kind=rp), allocatable :: batch_ang(:,:)
     real(kind=rp), allocatable :: temp_prescribed_vels(:,:)

     ! Checkpointing (managed by fsi_manager)
     real(kind=rp), allocatable :: global_disp_rel(:)
     real(kind=rp), allocatable :: global_body_vel(:)
     !> Newmark previous-acceleration checkpoint store (sub-iteration + Newmark).
     real(kind=rp), allocatable :: global_body_acc(:)
     !> Newmark prescribed-frame previous acceleration (a_f^n), packed for the
     !> checkpoint. Not reconstructable from the saved lags, so it is stored.
     real(kind=rp), allocatable :: global_frame_acc(:)
     real(kind=rp), allocatable :: global_body_vel_lag(:,:)
     real(kind=rp), allocatable :: global_moving_frame_presc_vel(:,:)

     ! Unused Green's-function storage (only present to satisfy fsi_manager_init).
     logical :: skip_greens_solve = .true.
     type(field_t), allocatable :: u_g(:), v_g(:), w_g(:), p_g(:)
     type(projection_t), allocatable :: proj_prs_green(:)
     type(projection_vel_t), allocatable :: proj_vel_green(:)

     integer :: max_subiter
     integer :: min_subiter
     real(kind=rp) :: relax

     integer :: accel_method = ACCEL_AITKEN !< constant | aitken | iqn
     real(kind=rp), allocatable :: accel_omega(:)   !< per-body omega (nbodies)
     real(kind=rp), allocatable :: accel_r_prev(:,:)!< prev increment (6,nbodies)
     real(kind=rp) :: aitken_max
     real(kind=rp) :: aitken_min
     character(len=16) :: conv_criterion !< velocity|force|both|either
     logical :: conv_relative
     real(kind=rp) :: conv_tol

     integer :: predictor_guess = PRED_NONE  !< none | ext
     integer :: predictor_order = 0  !< 0 = advection ramp; >0 = fixed

     ! fluid reference snapshot
     type(field_t) :: u_ref, v_ref, w_ref, p_ref
     type(field_t) :: fx_ref, fy_ref, fz_ref
     type(field_t) :: ue_ref, ve_ref, we_ref

     ! BDF-k position histories (frozen during sub-iterations)
     integer :: n_mesh = 0
     integer :: n_lag = 0
     type(field_t), allocatable :: mesh_x_lag(:), mesh_y_lag(:), mesh_z_lag(:)

     ! Rigid-body position histories
     real(kind=rp), allocatable :: pivot_hist(:,:,:)   ! (3, n_lag, nbodies)
     real(kind=rp), allocatable :: ghost_hist(:,:,:,:) ! (3, n_lag, 2, nbodies)
     real(kind=rp), allocatable :: disp_hist(:,:,:)    ! (6, n_lag, nbodies)

   contains
     procedure, pass(this) :: init => fluid_subiter_init
     procedure, pass(this) :: step => fluid_subiter_step
     procedure, pass(this) :: free => fluid_subiter_free
     procedure, pass(this) :: restart => fluid_subiter_restart
     procedure, pass(this) :: calc_fsi_terms => subiter_calc_fsi_terms
     procedure, pass(this) :: query_frame_prescribed_motion => &
          subiter_query_frame_prescribed_motion
     procedure, pass(this) :: log_fsi_results => subiter_log_results
     procedure, pass(this) :: snapshot_fluid => subiter_snapshot_fluid
     procedure, pass(this) :: restore_fluid  => subiter_restore_fluid
     procedure, pass(this) :: init_histories => subiter_init_histories
     procedure, pass(this) :: update_kinematics_bdfk => subiter_update_kinematics
     procedure, pass(this) :: commit_histories => subiter_commit_histories
     procedure, pass(this) :: add_spring_jacobian => subiter_add_spring_jacobian
  end type fluid_pnpn_fsi_subiter_t

contains

  subroutine fluid_subiter_init(this, msh, lx, params, user, chkp)
    class(fluid_pnpn_fsi_subiter_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    type(time_state_t) :: t_init
    integer :: i
    character(len=32) :: accel_str
    logical :: old_aitken
    character(:), allocatable :: tmp_str
    ! Buffers for the sub-iteration setup summary (logged at the end of init).
    character(len=128) :: log_buf
    character(len=16) :: scheme_name, ale_name, accel_name, pred_name, a0_name
  
    ! Initialize the base fluid_pnpn scheme (allocates the mesh and ALE scheme)
    call this%fluid_pnpn_t%init(msh, lx, params, user, chkp)

    ! Structure time integrator: 'bdf' (default) or 'newmark'. Selecting
    ! 'newmark' automatically switches the ALE mesh update to CN, so the user
    ! cannot pick an inconsistent structure/mesh pair. Read here (before the
    ! FSI manager / checkpoint registration) so the acceleration checkpoint
    ! store is registered only for the Newmark path.
    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.structure_scheme', tmp_str, 'bdf')
    select case (trim(tmp_str))
    case ('bdf')
       this%structure_newmark = .false.
    case ('newmark')
       this%structure_newmark = .true.
    case default
       call neko_error('Unknown structure_scheme: ' // trim(tmp_str) // &
            ' (use bdf|newmark)')
    end select

    ! Newmark initial acceleration a^0: 'zero' (default) or 'balance'.
    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.newmark_initial_acceleration', &
         tmp_str, 'zero')
    select case (trim(tmp_str))
    case ('zero')
       this%newmark_a0_balance = .false.
    case ('balance')
       this%newmark_a0_balance = .true.
    case default
       call neko_error('Unknown newmark_initial_acceleration: ' // &
            trim(tmp_str) // ' (use zero|balance)')
    end select

    ! Initialize the FSI manager
    call fsi_manager_init(params, msh, this%ale, this%c_Xh, this%dm_Xh, &
         this%if_fsi, this%nbodies_fsi, this%fsi_bodies, this%fsi_dof_map, &
         this%total_active_dofs, &
         this%M_global, this%B_global, this%X_sol, &
         this%u_g, this%v_g, this%w_g, this%p_g, &
         this%res_long_print, this%gravity_vec, this%proj_prs_green, &
         this%proj_vel_green, this%global_disp_rel, &
         this%global_body_vel, this%global_body_vel_lag, &
         this%global_moving_frame_presc_vel, this%skip_greens_solve, &
         this%non_linear_correction_term, &
         global_body_acc = this%global_body_acc, &
         global_frame_acc = this%global_frame_acc)
         
    this%skip_greens_solve = .true.

    ! Register the acceleration stores for checkpointing only for the Newmark
    ! path, so BDF checkpoints keep their exact byte layout.
    if (this%structure_newmark) then
       call this%chkp%add_fsi(this%global_disp_rel, this%global_body_vel, &
            this%global_body_vel_lag, this%global_moving_frame_presc_vel, &
            body_acc = this%global_body_acc, &
            frame_acc = this%global_frame_acc)
    else
       call this%chkp%add_fsi(this%global_disp_rel, this%global_body_vel, &
            this%global_body_vel_lag, this%global_moving_frame_presc_vel)
    end if

    ! Sub-iteration parameters.
    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.max_iterations', this%max_subiter, 20)
    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.min_iterations', this%min_subiter, 1)

    if (params%valid_path( &
         'case.fluid.fsi.subiteration.coupling_acceleration.method')) then
       call json_get_or_default(params, &
            'case.fluid.fsi.subiteration.coupling_acceleration.method', &
            tmp_str, 'aitken')
       accel_str = tmp_str
       call json_get_or_default(params, &
            'case.fluid.fsi.subiteration.coupling_acceleration.relaxation_value',&
            this%relax, 0.5_rp)
       call json_get_or_default(params, &
            'case.fluid.fsi.subiteration.coupling_acceleration.aitken_min', &
            this%aitken_min, 0.05_rp)
       call json_get_or_default(params, &
            'case.fluid.fsi.subiteration.coupling_acceleration.aitken_max', &
            this%aitken_max, 1.0_rp)
    end if

    select case (trim(accel_str))
    case ('constant')
       this%accel_method = ACCEL_CONSTANT
    case ('aitken')
       this%accel_method = ACCEL_AITKEN
    case ('iqn')
       this%accel_method = ACCEL_IQN
    case default
       call neko_error('Unknown coupling_acceleration method: ' // &
            trim(accel_str) // ' (use constant|aitken|iqn)')
    end select
    if (this%accel_method == ACCEL_IQN) then
       call neko_error('coupling_acceleration method "iqn" is not yet ' // &
            'implemented; use "constant" or "aitken".')
    end if

    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.criterion', tmp_str, 'velocity')
    this%conv_criterion = tmp_str
    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.relative', this%conv_relative, .true.)
    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.tolerance', this%conv_tol, 1.0e-6_rp)

    ! predictor_guess = none | ext ; predictor_order = 0 (advection ramp) or 1-3.
    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.predictor_guess', tmp_str, 'none')
    select case (trim(tmp_str))
    case ('none')
       this%predictor_guess = PRED_NONE
    case ('ext')
       this%predictor_guess = PRED_EXT
    case default
       call neko_error('Unknown predictor_guess: ' // trim(tmp_str) // &
            ' (use none|ext)')
    end select
    call json_get_or_default(params, &
         'case.fluid.fsi.subiteration.predictor_order', this%predictor_order, 0)
    if (this%predictor_order < 0 .or. this%predictor_order > 3) then
       call neko_error('predictor_order must be 0 (advection ramp) or 1-3.')
    end if

    ! -----------------------------------------------------------------------
    ! Resolved sub-iteration setup summary. Logged once so the run can be
    ! reproduced from the log file (each line maps to a JSON key under
    ! case.fluid.fsi.subiteration).
    ! -----------------------------------------------------------------------
    if (this%structure_newmark) then
       scheme_name = 'newmark'
       ale_name = 'cn'
    else
       scheme_name = 'bdf'
       ale_name = 'bdf'
    end if
    if (this%newmark_a0_balance) then
       a0_name = 'balance'
    else
       a0_name = 'zero'
    end if
    select case (this%accel_method)
    case (ACCEL_CONSTANT)
       accel_name = 'constant'
    case (ACCEL_AITKEN)
       accel_name = 'aitken'
    case (ACCEL_IQN)
       accel_name = 'iqn'
    case default
       accel_name = 'unknown'
    end select
    select case (this%predictor_guess)
    case (PRED_NONE)
       pred_name = 'none'
    case (PRED_EXT)
       pred_name = 'ext'
    case default
       pred_name = 'unknown'
    end select

    call neko_log%section('FSI sub-iteration setup')

    write(log_buf, '(A,A)') '   structure_scheme            : ', trim(scheme_name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') '   ale_scheme (auto-selected)  : ', trim(ale_name)
    call neko_log%message(log_buf)
    if (this%structure_newmark) then
       write(log_buf, '(A,A)') '   newmark_initial_acceleration: ', trim(a0_name)
       call neko_log%message(log_buf)
    end if

    write(log_buf, '(A,I0)') '   max_iterations              : ', this%max_subiter
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') '   min_iterations              : ', this%min_subiter
    call neko_log%message(log_buf)

    write(log_buf, '(A,A)') '   coupling_acceleration.method: ', trim(accel_name)
    call neko_log%message(log_buf)
    if (this%accel_method == ACCEL_AITKEN) then
       write(log_buf, '(A,ES12.5)') &
            '     relaxation_value          : ', this%relax
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES12.5)') &
            '     aitken_min                : ', this%aitken_min
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES12.5)') &
            '     aitken_max                : ', this%aitken_max
       call neko_log%message(log_buf)
    else if (this%accel_method == ACCEL_CONSTANT) then
       write(log_buf, '(A,ES12.5)') &
            '     relaxation_value          : ', this%relax
       call neko_log%message(log_buf)
    end if

    write(log_buf, '(A,A)') '   criterion                   : ', &
         trim(this%conv_criterion)
    call neko_log%message(log_buf)
    write(log_buf, '(A,L1)') '   relative                    : ', &
         this%conv_relative
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES12.5)') '   tolerance                   : ', &
         this%conv_tol
    call neko_log%message(log_buf)

    write(log_buf, '(A,A)') '   predictor_guess             : ', trim(pred_name)
    call neko_log%message(log_buf)
    if (this%predictor_guess .eq. PRED_EXT) then
       write(log_buf, '(A,I0)') '   predictor_order             : ', &
            this%predictor_order
       call neko_log%message(log_buf)
    end if

    call neko_log%end_section()

    ! Reference / scratch fluid fields.
    call this%u_ref%init(this%dm_Xh, 'u_ref')
    call this%v_ref%init(this%dm_Xh, 'v_ref')
    call this%w_ref%init(this%dm_Xh, 'w_ref')
    call this%p_ref%init(this%dm_Xh, 'p_ref')
    call this%fx_ref%init(this%dm_Xh, 'fx_ref')
    call this%fy_ref%init(this%dm_Xh, 'fy_ref')
    call this%fz_ref%init(this%dm_Xh, 'fz_ref')
    call this%ue_ref%init(this%dm_Xh, 'ue_ref')
    call this%ve_ref%init(this%dm_Xh, 've_ref')
    call this%we_ref%init(this%dm_Xh, 'we_ref')

    if (this%nbodies_fsi > 0) then
       allocate(this%batch_ids(this%nbodies_fsi))
       allocate(this%batch_trans(3, this%nbodies_fsi))
       allocate(this%batch_ang(3, this%nbodies_fsi))
       allocate(this%temp_prescribed_vels(6, this%nbodies_fsi))
       this%temp_prescribed_vels = 0.0_rp

       ! Per-body relaxation factor and previous-increment store.
       allocate(this%accel_omega(this%nbodies_fsi))
       allocate(this%accel_r_prev(6, this%nbodies_fsi))
       this%accel_omega = this%relax
       this%accel_r_prev = 0.0_rp

       allocate(this%rel_dof_v(6, this%nbodies_fsi))
       allocate(this%rel_dof_f(6, this%nbodies_fsi))
       this%rel_dof_v = -1.0_rp
       this%rel_dof_f = -1.0_rp

       ! Position-history allocation and initialisation.
       call this%init_histories()

       ! Register the position histories for checkpointing.
       call this%chkp%add_fsi_subiter(this%n_lag, &
            this%mesh_x_lag, this%mesh_y_lag, this%mesh_z_lag, &
            this%pivot_hist, this%ghost_hist, this%disp_hist)

       ! Select the implicit mesh-kinematics scheme that matches the structure
       ! integrator: BDF structure -> BDF ALE, Newmark structure -> CN ALE.
       ! The setter cross-checks params and errors if the sub-iteration scheme
       ! was not actually selected.
       if (this%structure_newmark) then
          call this%ale%select_scheme('cn', params)
       else
          call this%ale%select_scheme('bdf', params)
       end if

       ! Initial mesh velocity.
       if (.not. params%valid_path('case.restart_file')) then
          t_init%t = 0.0_rp
          t_init%tstep = 0
          t_init%dt = 0.0_rp
          do i = 1, this%nbodies_fsi
             this%batch_ids(i) = this%fsi_bodies(i)%ale_id
             this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
             this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
          end do
          call this%ale%update_mesh_velocity(this%c_Xh, t_init, &
               override_ids = this%batch_ids, &
               override_trans = this%batch_trans, &
               override_ang = this%batch_ang, &
               out_prescribed_vels = this%temp_prescribed_vels, mode = 0)
          do i = 1, this%nbodies_fsi
             this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) = &
                  this%temp_prescribed_vels(:, i)
          end do
       end if
    end if

  end subroutine fluid_subiter_init

  !> Allocate and seed the BDF-k position histories with the current state.
  subroutine subiter_init_histories(this)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    integer :: i, j, g, a, h
    character(len=16) :: idx

    this%n_mesh = this%dm_Xh%size()
    ! Number of BDF lag terms we may need.
    this%n_lag = max(this%ext_bdf%nadv, 3)

    if (.not. allocated(this%mesh_x_lag)) then
       allocate(this%mesh_x_lag(this%n_lag))
       allocate(this%mesh_y_lag(this%n_lag))
       allocate(this%mesh_z_lag(this%n_lag))
       do j = 1, this%n_lag
          write(idx, '(I0)') j
          call this%mesh_x_lag(j)%init(this%dm_Xh, 'fsi_mesh_x_lag'//trim(idx))
          call this%mesh_y_lag(j)%init(this%dm_Xh, 'fsi_mesh_y_lag'//trim(idx))
          call this%mesh_z_lag(j)%init(this%dm_Xh, 'fsi_mesh_z_lag'//trim(idx))
       end do
       allocate(this%pivot_hist(3, this%n_lag, this%nbodies_fsi))
       allocate(this%ghost_hist(3, this%n_lag, 2, this%nbodies_fsi))
       allocate(this%disp_hist(6, this%n_lag, this%nbodies_fsi))
    end if

    do j = 1, this%n_lag
       call copy(this%mesh_x_lag(j)%x, this%c_Xh%dof%x, this%n_mesh)
       call copy(this%mesh_y_lag(j)%x, this%c_Xh%dof%y, this%n_mesh)
       call copy(this%mesh_z_lag(j)%x, this%c_Xh%dof%z, this%n_mesh)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_copy(this%mesh_x_lag(j)%x_d, this%c_Xh%dof%x_d, this%n_mesh)
          call device_copy(this%mesh_y_lag(j)%x_d, this%c_Xh%dof%y_d, this%n_mesh)
          call device_copy(this%mesh_z_lag(j)%x_d, this%c_Xh%dof%z_d, this%n_mesh)
       end if
    end do

    do i = 1, this%nbodies_fsi
       a = this%fsi_bodies(i)%ale_id
       do j = 1, this%n_lag
          this%pivot_hist(:, j, i) = this%ale%ale_pivot(a)%pos
          this%disp_hist(:, j, i)  = this%fsi_bodies(i)%disp_rel
          do g = 1, 2
             h = this%ale%ghost_handles(g, a)
             this%ghost_hist(:, j, g, i) = this%ale%trackers(h)%pos
          end do
          this%fsi_bodies(i)%body_vel_lag(:, j) = this%fsi_bodies(i)%body_vel
       end do
    end do

  end subroutine subiter_init_histories

  ! Time step  (the strongly-coupled sub-iteration)
  subroutine fluid_subiter_step(this, time, dt_controller)
    class(fluid_pnpn_fsi_subiter_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller

    real(kind=rp) :: beta(0:3), dt
    integer :: nadv, n, i, k, row_g, kit
    real(kind=rp) :: F_fluid(6)
    real(kind=dp) :: t_sub0, t_sub

    ! Worst-DOF metrics and their (body, dof) index; per-DOF log helpers.
    real(kind=rp) :: vmetric, fmetric
    integer :: wv_body, wv_dof, wf_body, wf_dof
    logical :: show_v, show_f
    character(len=2) :: dof_lbl(6) = &
         ['Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz']
    character(len=LOG_SIZE) :: sub
    character(len=12) :: vtag, ftag
    character(len=24) :: tmp
    character(len=8) :: dv_lbl, fr_lbl  !< breakdown labels (rel/abs mode)

    real(kind=rp) :: omega
    logical :: converged
    character(len=128) :: msg
    real(kind=rp) :: r_body(6), dr_body(6)
    real(kind=rp) :: num, den

    real(kind=rp), allocatable :: M_linear(:,:), X_prev(:)
    real(kind=rp) :: nr_res
    integer :: nr_iter, max_nr

    dt = time%dt
    n = this%dm_Xh%size()
    nadv = this%ext_bdf%nadv

    ! beta(0) = diffusion_coeffs(1)  (= LHS coeff on v^{n+1})
    ! beta(j) = -diffusion_coeffs(j+1) (RHS-side coeffs)
    do i = 0, nadv
       beta(i) = this%ext_bdf%diffusion_coeffs%x(i+1)
       if (i .ge. 1) beta(i) = -beta(i)
    end do

    ! Freeze the previous-step mesh velocity wm^n for the CN reposition.
    ! At the top of the step this%ale%wm_x holds the converged wm from the
    ! previous step (or the initial mesh velocity on the first step). No-op
    ! for the BDF path.
    if (this%structure_newmark) call this%ale%snapshot_mesh_velocity()

    ! prescribed frame motion
    call this%query_frame_prescribed_motion(time, beta, nadv)

    ! Newmark: initialise the stored acceleration a^0 from the initial force
    ! balance on the first step (needs dt and the initial fluid force).
    if (this%structure_newmark .and. this%newmark_a0_balance &
         .and. .not. this%newmark_a0_done) then
       call subiter_newmark_init_accel(this, time, beta, nadv)
       this%newmark_a0_done = .true.
    end if

    ! predictor:  v_s := v^n (or EXT-k extrapolation)
    if (this%predictor_guess .eq. PRED_EXT) then
       call subiter_apply_ext_predictor(this, time)
    end if

    ! assemble the explicit RHS once on x^n
    ! adv%compute_ale uses the mesh velocity currently stored in the ALE manager
    ! (= last step's corrected velocity)
    ! DOUBLE CHECK TO BE SURE WM IS CORRECT HERE.
    call this%assemble_rhs(time)

    ! snapshot the fluid reference (u^n, f, u_e)
    call this%snapshot_fluid()

    ! Reset the per-body relaxation state for this step
    this%accel_omega = this%relax
    this%accel_r_prev = 0.0_rp
    omega = this%relax

    ! Convergence metric
    select case (trim(this%conv_criterion))
    case ('velocity')
       show_v = .true.  
       show_f = .false.
    case ('force')
       show_v = .false.
       show_f = .true.
    case default ! both / either
       show_v = .true.
       show_f = .true.
    end select
    vmetric = 0.0_rp
    fmetric = 0.0_rp
    wv_body = 0
    wv_dof = 0
    wf_body = 0
    wf_dof = 0
    if (this%conv_relative) then
       dv_lbl = 'dv_rel'
       fr_lbl = 'Fr_rel'
    else
       dv_lbl = 'dv_abs'
       fr_lbl = 'Fr_abs'
    end if

    call neko_log%message("---- FSI strong coupling (sub-iteration) ----")

    converged = .false.
    do kit = 1, this%max_subiter
       t_sub0 = MPI_Wtime()

       ! Implicit solve from the n-state.
       ! First sub-iteration (kit == 1): full restore to u^n so the
       ! (pressure) projection accelerates a clean solve.
       ! Later sub-iterations (kit >= 2): warm-start from the previous
       ! sub-iteration's u,v,w,p as the initial guess; only the frozen
       ! explicit RHS (f, u_e) is restored inside restore_fluid.
       call this%restore_fluid(warm_start = (kit /= 1))

       ! Implicit kinematics
       call this%update_kinematics_bdfk(time, beta, nadv)

       do i = 1, this%nbodies_fsi
          this%batch_ids(i) = this%fsi_bodies(i)%ale_id
          this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
          this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
       end do
       call this%ale%update_mesh_velocity(this%c_Xh, time, &
            override_ids = this%batch_ids, &
            override_trans = this%batch_trans, &
            override_ang = this%batch_ang, &
            out_prescribed_vels = this%temp_prescribed_vels, mode = 0)

       ! Implicit mesh position update.
       call this%ale%advance_mesh_implicit(this%c_Xh, time, beta, nadv, &
            this%mesh_x_lag, this%mesh_y_lag, this%mesh_z_lag)

       ! Implicit fluid solve on the repositioned mesh. No history commit,
       ! RHS already assembled.
       call this%step_ext(time, dt_controller, greens_function = .false., &
            skip_ale_msh_vel_update = .true., skip_rhs_assembly = .true., &
            skip_ale_advance = .true., skip_projection = (kit /= 1))

       ! Structural system RHS.
       call this%calc_fsi_terms(time, beta, nadv)
       call this%add_spring_jacobian(time, beta)

       do i = 1, this%nbodies_fsi
          call this%fsi_bodies(i)%force_monitor%compute_(time)
          F_fluid(1:3) = this%fsi_bodies(i)%force_monitor%total_force
          F_fluid(4:6) = this%fsi_bodies(i)%force_monitor%total_torque
          do k = 1, 6
             row_g = this%fsi_dof_map(i, k)
             if (row_g > 0) this%B_global(row_g) = this%B_global(row_g) + &
                  F_fluid(k)
          end do
       end do

       ! Force residual, computed here while B_global still
       ! holds the interface imbalance.
       call subiter_dof_metric(this, .true., &
            fmetric, wf_body, wf_dof)

       ! Solve the 6-DOF structural system  M * dv = B
       if (this%total_active_dofs > 0) then
          allocate(M_linear(this%total_active_dofs, this%total_active_dofs))
          allocate(X_prev(this%total_active_dofs))
          M_linear = this%M_global
          max_nr = 1
          if (this%non_linear_correction_term) then
             max_nr = 20
             call neko_log%message("  --- Fixed-point Iteration ---")
          end if
          do nr_iter = 1, max_nr
             this%M_global = M_linear
             if (this%non_linear_correction_term) then
                call add_fsi_non_linear_matrices(this%nbodies_fsi, &
                     this%fsi_bodies, this%fsi_dof_map, this%M_global, &
                     this%X_sol, this%ale%body_rot_matrices)
             end if
             X_prev = this%X_sol
             call linsolve_dense(this%total_active_dofs, this%M_global, &
                  this%B_global, this%X_sol)
             if (this%non_linear_correction_term) then
                nr_res = maxval(abs(this%X_sol - X_prev))

                write(msg, '(A, I2, A, ES13.6)') "    Iter: ", &
                     nr_iter, " | Max Residual: ", nr_res
                if (nr_res < 1.0e-14_rp) exit
             end if
          end do
          deallocate(M_linear, X_prev)
       end if

       ! Velocity increment; the un-relaxed structural
       ! correction X_sol, normalised per DOF by its own velocity scale.
       call subiter_dof_metric(this, .false., &
            vmetric, wv_body, wv_dof)

       ! Per-body relaxation factor.
       ! Each body gets its own factor, computed
       ! from that body's own increment vector.
       do i = 1, this%nbodies_fsi
          if (this%accel_method == ACCEL_AITKEN) then
             ! this body's increment vector (inactive DOFs stay 0).
             r_body = 0.0_rp
             do k = 1, 6
                row_g = this%fsi_dof_map(i, k)
                if (row_g > 0) r_body(k) = this%X_sol(row_g)
             end do
             if (kit == 1) then
                this%accel_omega(i) = this%relax
             else
                dr_body = r_body - this%accel_r_prev(:, i)
                den = dot_product(dr_body, dr_body)
                if (den > 1.0e-30_rp) then
                   num = dot_product(this%accel_r_prev(:, i), dr_body)
                   this%accel_omega(i) = -this%accel_omega(i) * num / den
                end if
                this%accel_omega(i) = max(this%aitken_min, &
                     min(this%aitken_max, this%accel_omega(i)))
             end if
             this%accel_r_prev(:, i) = r_body
          else
             ! ACCEL_CONSTANT
             this%accel_omega(i) = this%relax
          end if
       end do

       do i = 1, this%nbodies_fsi
          do k = 1, 6
             row_g = this%fsi_dof_map(i, k)
             if (row_g > 0) then
                this%fsi_bodies(i)%body_vel(k) = &
                     this%fsi_bodies(i)%body_vel(k) + &
                     this%accel_omega(i) * this%X_sol(row_g)
             end if
          end do
       end do

       ! representative factor for the log (least-relaxed body).
       omega = maxval(this%accel_omega)

       ! Worst velocity/force metric with the body: DOF that
       ! owns it, plus the representative relaxation factor.
       if (wv_body > 0) then
          write(vtag, '(A,I0,A,A)') 'b', wv_body, ' ', dof_lbl(wv_dof)
       else
          vtag = 'none'
       end if
       if (wf_body > 0) then
          write(ftag, '(A,I0,A,A)') 'b', wf_body, ' ', dof_lbl(wf_dof)
       else
          ftag = 'none'
       end if

       t_sub = MPI_Wtime() - t_sub0

       write(msg, '(A,I4,A,ES11.4,A,A,A,ES11.4,A,A,A,F6.3,A,ES15.7,A)') &
            "  subiter ", kit, &
            " | delta_v=", vmetric, " @" , trim(vtag), &
            " | F_res=", fmetric, " @", trim(ftag), &
            " | omega=", omega, &
            " | substep_time=", t_sub, " s"
       call neko_log%message(trim(msg))

       ! Per-body-per-DOF breakdown (only the criterion's metric(s), only
       ! active DOFs). Evidence that every DOF is individually converged.
       do i = 1, this%nbodies_fsi
          if (show_v) then
             sub = ''
             do k = 1, 6
                if (this%rel_dof_v(k, i) >= 0.0_rp) then
                   write(tmp, '(1X,A,A,ES10.3)') dof_lbl(k), '=', &
                        this%rel_dof_v(k, i)
                   sub = trim(sub) // trim(tmp)
                end if
             end do
             if (len_trim(sub) > 0) then
                write(msg, '(A,I0,3A)') "    body ", i, " ", &
                     trim(dv_lbl), " :" // trim(sub)
                call neko_log%message(trim(msg))
             end if
          end if
          if (show_f) then
             sub = ''
             do k = 1, 6
                if (this%rel_dof_f(k, i) >= 0.0_rp) then
                   write(tmp, '(1X,A,A,ES10.3)') dof_lbl(k), '=', &
                        this%rel_dof_f(k, i)
                   sub = trim(sub) // trim(tmp)
                end if
             end do
             if (len_trim(sub) > 0) then
                write(msg, '(A,I0,3A)') "    body ", i, " ", &
                     trim(fr_lbl), " :" // trim(sub)
                call neko_log%message(trim(msg))
             end if
          end if
       end do

       ! Convergence.
       if (kit >= this%min_subiter) then
          converged = subiter_is_converged(this, vmetric, fmetric)
          if (converged) then
             write(msg, '(A,I0,A)') "  converged in ", kit, " sub-iterations."
             call neko_log%message(trim(msg))
             exit
          end if
       end if
    end do

    if (.not. converged) then
       write(msg, '(A,I0,A)') "  WARNING: FSI sub-iteration did not converge in ", &
            this%max_subiter, " iterations."
       call neko_log%message(trim(msg))
    end if

    ! Newmark/CN consolidation: the Aitken update applies the final velocity
    ! increment *after* the last kinematics/mesh-velocity build, so recompute
    ! the geometry and the mesh velocity on the final converged body_vel. This
    ! makes the committed mesh position, pivot, ghost trackers, disp_rel
    ! (logged) and wm^{n+1} all mutually consistent to machine precision, not
    ! merely to the coupling tolerance. (BDF is left exactly as before.)
    if (this%structure_newmark) then
       call this%update_kinematics_bdfk(time, beta, nadv)
       do i = 1, this%nbodies_fsi
          this%batch_ids(i) = this%fsi_bodies(i)%ale_id
          this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
          this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
       end do
       call this%ale%update_mesh_velocity(this%c_Xh, time, &
            override_ids = this%batch_ids, &
            override_trans = this%batch_trans, &
            override_ang = this%batch_ang, &
            out_prescribed_vels = this%temp_prescribed_vels, mode = 0)
       do i = 1, this%nbodies_fsi
          this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) = &
               this%temp_prescribed_vels(:, i)
       end do
    end if

    ! The last trial solution (u,v,w,p) and the last mesh
    ! x^{n+1}(v_s) are kept. We only advance the geometry / structure
    ! histories now.
    call this%commit_histories(time, beta, nadv)

    ! Final mesh velocity = predictor for the next step.
    do i = 1, this%nbodies_fsi
       this%batch_ids(i) = this%fsi_bodies(i)%ale_id
       this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
       this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
    end do
    call this%ale%update_mesh_velocity(this%c_Xh, time, &
         override_ids = this%batch_ids, &
         override_trans = this%batch_trans, &
         override_ang = this%batch_ang, &
         out_prescribed_vels = this%temp_prescribed_vels, mode = 0)
    do i = 1, this%nbodies_fsi
       this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) = &
            this%temp_prescribed_vels(:, i)
    end do

    call fsi_prep_checkpoint(this%nbodies_fsi, this%fsi_bodies, &
         this%global_disp_rel, this%global_body_vel, &
         this%global_body_vel_lag, this%global_moving_frame_presc_vel, &
         global_body_acc = this%global_body_acc, &
         global_frame_acc = this%global_frame_acc)

    call this%log_fsi_results(time)
    call this%ale%log_pivot(time)
    call this%ale%log_rot_angles(time)

  end subroutine fluid_subiter_step


  ! EXT-k predictor for the sub-iteration trial-1 velcoity-guess.
  subroutine subiter_apply_ext_predictor(this, time)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i, k, row_g, eff_order
    real(kind=rp) :: ec(4), dt_loc(10)

    if (this%predictor_order <= 0) then
       ! Reuse the fluid's advection EXT-k coefficients.
       ec = this%ext_bdf%advection_coeffs%x
    else
       ! Fixed order, capped by the warmup ramp; coefficients for current dt.
       eff_order = max(1, min(this%predictor_order, this%ext_bdf%nadv))
       dt_loc = time%dt
       call this%ext_bdf%ext%compute_coeffs(ec, dt_loc, eff_order)
    end if

    do i = 1, this%nbodies_fsi
       do k = 1, 6
          row_g = this%fsi_dof_map(i, k)
          if (row_g <= 0) cycle
          this%fsi_bodies(i)%body_vel(k) = &
               ec(1) * this%fsi_bodies(i)%body_vel_lag(k, 1) + &
               ec(2) * this%fsi_bodies(i)%body_vel_lag(k, 2) + &
               ec(3) * this%fsi_bodies(i)%body_vel_lag(k, 3)
       end do
    end do
  end subroutine subiter_apply_ext_predictor

  ! Per-body-per-DOF convergence metric.
  !   is_force = .false. : numerator = velocity increment X_sol(row),
  !                        own scale = |body_vel(dof)|
  !   is_force = .true.  : numerator = force residual  B_global(row),
  !                        own scale = |total_force/torque(dof)|
  ! Each DOF is normalised (when conv_relative) by its own scale, floored at
  ! REL_FLOOR times the peak scale among active DOFs of the same *type*
  ! (translation k<=3, rotation k>=4).
  ! The returned metric is the
  ! worst DOF, with its (body, dof) index.
  subroutine subiter_dof_metric(this, is_force, metric, wbody, wdof)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    logical, intent(in) :: is_force
    real(kind=rp), intent(out) :: metric
    integer, intent(out) :: wbody, wdof
    integer :: i, k, row_g
    real(kind=rp) :: st, sr, ref, own, numer, denom, r
    real(kind=rp), parameter :: TINY = 1.0e-30_rp

    st = 0.0_rp
    sr = 0.0_rp
    do i = 1, this%nbodies_fsi
       do k = 1, 6
          if (this%fsi_dof_map(i, k) <= 0) cycle
          if (is_force) then
             if (k <= 3) then
                own = abs(this%fsi_bodies(i)%force_monitor%total_force(k))
             else
                own = abs(this%fsi_bodies(i)%force_monitor%total_torque(k-3))
             end if
          else
             own = abs(this%fsi_bodies(i)%body_vel(k))
          end if
          if (k <= 3) then
             st = max(st, own)
          else
             sr = max(sr, own)
          end if
       end do
    end do

    metric = 0.0_rp
    wbody = 0
    wdof = 0
    do i = 1, this%nbodies_fsi
       do k = 1, 6
          row_g = this%fsi_dof_map(i, k)
          if (row_g .le. 0) then
             if (is_force) then
                this%rel_dof_f(k, i) = -1.0_rp
             else
                this%rel_dof_v(k, i) = -1.0_rp
             end if
             cycle
          end if

          if (is_force) then
             numer = abs(this%B_global(row_g))
             if (k .le. 3) then
                own = abs(this%fsi_bodies(i)%force_monitor%total_force(k))
             else
                own = abs(this%fsi_bodies(i)%force_monitor%total_torque(k-3))
             end if
          else
             numer = abs(this%X_sol(row_g))
             own = abs(this%fsi_bodies(i)%body_vel(k))
          end if

          if (this%conv_relative) then
             if (k .le. 3) then
                ref = st
             else
                ref = sr
             end if
             denom = max(own, REL_FLOOR * ref)
             if (denom .lt. TINY) denom = TINY
             r = numer / denom
          else
             r = numer
          end if

          if (is_force) then
             this%rel_dof_f(k, i) = r
          else
             this%rel_dof_v(k, i) = r
          end if
          if (r > metric) then
             metric = r
             wbody = i
             wdof = k
          end if
       end do
    end do
  end subroutine subiter_dof_metric

  ! vmetric / fmetric are the worst-DOF relative (or absolute) errors already
  ! formed by subiter_dof_metric; this routine only applies the criterion.
  function subiter_is_converged(this, vmetric, fmetric) result(ok)
    class(fluid_pnpn_fsi_subiter_t), intent(in) :: this
    real(kind=rp), intent(in) :: vmetric, fmetric
    logical :: ok, v_ok, f_ok

    v_ok = (vmetric < this%conv_tol)
    f_ok = (fmetric < this%conv_tol)

    select case (trim(this%conv_criterion))
    case ('velocity')
       ok = v_ok
    case ('force')
       ok = f_ok
    case ('both')
       ok = v_ok .and. f_ok
    case ('either')
       ok = v_ok .or. f_ok
    case default
       ok = v_ok
    end select
  end function subiter_is_converged

  ! Fluid reference (u^n, f, u_e) snapshot for the sub-iteration.
  subroutine subiter_snapshot_fluid(this)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    call field_copy(this%u_ref, this%u)
    call field_copy(this%v_ref, this%v)
    call field_copy(this%w_ref, this%w)
    call field_copy(this%p_ref, this%p)
    call field_copy(this%fx_ref, this%f_x)
    call field_copy(this%fy_ref, this%f_y)
    call field_copy(this%fz_ref, this%f_z)
    call field_copy(this%ue_ref, this%u_e)
    call field_copy(this%ve_ref, this%v_e)
    call field_copy(this%we_ref, this%w_e)
  end subroutine subiter_snapshot_fluid

  subroutine subiter_restore_fluid(this, warm_start)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    ! warm_start = .true.  -> keep the current u,v,w,p (the previous
    !   sub-iteration's solution) as the initial guess for this trial's
    !   linear solve; only the frozen explicit RHS (f, u_e) is restored.
    ! warm_start = .false. (default) -> also restore u,v,w,p to u^n
    !   (used on the first sub-iteration).
    logical, intent(in), optional :: warm_start
    logical :: ws
    ! ===============================================================
    ws = .false.
    if (present(warm_start)) ws = warm_start

    ! Explicit RHS inputs: always restored.
    call field_copy(this%f_x, this%fx_ref)
    call field_copy(this%f_y, this%fy_ref)
    call field_copy(this%f_z, this%fz_ref)
    call field_copy(this%u_e, this%ue_ref)
    call field_copy(this%v_e, this%ve_ref)
    call field_copy(this%w_e, this%we_ref)

    if (.not. ws) then
       call field_copy(this%u, this%u_ref)
       call field_copy(this%v, this%v_ref)
       call field_copy(this%w, this%w_ref)
       call field_copy(this%p, this%p_ref)
    end if
  end subroutine subiter_restore_fluid

  !  Implicit (BDF-k) kinematics from the FROZEN n-histories, using v_s.
  !  Writes ale_pivot%pos, trackers%pos and disp_rel, and recomputes R^{n+1}.
  subroutine subiter_update_kinematics(this, time, beta, nadv)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    integer, intent(in) :: nadv
    integer :: i, a
    real(kind=rp) :: pvel(3), omega_tot(3), pvel_prev(3), omega_tot_prev(3)

    do i = 1, this%nbodies_fsi
       a = this%fsi_bodies(i)%ale_id

       ! total (frame + relative) rigid-body velocities of this body.
       pvel = this%fsi_bodies(i)%moving_frame_presc_vel(1:3, 0) + &
            this%fsi_bodies(i)%body_vel(1:3)
       omega_tot = this%fsi_bodies(i)%moving_frame_presc_vel(4:6, 0) + &
            this%fsi_bodies(i)%body_vel(4:6)

       if (this%structure_newmark) then
          ! Previous-step (n) total velocities for the CN (trapezoidal) rule:
          ! frame prescribed velocity at n = moving_frame_presc_vel(:,1),
          ! relative velocity at n = body_vel_lag(:,1).
          pvel_prev = this%fsi_bodies(i)%moving_frame_presc_vel(1:3, 1) + &
               this%fsi_bodies(i)%body_vel_lag(1:3, 1)
          omega_tot_prev = this%fsi_bodies(i)%moving_frame_presc_vel(4:6, 1) + &
               this%fsi_bodies(i)%body_vel_lag(4:6, 1)

          ! ALE geometry: pivot, rotation trackers, and rotation matrix (CN).
          call this%ale%step_body_kinematics_implicit(a, time, beta, nadv, &
               pvel, omega_tot, this%pivot_hist(:, :, i), &
               this%ghost_hist(:, :, :, i), &
               pvel_prev = pvel_prev, omega_tot_prev = omega_tot_prev)

          ! relative displacement (CN): disp_rel = disp_rel^n
          !   + (dt/2)(body_vel + body_vel^n).
          call this%ale%scheme%integrate_6dof(this%fsi_bodies(i)%disp_rel, &
               this%fsi_bodies(i)%body_vel, time, nadv, &
               hist6 = this%disp_hist(:, :, i), &
               v6_prev = this%fsi_bodies(i)%body_vel_lag(:, 1))
       else
          ! ALE geometry: pivot, rotation trackers, and rotation matrix (BDF)
          ! from the frozen (FSI-owned) pivot/ghost histories.
          call this%ale%step_body_kinematics_implicit(a, time, beta, nadv, &
               pvel, omega_tot, this%pivot_hist(:, :, i), &
               this%ghost_hist(:, :, :, i))

          ! relative displacement (BDF): d(disp_rel)/dt = body_vel.
          call this%ale%scheme%integrate_6dof(this%fsi_bodies(i)%disp_rel, &
               this%fsi_bodies(i)%body_vel, time, nadv, &
               beta = beta, hist6 = this%disp_hist(:, :, i))
       end if
    end do

  end subroutine subiter_update_kinematics

  ! Advance every geometry / structure history once, using
  ! the converged velocity field.
  subroutine subiter_commit_histories(this, time, beta, nadv)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    integer, intent(in) :: nadv
    integer :: i, j, g, a, h

    ! Rewind geometry to x^n (= lag 1) so that B = mass(x^n), shift the
    ! B-history, then restore the converged x^{n+1} and recompute the metrics.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%c_Xh%dof%x_d, this%mesh_x_lag(1)%x_d, this%n_mesh)
       call device_copy(this%c_Xh%dof%y_d, this%mesh_y_lag(1)%x_d, this%n_mesh)
       call device_copy(this%c_Xh%dof%z_d, this%mesh_z_lag(1)%x_d, this%n_mesh)
    else
       call copy(this%c_Xh%dof%x, this%mesh_x_lag(1)%x, this%n_mesh)
       call copy(this%c_Xh%dof%y, this%mesh_y_lag(1)%x, this%n_mesh)
       call copy(this%c_Xh%dof%z, this%mesh_z_lag(1)%x, this%n_mesh)
    end if
    call this%c_Xh%recompute_metrics()
    call this%c_Xh%update_B_history()

    ! reposition to the converged x^{n+1} and recompute.
    call this%ale%advance_mesh_implicit(this%c_Xh, time, beta, nadv, &
         this%mesh_x_lag, this%mesh_y_lag, this%mesh_z_lag)
    call this%c_Xh%recompute_metrics()
    call this%adv%recompute_metrics(this%c_Xh, .true.)

    do j = this%n_lag, 2, -1
       this%mesh_x_lag(j) = this%mesh_x_lag(j - 1)
       this%mesh_y_lag(j) = this%mesh_y_lag(j - 1)
       this%mesh_z_lag(j) = this%mesh_z_lag(j - 1)
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%mesh_x_lag(1)%x_d, this%c_Xh%dof%x_d, this%n_mesh)
       call device_copy(this%mesh_y_lag(1)%x_d, this%c_Xh%dof%y_d, this%n_mesh)
       call device_copy(this%mesh_z_lag(1)%x_d, this%c_Xh%dof%z_d, this%n_mesh)
    else
       call copy(this%mesh_x_lag(1)%x, this%c_Xh%dof%x, this%n_mesh)
       call copy(this%mesh_y_lag(1)%x, this%c_Xh%dof%y, this%n_mesh)
       call copy(this%mesh_z_lag(1)%x, this%c_Xh%dof%z, this%n_mesh)
    end if

    ! Shift the rigid-body position histories
    do i = 1, this%nbodies_fsi
       a = this%fsi_bodies(i)%ale_id
       do j = this%n_lag, 2, -1
          this%pivot_hist(:, j, i) = this%pivot_hist(:, j - 1, i)
          this%disp_hist(:, j, i)  = this%disp_hist(:, j - 1, i)
          do g = 1, 2
             this%ghost_hist(:, j, g, i) = this%ghost_hist(:, j - 1, g, i)
          end do
       end do
       this%pivot_hist(:, 1, i) = this%ale%ale_pivot(a)%pos
       this%disp_hist(:, 1, i)  = this%fsi_bodies(i)%disp_rel
       do g = 1, 2
          h = this%ale%ghost_handles(g, a)
          this%ghost_hist(:, 1, g, i) = this%ale%trackers(h)%pos
       end do
    end do

    ! Newmark: update the stored acceleration
    !   a^{n+1} = (2/dt)(v^{n+1} - v^n) - a^n,
    ! using v^n = body_vel_lag(:,1) and a^n = body_acc, BEFORE the velocity
    ! history is shifted below (which would overwrite v^n).
    if (this%structure_newmark) then
       do i = 1, this%nbodies_fsi
          this%fsi_bodies(i)%body_acc = &
               (2.0_rp / time%dt) * (this%fsi_bodies(i)%body_vel - &
               this%fsi_bodies(i)%body_vel_lag(:, 1)) &
               - this%fsi_bodies(i)%body_acc
       end do
    end if

    ! Shift the body velocity history
    do i = 1, this%nbodies_fsi
       do j = nadv, 2, -1
          this%fsi_bodies(i)%body_vel_lag(:, j) = &
               this%fsi_bodies(i)%body_vel_lag(:, j - 1)
       end do
       this%fsi_bodies(i)%body_vel_lag(:, 1) = this%fsi_bodies(i)%body_vel
    end do

  end subroutine subiter_commit_histories

  !  Prescribed frame motion.
  subroutine subiter_query_frame_prescribed_motion(this, time, beta, nadv)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    integer, intent(in) :: nadv
    integer :: i, k

    ! shift the prescribed-frame velocity history.
    do i = 1, this%nbodies_fsi
       do k = nadv, 1, -1
          this%fsi_bodies(i)%moving_frame_presc_vel(:, k) = &
               this%fsi_bodies(i)%moving_frame_presc_vel(:, k - 1)
       end do
       this%batch_ids(i) = this%fsi_bodies(i)%ale_id
    end do

    ! prescribed motion only (mode = 2 does not touch wm).
    call this%ale%update_mesh_velocity(this%c_Xh, time, &
         override_ids = this%batch_ids, &
         out_prescribed_vels = this%temp_prescribed_vels, mode = 2)

    do i = 1, this%nbodies_fsi
       this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) = &
            this%temp_prescribed_vels(:, i)
    end do

    ! prescribed-frame acceleration a_f from the prescribed velocity history,
    ! discretised consistently with the structure integrator.
    if (this%structure_newmark) then
       ! CN / Newmark (trapezoidal) recurrence, matching the structure:
       !   a_f^{n+1} = (2/dt)(v_f^{n+1} - v_f^n) - a_f^n.
       ! v_f^{n+1} = presc_vel(:,0), v_f^n = presc_vel(:,1), a_f^n =
       ! moving_frame_presc_acc_prev (0 at start / on restart). For the common
       ! inertial-frame case (presc_vel = 0) this is identically zero.
       do i = 1, this%nbodies_fsi
          this%fsi_bodies(i)%moving_frame_presc_acc = &
               (2.0_rp / time%dt) * &
               (this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) - &
               this%fsi_bodies(i)%moving_frame_presc_vel(:, 1)) - &
               this%fsi_bodies(i)%moving_frame_presc_acc_prev
          this%fsi_bodies(i)%moving_frame_presc_acc_prev = &
               this%fsi_bodies(i)%moving_frame_presc_acc
       end do
    else
       ! BDF-k of the prescribed velocity (unchanged).
       do i = 1, this%nbodies_fsi
          this%fsi_bodies(i)%moving_frame_presc_acc = 0.0_rp
          do k = 0, nadv
             this%fsi_bodies(i)%moving_frame_presc_acc = &
                  this%fsi_bodies(i)%moving_frame_presc_acc + &
                  (beta(k) * this%fsi_bodies(i)%moving_frame_presc_vel(:, k)) / &
                  time%dt
          end do
       end do
    end if
  end subroutine subiter_query_frame_prescribed_motion

  !  Assemble the structural inertial system (M_global, B_global).
  subroutine subiter_calc_fsi_terms(this, time, beta, nadv)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    integer, intent(in) :: nadv
    real(kind=rp) :: gamma, two_over_dt
    real(kind=rp), allocatable :: accel_hist(:,:)
    integer :: i

    if (this%structure_newmark) then
       ! Newmark (constant average acceleration): a^{n+1} = gamma*v_guess
       ! + a_hist, with gamma = 2/dt and a_hist = -g^n = -(2/dt)v^n - a^n,
       ! v^n = body_vel_lag(:,1), a^n = body_acc.
       two_over_dt = 2.0_rp / time%dt
       gamma = two_over_dt
       allocate(accel_hist(6, this%nbodies_fsi))
       do i = 1, this%nbodies_fsi
          accel_hist(:, i) = - two_over_dt * &
               this%fsi_bodies(i)%body_vel_lag(:, 1) &
               - this%fsi_bodies(i)%body_acc
       end do
       call assemble_structural_inertial_terms(this%nbodies_fsi, &
            this%fsi_bodies, this%fsi_dof_map, this%M_global, this%B_global, &
            this%ale%body_rot_matrices, time, gamma, beta, nadv, &
            this%gravity_vec, accel_hist = accel_hist)
       deallocate(accel_hist)
    else
       gamma = beta(0) / time%dt
       call assemble_structural_inertial_terms(this%nbodies_fsi, &
            this%fsi_bodies, this%fsi_dof_map, this%M_global, this%B_global, &
            this%ale%body_rot_matrices, time, gamma, beta, nadv, &
            this%gravity_vec)
    end if
  end subroutine subiter_calc_fsi_terms

  ! Newmark 'balance' initial acceleration:  a^0 = M_mass^{-1} F_net^0.
  !
  ! F_net^0 = (external + fluid forces) - (frame + velocity-dependent inertial
  ! terms) is obtained by reusing assemble_structural_inertial_terms with the
  ! relative acceleration forced to zero (accel_hist = -gamma*body_vel makes
  ! a_full = gamma*v_guess + accel_hist = 0), so its residual B equals exactly
  ! M_mass*a^0. M_mass is the pure generalized rigid-body mass matrix
  !   [  m I        -m [c]x ]
  !   [  m [c]x       I_P   ]
  ! built here from the current geometry (c = R*offset_com, I_P = R I_body R^T).
  ! M_global/B_global are overwritten but are reassembled by the sub-iteration
  ! loop afterwards.
  subroutine subiter_newmark_init_accel(this, time, beta, nadv)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    integer, intent(in) :: nadv
    real(kind=rp), allocatable :: accel_hist(:,:), Mmass(:,:), a0(:)
    real(kind=rp) :: gamma, F_fluid(6), R(3,3), c(3), I_body(3,3), I_P(3,3)
    real(kind=rp) :: Cs(3,3), M_loc(6,6), I3(3,3), m
    integer :: i, j, k, row_g, col_g, a

    if (this%total_active_dofs <= 0) return

    gamma = 2.0_rp / time%dt

    ! F_net^0 into B_global (relative acceleration forced to zero).
    allocate(accel_hist(6, this%nbodies_fsi))
    do i = 1, this%nbodies_fsi
       accel_hist(:, i) = - gamma * this%fsi_bodies(i)%body_vel
    end do
    call assemble_structural_inertial_terms(this%nbodies_fsi, &
         this%fsi_bodies, this%fsi_dof_map, this%M_global, this%B_global, &
         this%ale%body_rot_matrices, time, gamma, beta, nadv, &
         this%gravity_vec, accel_hist = accel_hist)
    deallocate(accel_hist)

    ! Add the initial fluid force/torque to the net force.
    do i = 1, this%nbodies_fsi
       call this%fsi_bodies(i)%force_monitor%compute_(time)
       F_fluid(1:3) = this%fsi_bodies(i)%force_monitor%total_force
       F_fluid(4:6) = this%fsi_bodies(i)%force_monitor%total_torque
       do k = 1, 6
          row_g = this%fsi_dof_map(i, k)
          if (row_g > 0) this%B_global(row_g) = this%B_global(row_g) + F_fluid(k)
       end do
    end do

    ! Build the generalized mass matrix M_mass.
    I3 = 0.0_rp
    I3(1,1) = 1.0_rp; I3(2,2) = 1.0_rp; I3(3,3) = 1.0_rp
    allocate(Mmass(this%total_active_dofs, this%total_active_dofs))
    allocate(a0(this%total_active_dofs))
    Mmass = 0.0_rp
    do i = 1, this%nbodies_fsi
       a = this%fsi_bodies(i)%ale_id
       m = this%fsi_bodies(i)%mass
       R = this%ale%body_rot_matrices(:, :, a)
       I_body = this%fsi_bodies(i)%I_body_tensor
       c = matmul(R, this%fsi_bodies(i)%local_offset_com)
       I_P = matmul(R, matmul(I_body, transpose(R)))
       ! Cs = skew(c)
       Cs = 0.0_rp
       Cs(1,2) = -c(3); Cs(1,3) =  c(2)
       Cs(2,1) =  c(3); Cs(2,3) = -c(1)
       Cs(3,1) = -c(2); Cs(3,2) =  c(1)
       M_loc = 0.0_rp
       M_loc(1:3, 1:3) = m * I3
       M_loc(1:3, 4:6) = - m * Cs
       M_loc(4:6, 1:3) =   m * Cs
       M_loc(4:6, 4:6) = I_P
       do j = 1, 6
          row_g = this%fsi_dof_map(i, j)
          if (row_g > 0) then
             do k = 1, 6
                col_g = this%fsi_dof_map(i, k)
                if (col_g > 0) Mmass(row_g, col_g) = Mmass(row_g, col_g) + &
                     M_loc(j, k)
             end do
          end if
       end do
    end do

    ! Solve M_mass a^0 = F_net^0 and scatter into body_acc.
    call linsolve_dense(this%total_active_dofs, Mmass, this%B_global, a0)
    do i = 1, this%nbodies_fsi
       this%fsi_bodies(i)%body_acc = 0.0_rp
       do k = 1, 6
          row_g = this%fsi_dof_map(i, k)
          if (row_g > 0) this%fsi_bodies(i)%body_acc(k) = a0(row_g)
       end do
    end do
    deallocate(Mmass, a0)

    call neko_log%message("  Newmark: initial acceleration a^0 set from " // &
         "the initial force balance.")
  end subroutine subiter_newmark_init_accel

  ! Add the spring contribution to the LHS Jacobian M_global.
  !
  ! calc_fsi_terms (via assemble_structural_inertial_terms, shared with the
  ! Green's scheme) places the spring only in the residual B, as
  ! B -= K*(disp_rel - pos_eq).
  !
  ! That is correct for the Green's scheme, where disp_rel is advanced
  ! explicitly (Adams-Bashforth) and is therefore independent of the velocity
  ! being solved for, so d(spring)/dv = 0.
  !
  ! In the sub-iteration, disp_rel is implicit in the current velocity iterate,
  ! with a scheme-dependent sensitivity dvdisp = d(disp_rel)/d(body_vel):
  !   BDF: disp_rel = (dt*v - sum_hist)/beta(0)   => dvdisp = dt/beta(0)
  !   CN : disp_rel = disp_rel^n + (dt/2)(v+v^n)  => dvdisp = dt/2
  ! K_lin/K_ang are diagonal per DOF and disp_rel(j) depends only on
  ! body_vel(j), so the Jacobian is diagonal: M(row,row) += K*dvdisp.
  subroutine subiter_add_spring_jacobian(this, time, beta)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    integer :: i, j, row_g
    real(kind=rp) :: dvdisp, k_dof

    if (this%total_active_dofs <= 0) return
    if (this%structure_newmark) then
       dvdisp = 0.5_rp * time%dt
    else
       dvdisp = time%dt / beta(0)
    end if

    do i = 1, this%nbodies_fsi
       do j = 1, 6
          row_g = this%fsi_dof_map(i, j)
          if (row_g > 0) then
             if (j <= 3) then
                k_dof = this%fsi_bodies(i)%K_lin(j) ! translational spring
             else
                k_dof = this%fsi_bodies(i)%K_ang(j - 3) ! torsional spring
             end if
             this%M_global(row_g, row_g) = this%M_global(row_g, row_g) &
                  + k_dof * dvdisp
          end if
       end do
    end do
  end subroutine subiter_add_spring_jacobian

  !  Logging
  subroutine subiter_log_results(this, time)
    class(fluid_pnpn_fsi_subiter_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    character(len=1024) :: msg
    character(len=128) :: fmt_res
    integer :: i

    if (this%nbodies_fsi == 0) return
    call neko_log%message("--------- FSI Results ----------")
    if (this%res_long_print) then
       fmt_res = '(A, I0, A, ES23.15, A, A, A, 3(ES22.15, :, 2X))'
    else
       fmt_res = '(A, I0, A, ES17.10, A, A, A, 3(ES17.10, :, 2X))'
    end if
    call neko_log%message("variable, step, time, body, x, y, z")
    do i = 1, this%nbodies_fsi
       write(msg, fmt_res) "FSI_DISP_L  ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", this%fsi_bodies(i)%disp_rel(1:3)
       call neko_log%message(trim(msg))
       write(msg, fmt_res) "FSI_DISP_A  ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", this%fsi_bodies(i)%disp_rel(4:6)
       call neko_log%message(trim(msg))
       write(msg, fmt_res) "FSI_VEL_L   ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", this%fsi_bodies(i)%body_vel(1:3)
       call neko_log%message(trim(msg))
       write(msg, fmt_res) "FSI_VEL_A   ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", this%fsi_bodies(i)%body_vel(4:6)
       call neko_log%message(trim(msg))
    end do
    call neko_log%message(" ")
  end subroutine subiter_log_results

  !  Restart
  subroutine fluid_subiter_restart(this, chkp)
    class(fluid_pnpn_fsi_subiter_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    type(time_state_t) :: t_restart
    integer :: i

    ! Fluid_pnpn scheme restart (fluid, ALE, adv, etc.)
    call this%fluid_pnpn_t%restart(chkp)

    if (this%if_fsi .and. this%nbodies_fsi > 0) then
       call fsi_restart_restore(this%nbodies_fsi, this%fsi_bodies, &
            this%global_disp_rel, this%global_body_vel, &
            this%global_body_vel_lag, this%global_moving_frame_presc_vel, &
            global_body_acc = this%global_body_acc, &
            global_frame_acc = this%global_frame_acc)

       t_restart%t = chkp%t
       t_restart%tstep = 0
       t_restart%dt = chkp%dtlag(1)

       do i = 1, this%nbodies_fsi
          this%batch_ids(i) = this%fsi_bodies(i)%ale_id
          this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
          this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
       end do
       call this%ale%update_mesh_velocity(this%c_Xh, t_restart, &
            override_ids = this%batch_ids, &
            override_trans = this%batch_trans, &
            override_ang = this%batch_ang, mode = 0)

       if (this%structure_newmark) then
          ! Seed the frozen wm^n for the CN reposition from the rebuilt mesh
          ! velocity above (the first post-restart step re-snapshots it too).
          call this%ale%snapshot_mesh_velocity()
          ! body_acc was restored from the checkpoint (or left at zero for an
          ! older checkpoint), so it is the correct a^n. Skip the first-step
          ! a^0 (re)initialisation.
          this%newmark_a0_done = .true.
          ! moving_frame_presc_acc_prev (a_f^n) was restored from the
          ! checkpoint above (bit 512). It is not reconstructable from the
          ! saved velocity lags, so it is stored explicitly to keep the
          ! restart bit-exact for an accelerating prescribed frame. For an
          ! older checkpoint without the field it stays at zero, matching the
          ! previous behaviour.
       end if

       if (chkp%fsi_subiter_restored) then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             do i = 1, this%n_lag
                call this%mesh_x_lag(i)%copy_from(HOST_TO_DEVICE, sync = .false.)
                call this%mesh_y_lag(i)%copy_from(HOST_TO_DEVICE, sync = .false.)
                call this%mesh_z_lag(i)%copy_from(HOST_TO_DEVICE, sync = .false.)
             end do
             call device_sync()
          end if
       end if
    end if
  end subroutine fluid_subiter_restart

  !  Destruction
  subroutine fluid_subiter_free(this)
    class(fluid_pnpn_fsi_subiter_t), intent(inout) :: this
    integer :: i, j, k

    if (allocated(this%fsi_bodies)) then
       do i = 1, this%nbodies_fsi
          call this%fsi_bodies(i)%force_monitor%free()
       end do
       deallocate(this%fsi_bodies)
    end if

    if (allocated(this%M_global)) deallocate(this%M_global)
    if (allocated(this%B_global)) deallocate(this%B_global)
    if (allocated(this%X_sol)) deallocate(this%X_sol)
    if (allocated(this%fsi_dof_map)) deallocate(this%fsi_dof_map)

    if (allocated(this%accel_omega)) deallocate(this%accel_omega)
    if (allocated(this%accel_r_prev)) deallocate(this%accel_r_prev)

    if (allocated(this%rel_dof_v)) deallocate(this%rel_dof_v)
    if (allocated(this%rel_dof_f)) deallocate(this%rel_dof_f)

    if (allocated(this%u_g)) then
       do k = 1, size(this%u_g)
          call this%u_g(k)%free()
          call this%v_g(k)%free()
          call this%w_g(k)%free()
          call this%p_g(k)%free()
       end do
       deallocate(this%u_g, this%v_g, this%w_g, this%p_g)
    end if
    if (allocated(this%proj_prs_green)) then
       do k = 1, size(this%proj_prs_green)
          call this%proj_prs_green(k)%free()
          call this%proj_vel_green(k)%free()
       end do
       deallocate(this%proj_prs_green, this%proj_vel_green)
    end if

    if (allocated(this%batch_ids)) deallocate(this%batch_ids)
    if (allocated(this%batch_trans)) deallocate(this%batch_trans)
    if (allocated(this%batch_ang)) deallocate(this%batch_ang)
    if (allocated(this%temp_prescribed_vels)) deallocate(this%temp_prescribed_vels)

    if (allocated(this%mesh_x_lag)) then
       do j = 1, size(this%mesh_x_lag)
          call this%mesh_x_lag(j)%free()
          call this%mesh_y_lag(j)%free()
          call this%mesh_z_lag(j)%free()
       end do
       deallocate(this%mesh_x_lag, this%mesh_y_lag, this%mesh_z_lag)
    end if
    if (allocated(this%pivot_hist)) deallocate(this%pivot_hist)
    if (allocated(this%ghost_hist)) deallocate(this%ghost_hist)
    if (allocated(this%disp_hist)) deallocate(this%disp_hist)

    call this%u_ref%free()
    call this%v_ref%free()
    call this%w_ref%free()
    call this%p_ref%free()
    call this%fx_ref%free()
    call this%fy_ref%free()
    call this%fz_ref%free()
    call this%ue_ref%free()
    call this%ve_ref%free()
    call this%we_ref%free()

    call this%fluid_pnpn_t%free()
  end subroutine fluid_subiter_free

end module fluid_pnpn_fsi_subiteration