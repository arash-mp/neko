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
!> PnPn fluid scheme coupled to rigid bodies with the Green's function method.
module fluid_pnpn_fsi_greens
  use fsi_dynamics, only : fsi_body_t, assemble_structural_inertial_terms
  use fsi_manager, only : fsi_manager_init, fsi_solve_structure, &
       fsi_register_checkpoint, fsi_prep_checkpoint, fsi_restart_restore
  use fluid_pnpn, only : fluid_pnpn_t
  use field, only : field_t
  use field_math, only : field_copy, field_add2s2, field_cfill
  use num_types, only : rp, dp
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use projection, only : projection_t
  use projection_vel, only : projection_vel_t
  use bc, only : bc_t, BC_DIRICHLET
  use bc_list, only : bc_list_t
  use zero_dirichlet, only : zero_dirichlet_t
  use no_slip, only : no_slip_t
  use inflow, only : inflow_t
  use field_dirichlet_vector, only : field_dirichlet_vector_t
  use wall_model_bc, only : wall_model_bc_t
  use symmetry_aligned, only : symmetry_aligned_t
  use symmetry, only : symmetry_t
  use non_normal, only : non_normal_t
  use shear_stress, only : shear_stress_t
  use gs_ops, only : GS_OP_MIN, GS_OP_MAX
  use operators, only : rotate_cyc
  use device, only : device_event_sync, glb_cmd_event
  use profiler, only : profiler_start_region, profiler_end_region
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  use utils, only : neko_error
  use logger, only : neko_log
  use mesh, only : mesh_t
  use user_intf, only : user_t, user_fsi_structural_parameters_intf, &
       user_fsi_structural_terms_intf, dummy_fsi_structural_terms
  use checkpoint, only : chkp_t
  use mpi_f08, only : MPI_Wtime
  use ab_time_scheme, only : ab_time_scheme_t
  use math, only : rzero
  use ale_manager, only : ALE_VEL_SUPERPOSE, ALE_VEL_EXCLUSIVE, ALE_VEL_QUERY
  implicit none
  private

  !> PnPn fluid scheme coupled to rigid bodies with the Green's function
  !! method of Fischer, P., Schmitt, M., & Tomboulides, A. (2017). Recent
  !! developments in spectral element simulations of moving-domain problems.
  !! Recent progress and modern challenges in applied mathematics, modeling
  !! and computational science, 213-244.
  type, public, extends(fluid_pnpn_t) :: fluid_pnpn_fsi_greens_t
     logical :: if_fsi = .false.
     logical :: skip_greens_solve = .false.
     logical :: non_linear_correction_term = .false.
     logical :: res_long_print = .false.
     real(kind=rp) :: gravity_vec(3) = 0.0_rp

     ! Storage for the Standard solution (u_s)
     type(field_t) :: u_s, v_s, w_s, p_s

     ! Storage for Green's function fields
     type(field_t), allocatable :: u_g(:)
     type(field_t), allocatable :: v_g(:)
     type(field_t), allocatable :: w_g(:)
     type(field_t), allocatable :: p_g(:)

     ! N-Body Storage
     integer :: nbodies_fsi = 0
     type(fsi_body_t), allocatable :: fsi_bodies(:)

     ! Mapping
     integer, allocatable :: fsi_dof_map(:,:)
     integer :: total_active_dofs = 0

     ! FSI global system matrices
     real(kind=rp), allocatable :: M_global(:,:)
     real(kind=rp), allocatable :: B_global(:)
     real(kind=rp), allocatable :: X_sol(:)
     type(projection_t), allocatable :: proj_prs_green(:)
     type(projection_vel_t), allocatable :: proj_vel_green(:)
     !> Boundary conditions of the Green's function problems.
     type(bc_list_t) :: bcs_vel_green
     type(bc_list_t) :: bcs_prs_green
     !> True while a Green's function problem is being solved.
     logical :: greens_mode = .false.
     ! Batch Arrays for ALE Override
     integer, allocatable :: batch_ids(:)
     real(kind=rp), allocatable :: batch_trans(:,:)
     real(kind=rp), allocatable :: batch_ang(:,:)
     real(kind=rp), allocatable :: temp_prescribed_vels(:,:)

     real(kind=rp), allocatable :: global_disp_rel(:)
     real(kind=rp), allocatable :: global_body_vel(:)
     real(kind=rp), allocatable :: global_body_vel_lag(:,:)
     real(kind=rp), allocatable :: global_moving_frame_presc_vel(:,:)

     !> User hook: runtime modification of FSI body parameters. Always
     !> associated (dummy if not registered in the user file).
     procedure(user_fsi_structural_parameters_intf), nopass, pointer :: &
          user_fsi_body_params => null()
     !> User hook: extra structural equation terms.
     procedure(user_fsi_structural_terms_intf), nopass, pointer :: &
          user_fsi_structural_terms => null()
     !> True when the user registered structural terms.
     logical :: has_user_structural_terms = .false.
     !> Accumulated wall time of the standard and of the Green's solves.
     real(kind=dp) :: total_elapsed_s = 0.0_dp
     real(kind=dp) :: total_elapsed_g = 0.0_dp
     !> Maximum number of pressure iterations in a Green's function solve.
     integer :: greens_prs_max_iter
   contains
     procedure, pass(this) :: init => fluid_fsi_init
     procedure, pass(this) :: step => fluid_fsi_step
     procedure, pass(this) :: free => fluid_fsi_free
     !> Restart from a previous solution.
     procedure, pass(this) :: restart => fluid_fsi_restart
     procedure, pass(this) :: bc_apply_vel => fluid_fsi_bc_apply_vel
     procedure, pass(this) :: bc_apply_prs => fluid_fsi_bc_apply_prs
     procedure, pass(this) :: calc_fsi_terms => &
          assemble_fsi_structural_inertial_terms
     procedure, pass(this) :: log_fsi_results => fluid_fsi_log_results
     procedure, pass(this) :: query_frame_prescribed_motion => &
          fluid_fsi_query_frame_prescribed_motion
  end type fluid_pnpn_fsi_greens_t

contains

  !> Initialise the fluid scheme and the FSI bodies.
  subroutine fluid_fsi_init(this, msh, lx, params, user, chkp)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    type(time_state_t) :: t_init
    integer :: i

    ! Initialize the base PnPn solver
    call this%fluid_pnpn_t%init(msh, lx, params, user, chkp)

    if (this%freeze) then
       call neko_error("FSI does not support case.fluid.freeze.")
    end if

    call fluid_fsi_setup_green_bcs(this)

    ! User hook for FSI body parameters (never null after user%init).
    this%user_fsi_body_params => user%fsi_structural_parameters

    ! User hook for extra structural terms (never null after user%init).
    this%user_fsi_structural_terms => user%fsi_structural_terms
    this%has_user_structural_terms = .not. associated( &
         this%user_fsi_structural_terms, dummy_fsi_structural_terms)

    ! Initialize Standard Fields locally
    call this%u_s%init(this%dm_Xh, 'u_s')
    call this%v_s%init(this%dm_Xh, 'v_s')
    call this%w_s%init(this%dm_Xh, 'w_s')
    call this%p_s%init(this%dm_Xh, 'p_s')

    ! Init fsi_manager
    call fsi_manager_init(params, this%ale, this%c_Xh, this%dm_Xh, &
         this%if_fsi, this%nbodies_fsi, this%fsi_bodies, this%fsi_dof_map, &
         this%total_active_dofs, &
         this%M_global, this%B_global, this%X_sol, &
         this%u_g, this%v_g, this%w_g, this%p_g, &
         this%res_long_print, this%gravity_vec, this%proj_prs_green, &
         this%proj_vel_green, this%global_disp_rel, &
         this%global_body_vel, this%global_body_vel_lag, &
         this%global_moving_frame_presc_vel, this%skip_greens_solve, &
         this%non_linear_correction_term)

    ! The main pressure solver's limit unless set
    call json_get_or_default(params, &
         'case.fluid.fsi.pressure_solver.max_iterations', &
         this%greens_prs_max_iter, this%ksp_prs%max_iter)
    if (this%greens_prs_max_iter .lt. 1) then
       call neko_error("case.fluid.fsi.pressure_solver.max_iterations " // &
            "must be at least 1")
    end if

    call fsi_register_checkpoint(this%chkp, this%global_disp_rel, &
         this%global_body_vel, this%global_body_vel_lag, &
         this%global_moving_frame_presc_vel)

    if (this%nbodies_fsi > 0) then
       allocate(this%batch_ids(this%nbodies_fsi))
       allocate(this%batch_trans(3, this%nbodies_fsi))
       allocate(this%batch_ang(3, this%nbodies_fsi))

       allocate(this%temp_prescribed_vels(6, this%nbodies_fsi))
       this%temp_prescribed_vels = 0.0_rp

    end if

    ! For FSI, we calculate the inital mesh velocity here.
    ! In case of restart, we skip this.
    if (this%nbodies_fsi > 0 .and. &
         (.not. params%valid_path('case.restart_file'))) then
       t_init%t = 0.0_rp
       t_init%tstep = 0
       t_init%dt = 0.0_rp
       do i = 1, this%nbodies_fsi
          this%batch_ids(i) = this%fsi_bodies(i)%ale_id
          this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
          this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
       end do
       ! Apply Initial Guess + Prescribed Motion
       call this%ale%update_mesh_velocity(this%c_Xh, t_init, &
            override_ids = this%batch_ids, &
            override_trans = this%batch_trans, &
            override_ang = this%batch_ang, &
            out_prescribed_vels = this%temp_prescribed_vels, &
            mode = ALE_VEL_SUPERPOSE)
       do i = 1, this%nbodies_fsi
          this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) = &
               this%temp_prescribed_vels(:, i)
       end do
    end if

  end subroutine fluid_fsi_init

  !> Advance the fluid and the FSI bodies one time step.
  subroutine fluid_fsi_step(this, time, dt_controller)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    type(ab_time_scheme_t) :: ab_scheme_obj
    character(len=1000) :: msg
    real(kind=rp) :: ab_coeffs(4), dt_history(10)
    real(kind=rp) :: beta(0:3)
    integer :: nadv, i, k, row_g
    real(kind=rp) :: F_fluid(6)
    real(kind=dp) :: start_time_s, end_time_s, step_time_s
    logical :: iter_verbose

    nadv = this%ext_bdf%nadv

    do i = 0, nadv
       beta(i) = this%ext_bdf%diffusion_coeffs%x(i+1)
       if (i .ge. 1) beta(i) = -beta(i)
    end do

    call rzero(ab_coeffs, 4)
    dt_history(1) = time%dt
    dt_history(2) = time%dtlag(1)
    dt_history(3) = time%dtlag(2)
    call ab_scheme_obj%compute_coeffs(ab_coeffs, dt_history, nadv)

    ! Calculate the new displacements at the current time-step
    ! using AB-k.
    ! This is the "relative displacement" w.r.t the body.
    ! I have used AB instead of extrapolation for force to
    ! have consistency with how the pivot and mesh is updated.
    ! Otherwise (I think) we will have drift between displacement
    ! and pivot point.
    ! disp_rel <- disp_rel + dt*ab_coeffs(1)*body_vel
    !      + sum_{j=2}^{nadv} dt*ab_coeffs(j)*body_vel_lag(:,j)
    do i = 1, this%nbodies_fsi
       call this%ale%scheme%integrate_6dof(this%fsi_bodies(i)%disp_rel, &
            this%fsi_bodies(i)%body_vel, time, nadv, &
            v6_lag = this%fsi_bodies(i)%body_vel_lag, ab_coeffs = ab_coeffs)
    end do

    ! Add prescribed frame motion terms.
    ! I should double check if it is the right place to put this. I think it is.
    ! But, should become sure.
    call this%query_frame_prescribed_motion(time, beta, nadv)

    ! Standard fluid step
    start_time_s = MPI_WTIME()

    ! Fluid standard step, using final velocity from previous step,
    ! which includes the FSI correction and prescribed motions.
    ! We advance the mesh only here.
    ! We should not skip the mesh velocity here, since we need to update the
    ! rotation matrix.
    ! I think I can move compute_rotationm matrics inside
    ! advance_mesh_explicit. I need to remember if there was any reason that
    ! I put it there at first place.
    ! Calling update_mesh_velocity here should be totally harmless.
    call this%fluid_pnpn_t%step(time, dt_controller)

    ! Add FSI structural terms at current time step.
    ! Rotation matrix etc should be updated at this point.
    ! I moved this after the above step, since the rotation matrix should for
    ! the current time step.
    ! Need to verify more.
    call this%calc_fsi_terms(time, beta, nadv)

    end_time_s = MPI_WTIME()
    step_time_s = end_time_s - start_time_s
    this%total_elapsed_s = this%total_elapsed_s + step_time_s

    write(msg, '(A, E15.7, A, I0, A, E15.7)') "Standard step time (s):  ", &
         step_time_s, "  Step: ", time%tstep, "  time: ", time%t
    call neko_log%message(trim(msg))
    call neko_log%message(' ')

    ! Store the standard solution.
    call field_copy(this%u_s, this%u)
    call field_copy(this%v_s, this%v)
    call field_copy(this%w_s, this%w)
    call field_copy(this%p_s, this%p)

    ! Compute fluid forces/torques
    do i = 1, this%nbodies_fsi
       call this%fsi_bodies(i)%force_monitor%compute_(time)
       F_fluid(1:3) = this%fsi_bodies(i)%force_monitor%total_force
       F_fluid(4:6) = this%fsi_bodies(i)%force_monitor%total_torque

       ! Fill B_global with fluid forces (F_s)
       do k = 1, 6
          row_g = this%fsi_dof_map(i, k)
          if (row_g > 0) then
             this%B_global(row_g) = this%B_global(row_g) + F_fluid(k)
          end if
       end do
    end do

    if (.not. this%skip_greens_solve) then
       call fluid_fsi_compute_greens_functions(this, time, dt_controller)
    else
       call neko_log%message("Weak coupling enabled: " // &
            "Skipping Green's function fluid feedback.")
    end if

    call neko_log%message(' ')
    write(msg, '(A, E15.7, A, I0, A, E15.7)') &
         "Standard's step total elapsed time (s):  ", &
         this%total_elapsed_s, "  Step: ", time%tstep, "  time: ", time%t
    call neko_log%message(trim(msg))
    write(msg, '(A, E15.7, A, I0, A, E15.7)') &
         "Green's step total elapsed time (s):  ", &
         this%total_elapsed_g, "  Step: ", time%tstep, "  time: ", time%t
    call neko_log%message(trim(msg))
    call neko_log%message(' ')

    ! Calculate all FSI corrections: M_global * X_sol = B_global.
    ! With current algorithm, the correction is actually just
    ! the velocity change compared to the previous time step.
    ! X_sol enters holding the previous step's solution, which serves as
    ! the initial guess of the fixed-point loop below.
    if (this%total_active_dofs > 0) then
       iter_verbose = this%non_linear_correction_term .or. &
            this%has_user_structural_terms
       call fsi_solve_structure(this%nbodies_fsi, this%fsi_bodies, &
            this%fsi_dof_map, this%total_active_dofs, this%M_global, &
            this%B_global, this%X_sol, this%ale%body_rot_matrices, time, &
            this%gravity_vec, this%non_linear_correction_term, &
            this%user_fsi_structural_terms, iter_verbose)
       if (iter_verbose) call neko_log%message(' ')
    end if

    call fluid_fsi_apply_correction(this)
    call fluid_fsi_update_bodies(this, nadv)

    ! Calculate Final Velocity (FSI + Prescribed)
    ! This velocity will be used as the "guessed" velocity
    ! for the next time step, and also for the ALE mesh update.
    call this%ale%update_mesh_velocity(this%c_Xh, time, &
         override_ids = this%batch_ids, &
         override_trans = this%batch_trans, &
         override_ang = this%batch_ang, &
         out_prescribed_vels = this%temp_prescribed_vels, &
         mode = ALE_VEL_SUPERPOSE)
    do i = 1, this%nbodies_fsi
       this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) = &
            this%temp_prescribed_vels(:, i)
    end do

    call fsi_prep_checkpoint(this%nbodies_fsi, this%fsi_bodies, &
         this%global_disp_rel, &
         this%global_body_vel, &
         this%global_body_vel_lag, &
         this%global_moving_frame_presc_vel)

    call this%log_fsi_results(time)
    call this%ale%log_pivot(time)
    call this%ale%log_rot_angles(time)

  end subroutine fluid_fsi_step

  !> Solve the Green's function problem of every active DOF and add the
  !! resulting forces to the structural system.
  subroutine fluid_fsi_compute_greens_functions(this, time, dt_controller)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    integer :: i, j, k, k_row, row_g, col_g
    real(kind=dp) :: start_time_g, end_time_g, step_time_g
    character(len=1000) :: msg

    do j = 1, this%nbodies_fsi
       do k = 1, 6
          ! Solve Green's function for active DOFs only.
          col_g = this%fsi_dof_map(j, k)
          if (col_g == 0) cycle

          ! Use the last impulse response fields for initial guess.
          call field_copy(this%u, this%u_g(col_g))
          call field_copy(this%v, this%v_g(col_g))
          call field_copy(this%w, this%w_g(col_g))
          call field_copy(this%p, this%p_g(col_g))

          ! Setup Perturbation
          this%batch_ids(1) = this%fsi_bodies(j)%ale_id
          this%batch_trans(:,1) = 0.0_rp
          this%batch_ang(:,1) = 0.0_rp

          if (k <= 3) then
             ! translational DOF
             this%batch_trans(k, 1) = 1.0_rp
          else
             ! rotational DOF
             this%batch_ang(k-3, 1) = 1.0_rp
          end if

          ! Mode 1: Set rigid body vels to zero, then apply impulse.
          call this%ale%update_mesh_velocity(this%c_Xh, time, &
               override_ids = this%batch_ids(1:1), &
               override_trans = this%batch_trans(:,1:1), &
               override_ang = this%batch_ang(:,1:1), &
               mode = ALE_VEL_EXCLUSIVE)

          start_time_g = MPI_WTIME()
          call fluid_fsi_greens_solve(this, time, dt_controller, col_g)
          end_time_g = MPI_WTIME()
          step_time_g = end_time_g - start_time_g
          this%total_elapsed_g = this%total_elapsed_g + step_time_g

          write(msg, '(A, E15.7, A, I0, A, E15.7)') &
               "Green's step time (s):  ", step_time_g, "  Step: ", &
               time%tstep, "  time: ", time%t
          call neko_log%message(trim(msg))
          call neko_log%message(' ')

          ! Save Green's function response.
          call field_copy(this%u_g(col_g), this%u)
          call field_copy(this%v_g(col_g), this%v)
          call field_copy(this%w_g(col_g), this%w)
          call field_copy(this%p_g(col_g), this%p)

          ! Fill M matrix with Impulse forces/torques (F_g)
          ! Here, we add the cross-coupling forces on all bodies, on all DOFs.
          do i = 1, this%nbodies_fsi
             call this%fsi_bodies(i)%force_monitor%compute_(time)
             do k_row = 1, 6
                row_g = this%fsi_dof_map(i, k_row)
                if (row_g > 0) then
                   if (k_row <= 3) then
                      this%M_global(row_g, col_g) = &
                           this%M_global(row_g, col_g) - &
                           this%fsi_bodies(i)%&
                           force_monitor%total_force(k_row)
                   else
                      this%M_global(row_g, col_g) = &
                           this%M_global(row_g, col_g) - &
                           this%fsi_bodies(i)%&
                           force_monitor%total_torque(k_row-3)
                   end if
                end if
             end do
          end do
       end do
    end do
  end subroutine fluid_fsi_compute_greens_functions

  !> Restore the standard solution and, for strong coupling, add the
  !! Green's function correction: u = u_s + sum( X_sol(k) * u_g(k) ).
  subroutine fluid_fsi_apply_correction(this)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    integer :: n, idx_g

    n = this%dm_Xh%size()

    ! Restore standard fields (this is the final state for weak coupling)
    call field_copy(this%u, this%u_s)
    call field_copy(this%v, this%v_s)
    call field_copy(this%w, this%w_s)
    call field_copy(this%p, this%p_s)

    ! Add FSI correction to the fluid solution only if strong coupling
    ! u = u_s + sum( X_sol(k) * u_g(k) )
    if (.not. this%skip_greens_solve) then
       if (this%total_active_dofs > 0) then
          do idx_g = 1, this%total_active_dofs
             call field_add2s2(this%u, this%u_g(idx_g), this%X_sol(idx_g), n)
             call field_add2s2(this%v, this%v_g(idx_g), this%X_sol(idx_g), n)
             call field_add2s2(this%w, this%w_g(idx_g), this%X_sol(idx_g), n)
             call field_add2s2(this%p, this%p_g(idx_g), this%X_sol(idx_g), n)
          end do
       end if
    end if
  end subroutine fluid_fsi_apply_correction

  !> Update the body velocities with the structural solution and shift
  !! their history.
  subroutine fluid_fsi_update_bodies(this, nadv)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    integer, intent(in) :: nadv
    integer :: i, k, row_g

    do i = 1, this%nbodies_fsi
       do k = 1, 6
          row_g = this%fsi_dof_map(i, k)
          if (row_g > 0) then
             ! Corrected fsi_body velocity.
             ! Note that body_vel_lag(k, 1) = body_vel(k).
             this%fsi_bodies(i)%body_vel(k) = this%X_sol(row_g) + &
                  this%fsi_bodies(i)%body_vel_lag(k, 1)
          end if
       end do

       ! Update History
       do k = nadv, 2, -1
          this%fsi_bodies(i)%body_vel_lag(:, k) = &
               this%fsi_bodies(i)%body_vel_lag(:, k-1)
       end do
       this%fsi_bodies(i)%body_vel_lag(:, 1) = this%fsi_bodies(i)%body_vel

       ! Fill Batch Arrays
       this%batch_ids(i) = this%fsi_bodies(i)%ale_id
       this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
       this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
    end do
  end subroutine fluid_fsi_update_bodies

  !> Solve the Green's function problem of one active DOF.
  subroutine fluid_fsi_greens_solve(this, time, dt_controller, col_g)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    integer, intent(in) :: col_g
    integer :: prs_max_iter

    call profiler_start_region('Fluid', 1)

    ! No explicit forcing or history in the Green's function problem.
    call field_cfill(this%f_x, 0.0_rp)
    call field_cfill(this%f_y, 0.0_rp)
    call field_cfill(this%f_z, 0.0_rp)
    call field_cfill(this%u_e, 0.0_rp)
    call field_cfill(this%v_e, 0.0_rp)
    call field_cfill(this%w_e, 0.0_rp)

    call neko_log%message(" ")
    call neko_log%message("--------Green's Step----------")

    prs_max_iter = this%ksp_prs%max_iter
    this%ksp_prs%max_iter = this%greens_prs_max_iter
    this%greens_mode = .true.
    call this%solve(time, dt_controller, this%proj_prs_green(col_g), &
         this%proj_vel_green(col_g))
    this%greens_mode = .false.
    this%ksp_prs%max_iter = prs_max_iter

    call profiler_end_region('Fluid', 1)
  end subroutine fluid_fsi_greens_solve

  !> Apply the velocity boundary conditions; those of the Green's function
  !! problem while one is being solved.
  subroutine fluid_fsi_bc_apply_vel(this, time, strong)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    logical, intent(in) :: strong
    class(bc_t), pointer :: b
    integer :: i

    if (.not. this%greens_mode) then
       call this%fluid_pnpn_t%bc_apply_vel(time, strong)
       return
    end if

    b => null()

    call this%bcs_vel_green%apply_vector(this%u%x, this%v%x, this%w%x, &
         this%dm_Xh%size(), time, strong = .true.)

    call rotate_cyc(this%u, this%v, this%w, 1, this%c_Xh)
    call this%gs_Xh%op(this%u, GS_OP_MIN, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call this%gs_Xh%op(this%v, GS_OP_MIN, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call this%gs_Xh%op(this%w, GS_OP_MIN, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call rotate_cyc(this%u, this%v, this%w, 0, this%c_Xh)

    call this%bcs_vel_green%apply_vector(this%u%x, this%v%x, this%w%x, &
         this%dm_Xh%size(), time, strong = .true.)

    call rotate_cyc(this%u%x, this%v%x, this%w%x, 1, this%c_Xh)
    call this%gs_Xh%op(this%u, GS_OP_MAX, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call this%gs_Xh%op(this%v, GS_OP_MAX, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call this%gs_Xh%op(this%w, GS_OP_MAX, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call rotate_cyc(this%u%x, this%v%x, this%w%x, 0, this%c_Xh)

    do i = 1, this%bcs_vel_green%size()
       b => this%bcs_vel_green%get(i)
       b%updated = .false.
    end do
    nullify(b)
  end subroutine fluid_fsi_bc_apply_vel

  !> Apply the pressure boundary conditions; those of the Green's function
  !! problem while one is being solved.
  subroutine fluid_fsi_bc_apply_prs(this, time)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    class(bc_t), pointer :: b
    integer :: i

    if (.not. this%greens_mode) then
       call this%fluid_pnpn_t%bc_apply_prs(time)
       return
    end if

    b => null()

    call this%bcs_prs_green%apply(this%p, time)

    call this%gs_Xh%op(this%p, GS_OP_MIN, glb_cmd_event)
    call device_event_sync(glb_cmd_event)

    call this%bcs_prs_green%apply(this%p, time)

    call this%gs_Xh%op(this%p, GS_OP_MAX, glb_cmd_event)
    call device_event_sync(glb_cmd_event)

    do i = 1, this%bcs_prs_green%size()
       b => this%bcs_prs_green%get(i)
       b%updated = .false.
    end do
    nullify(b)
  end subroutine fluid_fsi_bc_apply_prs

  !> Build the boundary conditions of the Green's function problems from
  !! the fluid ones.
  subroutine fluid_fsi_setup_green_bcs(this)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    class(bc_t), pointer :: bc_i
    integer :: i

    call this%bcs_vel_green%init()
    do i = 1, this%bcs_vel%size()
       bc_i => this%bcs_vel%get(i)
       call fluid_fsi_green_vel_bc(this, bc_i)
    end do

    call this%bcs_prs_green%init()
    do i = 1, this%bcs_prs%size()
       bc_i => this%bcs_prs%get(i)
       call fluid_fsi_green_prs_bc(this, bc_i)
    end do
  end subroutine fluid_fsi_setup_green_bcs

  !> Subroutine to setup Green's function Velocity BCs
  subroutine fluid_fsi_green_vel_bc(this, bc_i)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    class(bc_t), intent(inout) :: bc_i
    class(bc_t), pointer :: bc_green

    bc_green => null()

    select type (orig => bc_i)

       ! Moving Wall (FSI Body)
       ! MUST remain a 'no_slip_t' so it tracks mesh velocity.
    type is (no_slip_t)
       if (orig%is_moving) then
          allocate(no_slip_t :: bc_green)
          ! We must verify if no_slip needs init via JSON or manual components
          select type (n => bc_green)
          type is (no_slip_t)
             call n%zero_dirichlet_t%init_from_components(this%c_Xh)
             n%is_moving = .true.

             if (associated(orig%wx)) n%wx => orig%wx
             if (associated(orig%wy)) n%wy => orig%wy
             if (associated(orig%wz)) n%wz => orig%wz
          end select
       else
          ! Stationary Wall -> Zero Dirichlet
          allocate(zero_dirichlet_t :: bc_green)
       end if

       ! Inlets / User Dirichlet -> Zero Dirichlet (0.0)
    type is (inflow_t)
       allocate(zero_dirichlet_t :: bc_green)
    type is (field_dirichlet_vector_t)
       allocate(zero_dirichlet_t :: bc_green)
    type is (wall_model_bc_t)
       allocate(zero_dirichlet_t :: bc_green)

       ! Constraints (Symmetry, etc.) -> Keep Same Pointer
       !
       ! Only the axis-aligned constraint bcs can be cloned here. The mixed
       ! bcs (symmetry_t, non_normal_t, shear_stress_t and its descendants)
       ! are resolved globally by the coupled_vector_bc_projector_t and carry
       ! no usable constraint of their own: an unregistered clone has an empty
       ! resolved_msk and would silently apply nothing.
    type is (symmetry_aligned_t)
       allocate(symmetry_aligned_t :: bc_green)
       select type (s => bc_green)
       type is (symmetry_aligned_t)
          call s%init_from_components(this%c_Xh)
       end select

    type is (symmetry_t)
       call neko_error("FSI Green's functions do not support the coupled " // &
            "symmetry boundary condition. Disable the full stress " // &
            "formulation to get the axis-aligned variant.")

    type is (non_normal_t)
       call neko_error("FSI Green's functions do not support the coupled " // &
            "normal_outflow boundary condition. Disable the full stress " // &
            "formulation to get the axis-aligned variant.")

    type is (shear_stress_t)
       call neko_error("FSI Green's functions do not support the " // &
            "shear_stress boundary condition, which is now a mixed bc " // &
            "resolved by the coupled velocity projector.")

       ! Outflow -> Natural (Null)
    class default
       ! bc_green => null()
    end select

    ! If we created a valid Green's BC, append it
    if (associated(bc_green)) then

       ! Initialize if it's a new object
       select type (z => bc_green)
       type is (zero_dirichlet_t)
          call z%init_from_components(this%c_Xh)
       end select

       call bc_green%mark_facets(bc_i%marked_facet)
       call bc_green%finalize()

       call this%bcs_vel_green%append(bc_green)
    end if
  end subroutine fluid_fsi_green_vel_bc

  !> Helper subroutine to setup Green's function Pressure BCs
  subroutine fluid_fsi_green_prs_bc(this, bc_i)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    class(bc_t), intent(inout) :: bc_i
    class(bc_t), pointer :: bc_p_green

    bc_p_green => null()

    if (bc_i%bc_type .eq. BC_DIRICHLET) then
       allocate(zero_dirichlet_t :: bc_p_green)
       select type (z => bc_p_green)
       type is (zero_dirichlet_t)
          call z%init_from_components(this%c_Xh)
          call z%mark_facets(bc_i%marked_facet)
          call z%finalize()
       end select
       call this%bcs_prs_green%append(bc_p_green)
    end if

  end subroutine fluid_fsi_green_prs_bc

  !> Destructor.
  subroutine fluid_fsi_free(this)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    class(bc_t), pointer :: bc
    integer :: i, k

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

    if (allocated(this%u_g)) then
       do k = 1, this%total_active_dofs
          call this%u_g(k)%free()
          call this%v_g(k)%free()
          call this%w_g(k)%free()
          call this%p_g(k)%free()
       end do
       deallocate(this%u_g, this%v_g, this%w_g, this%p_g)
    end if

    if (allocated(this%proj_prs_green)) then
       do k = 1, this%total_active_dofs
          call this%proj_prs_green(k)%free()
          call this%proj_vel_green(k)%free()
       end do
       deallocate(this%proj_prs_green)
       deallocate(this%proj_vel_green)
    end if

    if (allocated(this%batch_ids)) deallocate(this%batch_ids)
    if (allocated(this%batch_trans)) deallocate(this%batch_trans)
    if (allocated(this%batch_ang)) deallocate(this%batch_ang)
    if (allocated(this%temp_prescribed_vels)) &
         deallocate(this%temp_prescribed_vels)

    call this%u_s%free()
    call this%v_s%free()
    call this%w_s%free()
    call this%p_s%free()

    do i = 1, this%bcs_vel_green%size()
       bc => this%bcs_vel_green%get(i)
       call bc%free()
    end do
    call this%bcs_vel_green%free()

    do i = 1, this%bcs_prs_green%size()
       bc => this%bcs_prs_green%get(i)
       call bc%free()
    end do
    call this%bcs_prs_green%free()

    ! Free the base fluid_pnpn_t
    call this%fluid_pnpn_t%free()

  end subroutine fluid_fsi_free

  !> Get the prescribed velocity of the moving frame of each body and
  !! compute its acceleration.
  subroutine fluid_fsi_query_frame_prescribed_motion(this, time, beta, nadv)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    integer, intent(in) :: nadv
    integer :: i, k

    ! Shift history back
    ! index 0 is the current time prescribed velocity, 1-3 are the history.
    do i = 1, this%nbodies_fsi
       do k = nadv, 1, -1
          this%fsi_bodies(i)%moving_frame_presc_vel(:, k) = &
               this%fsi_bodies(i)%moving_frame_presc_vel(:, k-1)
       end do
    end do

    do i = 1, this%nbodies_fsi
       this%batch_ids(i) = this%fsi_bodies(i)%ale_id
    end do

    ! Here we only get the prescribed motion for frame of movment.
    call this%ale%update_mesh_velocity(this%c_Xh, time, &
         override_ids = this%batch_ids, &
         out_prescribed_vels = this%temp_prescribed_vels, &
         mode = ALE_VEL_QUERY)

    ! current velocity of the moving frame
    do i = 1, this%nbodies_fsi
       this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) = &
            this%temp_prescribed_vels(:, i)
    end do

    do i = 1, this%nbodies_fsi
       this%fsi_bodies(i)%moving_frame_presc_acc = 0.0_rp
       do k = 0, nadv
          this%fsi_bodies(i)%moving_frame_presc_acc = &
               this%fsi_bodies(i)%moving_frame_presc_acc + &
               (beta(k) * &
               this%fsi_bodies(i)%moving_frame_presc_vel(:, k)) / time%dt
       end do
    end do
  end subroutine fluid_fsi_query_frame_prescribed_motion

  !> Assemble the structural and inertial terms of the FSI system.
  subroutine assemble_fsi_structural_inertial_terms(this, time, beta, nadv)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    integer, intent(in) :: nadv
    real(kind=rp) :: gamma

    gamma = beta(0) / time%dt

    ! Here we fill the M_global and B_global
    ! using the contribution from structure and also
    ! bodies' inertial motion from previous time steps.
    call assemble_structural_inertial_terms(this%nbodies_fsi, &
         this%fsi_bodies, &
         this%fsi_dof_map, this%M_global, this%B_global, &
         this%ale%body_rot_matrices, &
         time, gamma, beta, nadv, this%gravity_vec, &
         this%user_fsi_body_params)

  end subroutine assemble_fsi_structural_inertial_terms

  !> Log the state of the FSI bodies.
  subroutine fluid_fsi_log_results(this, time)
    class(fluid_pnpn_fsi_greens_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    character(len=1024) :: msg
    character(len=128) :: fmt_res
    integer :: i, k, row_g
    real(kind=rp) :: corr_coef(6)

    if (this%nbodies_fsi == 0) return

    call neko_log%message("---------FSI Results----------")

    if (this%res_long_print) then
       fmt_res = '(A, I0, A, ES23.15, A, A, A, 3(ES22.15, :, 2X))'
    else
       fmt_res = '(A, I0, A, ES17.10, A, A, A, 3(ES17.10, :, 2X))'
    end if

    call neko_log%message("variable, time step, time, body, x_val, " // &
         "y_val, z_val")

    do i = 1, this%nbodies_fsi
       ! Correction Coefficients (X_sol) for this body
       do k = 1, 6
          row_g = this%fsi_dof_map(i, k)
          if (row_g > 0) then
             corr_coef(k) = this%X_sol(row_g)
          else
             corr_coef(k) = 0.0_rp
          end if
       end do

       ! Linear Displacement (x, y, z)
       write(msg, fmt_res) &
            "FSI_DISP_L  ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", &
            this%fsi_bodies(i)%disp_rel(1:3)
       call neko_log%message(trim(msg))

       ! Angular Displacement (rx, ry, rz)
       write(msg, fmt_res) &
            "FSI_DISP_A  ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", &
            this%fsi_bodies(i)%disp_rel(4:6)
       call neko_log%message(trim(msg))

       ! Linear Velocity (x, y, z)
       write(msg, fmt_res) &
            "FSI_VEL_L   ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", &
            this%fsi_bodies(i)%body_vel(1:3)
       call neko_log%message(trim(msg))

       ! Angular Velocity (rx, ry, rz)
       write(msg, fmt_res) &
            "FSI_VEL_A   ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", &
            this%fsi_bodies(i)%body_vel(4:6)
       call neko_log%message(trim(msg))

       ! Linear Correction Coef
       write(msg, fmt_res) &
            "FSI_CORR_L  ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", &
            corr_coef(1:3)
       call neko_log%message(trim(msg))

       ! Angular Correction Coef
       write(msg, fmt_res) &
            "FSI_CORR_A  ", time%tstep, "  ", time%t, "  ", &
            trim(this%fsi_bodies(i)%name), "  ", &
            corr_coef(4:6)
       call neko_log%message(trim(msg))

    end do
    call neko_log%message(" ")

  end subroutine fluid_fsi_log_results

  !> Restart from a previous solution.
  subroutine fluid_fsi_restart(this, chkp)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    type(time_state_t) :: t_restart
    real(kind=dp), pointer :: tlag(:), dtlag(:)
    integer :: i, n

    call chkp%get_time_history(tlag, dtlag)

    n = this%u%dof%size()

    ! Restart the base fluid_pnpn_t
    call this%fluid_pnpn_t%restart(chkp)

    ! Restore FSI specific arrays
    if (this%if_fsi .and. this%nbodies_fsi > 0) then
       call fsi_restart_restore(this%nbodies_fsi, this%fsi_bodies, &
            this%global_disp_rel, this%global_body_vel, &
            this%global_body_vel_lag, &
            this%global_moving_frame_presc_vel)

       t_restart%t = chkp%t
       t_restart%tstep = 0
       t_restart%dt = dtlag(1)

       do i = 1, this%nbodies_fsi
          this%batch_ids(i) = this%fsi_bodies(i)%ale_id
          this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
          this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
       end do

       call this%ale%update_mesh_velocity(this%c_Xh, t_restart, &
            override_ids = this%batch_ids, &
            override_trans = this%batch_trans, &
            override_ang = this%batch_ang, &
            mode = ALE_VEL_SUPERPOSE)
    end if
  end subroutine fluid_fsi_restart

end module fluid_pnpn_fsi_greens
