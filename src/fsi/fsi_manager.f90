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
!> Setup of the FSI rigid bodies, the dense structural solve and the
!! checkpointing of the FSI state.
module fsi_manager
  use fsi_dynamics, only : fsi_body_t, add_fsi_non_linear_matrices, &
       add_fsi_user_structural_terms
  use fsi_body_params, only : FSI_INERTIA_ABOUT_PIVOT, &
       FSI_INERTIA_ABOUT_COM, params_inertia_about_pivot, validate_body_params
  use field, only : field_t
  use field_math, only : field_rzero
  use num_types, only : rp
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use utils, only : neko_error
  use logger, only : neko_log
  use time_state, only : time_state_t
  use user_intf, only : user_fsi_structural_terms_intf
  use checkpoint, only : chkp_t
  use checkpoint_payload, only : checkpoint_payload_t
  use ale_manager, only : ale_manager_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use projection, only : projection_t
  use projection_vel, only : projection_vel_t
  implicit none
  private

  public :: fsi_manager_init
  public :: linsolve_dense
  public :: fsi_solve_structure
  public :: fsi_register_checkpoint
  public :: fsi_prep_checkpoint
  public :: fsi_restart_restore

contains

  !> Read the FSI settings and bodies from the case file and set up the FSI
  !! data structures.
  subroutine fsi_manager_init(params, ale, c_Xh, dm_Xh, &
       if_fsi, nbodies_fsi, bodies, fsi_dof_map, &
       total_active_dofs, M_global, B_global, X_sol, &
       u_g, v_g, w_g, p_g, res_long_print, gravity_vec, &
       proj_prs_green, proj_vel_green, global_disp_rel, &
       global_body_vel, global_body_vel_lag, &
       global_moving_frame_presc_vel, skip_greens_solve, &
       non_linear_correction_term, global_body_acc, global_frame_acc)
    type(json_file), target, intent(inout) :: params
    type(ale_manager_t), target, intent(inout) :: ale
    type(coef_t), intent(inout) :: c_Xh
    type(dofmap_t), intent(in) :: dm_Xh
    logical, intent(out) :: if_fsi
    integer, intent(out) :: nbodies_fsi
    type(fsi_body_t), allocatable, intent(out) :: bodies(:)
    integer, allocatable, intent(out) :: fsi_dof_map(:,:)
    integer, intent(out) :: total_active_dofs
    real(kind=rp), allocatable, intent(out) :: M_global(:,:), B_global(:)
    real(kind=rp), allocatable, intent(out) :: X_sol(:)
    type(field_t), allocatable, intent(out) :: u_g(:), v_g(:), w_g(:)
    type(field_t), allocatable, intent(out) :: p_g(:)
    logical, intent(out) :: res_long_print
    real(kind=rp), intent(out) :: gravity_vec(3)
    type(projection_t), allocatable, intent(out) :: proj_prs_green(:)
    type(projection_vel_t), allocatable, intent(out) :: proj_vel_green(:)
    real(kind=rp), allocatable, intent(out) :: global_disp_rel(:)
    real(kind=rp), allocatable, intent(out) :: global_body_vel(:)
    real(kind=rp), allocatable, intent(out) :: global_body_vel_lag(:,:)
    real(kind=rp), allocatable, intent(out) :: &
         global_moving_frame_presc_vel(:,:)
    logical, intent(out) :: skip_greens_solve
    logical, intent(out) :: non_linear_correction_term
    !> Newmark previous-acceleration.
    real(kind=rp), allocatable, intent(out), optional :: global_body_acc(:)
    !> Newmark prescribed-frame previous-acceleration.
    real(kind=rp), allocatable, intent(out), optional :: global_frame_acc(:)
    ! Projection settings of the Green's function solves
    integer :: fsi_pr_projection_dim, fsi_vel_projection_dim
    integer :: fsi_pr_projection_activ_step, fsi_vel_projection_activ_step
    logical :: fsi_pr_projection_reortho
    integer :: i, k, n_bodies
    character(len=200) :: name_buf, name_buf2, log_buf, field_name
    real(kind=rp) :: center_dummy(3), force_scale
    real(kind=rp) :: com_position(3), cob_position(3)
    real(kind=rp), allocatable :: temp_vec(:)
    character(len=:), allocatable :: temp_str
    character(len=20) :: coupling_mode
    logical :: log_forces, long_print

    center_dummy = 0.0_rp

    call neko_log%section("Fluid-Structure Interaction")

    if (.not. params%valid_path('case.fluid.fsi')) then
       call neko_error("Missing 'case.fluid.fsi' block in case.json")
    end if

    call json_get(params, 'case.fluid.fsi.enabled', if_fsi)
    if (.not. if_fsi) then
       call neko_error("Scheme pnpn_fsi: FSI block present, but " // &
            "'fsi.enabled': false in case file")
    end if

    call json_get_or_default(params, 'case.fluid.fsi.log_forces', log_forces, &
         .true.)
    call json_get_or_default(params, 'case.fluid.fsi.long_print', long_print, &
         .false.)
    call json_get_or_default(params, 'case.fluid.fsi.results_long_print', &
         res_long_print, .false.)
    call json_get_or_default(params, 'case.fluid.fsi.force_scale', &
         force_scale, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.fsi.skip_greens_solve', &
         skip_greens_solve, .false.)
    call json_get_or_default(params, &
         'case.fluid.fsi.non_linear_correction_term', &
         non_linear_correction_term, .false.)

    call json_get(params, 'case.fluid.fsi.coupling', temp_str)
    if (trim(temp_str) .eq. 'subiteration') then
       skip_greens_solve = .true.
       coupling_mode = "subiteration"
    else if (trim(temp_str) .eq. 'greens') then
       coupling_mode = "greens"
    else
       call neko_error("FSI: Unknown coupling mode: " // trim(temp_str) // &
            ". Must be 'subiteration' or 'greens'.")
    end if

    if ((.not. skip_greens_solve) .and. coupling_mode == "greens") then
       call neko_log%message(" ")
       call neko_log%message( &
            "-----------------------------------------------------")
       call neko_log%message("FSI Coupling: Greens + Strong Coupling " // &
            "(Implicit FSI)")
       call neko_log%message( &
            "-----------------------------------------------------")
       call neko_log%message(" ")
    else if ((skip_greens_solve) .and. coupling_mode == "greens") then
       call neko_log%message(" ")
       call neko_log%message( &
            "---------------------------------------------------")
       call neko_log%message("FSI Coupling: Greens + Weak Coupling " // &
            "(Explicit FSI)")
       call neko_log%message( &
            "---------------------------------------------------")
       call neko_log%message(" ")
    else if (coupling_mode == "subiteration") then
       call neko_log%message(" ")
       call neko_log%message( &
            "-----------------------------------------------------")
       call neko_log%message("FSI Coupling: Subiteration")
       call neko_log%message( &
            "-----------------------------------------------------")
       call neko_log%message(" ")
    end if

    call json_get_or_default(params, &
         'case.fluid.fsi.pressure_solver.projection_space_size', &
         fsi_pr_projection_dim, 1)

    call json_get_or_default(params, &
         'case.fluid.fsi.pressure_solver.projection_hold_steps', &
         fsi_pr_projection_activ_step, 5)

    call json_get_or_default(params, &
         'case.fluid.fsi.pressure_solver.projection_reorthogonalize_basis', &
         fsi_pr_projection_reortho, .false.)

    call json_get_or_default(params, &
         'case.fluid.fsi.velocity_solver.projection_space_size', &
         fsi_vel_projection_dim, 0)

    call json_get_or_default(params, &
         'case.fluid.fsi.velocity_solver.projection_hold_steps', &
         fsi_vel_projection_activ_step, 5)

    gravity_vec = 0.0_rp
    if (params%valid_path('case.fluid.fsi.gravity_vec')) then
       call json_get(params, 'case.fluid.fsi.gravity_vec', temp_vec, &
            expected_size = 3)
       gravity_vec = temp_vec
    end if

    call params%info('case.fluid.fsi.bodies', n_children = n_bodies)
    nbodies_fsi = n_bodies
    if (nbodies_fsi == 0) then
       call neko_error("FSI: 'case.fluid.fsi.bodies' is empty (no FSI " // &
            "bodies defined)")
    end if

    ! Global FSI Logging
    write(log_buf, '(A,3(ES13.6,1X))') ' Gravity Vector  : ', gravity_vec
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') ' Force Scale     : ', force_scale
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') ' Number of Bodies: ', nbodies_fsi
    call neko_log%message(log_buf)
    write(log_buf, '(A,L1)') ' Non-Linear Corr.: ', non_linear_correction_term
    call neko_log%message(log_buf)
    call neko_log%message(' ')

    allocate(bodies(nbodies_fsi))
    allocate(fsi_dof_map(nbodies_fsi, 6))
    allocate(global_disp_rel(6 * nbodies_fsi))
    allocate(global_body_vel(6 * nbodies_fsi))

    if (present(global_body_acc)) then
       allocate(global_body_acc(6 * nbodies_fsi))
       global_body_acc = 0.0_rp
    end if

    if (present(global_frame_acc)) then
       allocate(global_frame_acc(6 * nbodies_fsi))
       global_frame_acc = 0.0_rp
    end if

    allocate(global_body_vel_lag(6 * nbodies_fsi, &
         lbound(bodies(1)%body_vel_lag, 2) : ubound(bodies(1)%body_vel_lag, 2)))

    allocate(global_moving_frame_presc_vel(6 * nbodies_fsi, &
         lbound(bodies(1)%moving_frame_presc_vel, 2) : &
         ubound(bodies(1)%moving_frame_presc_vel, 2)))

    global_disp_rel = 0.0_rp
    global_body_vel = 0.0_rp
    global_body_vel_lag = 0.0_rp
    global_moving_frame_presc_vel = 0.0_rp

    total_active_dofs = 0

    do i = 1, nbodies_fsi
       call read_body(params, ale, i, bodies(i), com_position, cob_position)

       do k = 1, 6
          if (bodies(i)%active_dofs(k) == 1) then
             total_active_dofs = total_active_dofs + 1
             fsi_dof_map(i, k) = total_active_dofs
          else
             fsi_dof_map(i, k) = 0
          end if
       end do

       call log_body(bodies(i), com_position, cob_position, gravity_vec)

       write(name_buf, '(A,I0)') "fsi_force_monitor_", i
       write(name_buf2, '(A,A)') "fsi_body: ", trim(bodies(i)%name)

       call bodies(i)%force_monitor%init_common(name = trim(name_buf), &
            fluid_name = "fluid", zone_id = bodies(i)%zone_id, &
            zone_name = name_buf2, center = center_dummy, &
            scale = force_scale, coef = c_Xh, long_print = long_print, &
            center_type = 'pivot', full_log = log_forces)
    end do

    if (total_active_dofs .le. 0) then
       call neko_error("FSI: No active DOFs found!")
    end if

    allocate(M_global(total_active_dofs, total_active_dofs))
    allocate(B_global(total_active_dofs))
    allocate(X_sol(total_active_dofs))
    X_sol = 0.0_rp

    ! The Green's function storage is only needed when they are solved for
    if (.not. skip_greens_solve) then
       allocate(u_g(total_active_dofs))
       allocate(v_g(total_active_dofs))
       allocate(w_g(total_active_dofs))
       allocate(p_g(total_active_dofs))
       allocate(proj_prs_green(total_active_dofs))
       allocate(proj_vel_green(total_active_dofs))

       do k = 1, total_active_dofs
          write(field_name, '(A,I0)') 'u_g_', k
          call u_g(k)%init(dm_Xh, trim(field_name))
          write(field_name, '(A,I0)') 'v_g_', k
          call v_g(k)%init(dm_Xh, trim(field_name))
          write(field_name, '(A,I0)') 'w_g_', k
          call w_g(k)%init(dm_Xh, trim(field_name))
          write(field_name, '(A,I0)') 'p_g_', k
          call p_g(k)%init(dm_Xh, trim(field_name))

          call proj_prs_green(k)%init(dm_Xh%size(), fsi_pr_projection_dim, &
               fsi_pr_projection_activ_step, fsi_pr_projection_reortho)

          call proj_vel_green(k)%init(dm_Xh%size(), fsi_vel_projection_dim, &
               fsi_vel_projection_activ_step)
          call field_rzero(u_g(k))
          call field_rzero(v_g(k))
          call field_rzero(w_g(k))
          call field_rzero(p_g(k))
       end do
    end if

    call neko_log%end_section()

  end subroutine fsi_manager_init

  !> Read the parameters and initial state of body `i` from the case file.
  subroutine read_body(params, ale, i, body, com_position, cob_position)
    type(json_file), intent(inout) :: params
    type(ale_manager_t), intent(in) :: ale
    integer, intent(in) :: i
    type(fsi_body_t), intent(inout) :: body
    real(kind=rp), intent(out) :: com_position(3), cob_position(3)
    type(json_file) :: body_sub
    integer :: j, k, m
    real(kind=rp) :: pivot_init(3)
    real(kind=rp), allocatable :: temp_vec(:)
    integer, allocatable :: temp_vec_int(:)
    character(len=:), allocatable :: temp_str

    call json_extract_item(params, 'case.fluid.fsi.bodies', i, body_sub)

    if (body_sub%valid_path('name')) then
       call json_get(body_sub, 'name', temp_str)
       body%name = temp_str
    else
       write(body%name, '(A,I0)') 'fsi_body_', i
    end if

    call json_get(body_sub, 'zone_id', body%zone_id)

    body%ale_id = -1
    do j = 1, ale%config%nbodies
       if (any(ale%config%bodies(j)%zone_indices == body%zone_id)) then
          body%ale_id = j
          exit
       end if
    end do

    if (body%ale_id == -1) then
       call neko_error("FSI: Body " // trim(body%name) // &
            " zone_id not found in ALE bodies!")
    end if

    call json_get(body_sub, 'mass', body%prm%mass)

    body%prm%inertia = 0.0_rp
    if (body_sub%valid_path('I_xx_xy_xz_p') .or. &
         body_sub%valid_path('I_yx_yy_yz_p') .or. &
         body_sub%valid_path('I_zx_zy_zz_p')) then

       call json_get(body_sub, 'I_xx_xy_xz_p', temp_vec, expected_size = 3)
       body%prm%inertia(1,:) = temp_vec

       call json_get(body_sub, 'I_yx_yy_yz_p', temp_vec, expected_size = 3)
       body%prm%inertia(2,:) = temp_vec

       call json_get(body_sub, 'I_zx_zy_zz_p', temp_vec, expected_size = 3)
       body%prm%inertia(3,:) = temp_vec

       body%prm%inertia_ref = FSI_INERTIA_ABOUT_PIVOT

    else if (body_sub%valid_path('I_xx_xy_xz_com') .or. &
         body_sub%valid_path('I_yx_yy_yz_com') .or. &
         body_sub%valid_path('I_zx_zy_zz_com')) then

       call json_get(body_sub, 'I_xx_xy_xz_com', temp_vec, expected_size = 3)
       body%prm%inertia(1,:) = temp_vec

       call json_get(body_sub, 'I_yx_yy_yz_com', temp_vec, expected_size = 3)
       body%prm%inertia(2,:) = temp_vec

       call json_get(body_sub, 'I_zx_zy_zz_com', temp_vec, expected_size = 3)
       body%prm%inertia(3,:) = temp_vec

       body%prm%inertia_ref = FSI_INERTIA_ABOUT_COM

    else
       call neko_error("FSI: Body " // trim(body%name) // &
            " missing full inertia tensor inputs (I_xx_xy_xz_p OR " // &
            "I_xx_xy_xz_com)")
    end if

    call json_get(body_sub, 'active_dofs', temp_vec_int, expected_size = 6)
    body%active_dofs = temp_vec_int

    ! Absolute COM position is used only to compute the pivot offset.
    call json_get(body_sub, 'center_of_mass', temp_vec, expected_size = 3)
    com_position = temp_vec

    pivot_init = ale%ale_pivot(body%ale_id)%pos

    body%prm%offset_com = com_position - pivot_init

    cob_position = 0.0_rp
    body%prm%offset_cob = 0.0_rp
    body%prm%mass_disp = 0.0_rp
    if (params%valid_path('case.fluid.fsi.gravity_vec')) then
       ! IF Gravity exists, buoyancy is mandatory;
       ! otherwise, the physics is wrong!
       call json_get(body_sub, 'center_of_buoyancy', temp_vec, &
            expected_size = 3)
       cob_position = temp_vec

       body%prm%offset_cob = cob_position - pivot_init
       call json_get(body_sub, 'mass_disp', body%prm%mass_disp)
    end if

    ! Prescribed constant forcing at pivot location
    body%prm%F_prescribed_pivot = 0.0_rp
    if (body_sub%valid_path('F0_p')) then
       call json_get(body_sub, 'F0_p', temp_vec, expected_size = 6)
       body%prm%F_prescribed_pivot = temp_vec
    end if

    call json_get(body_sub, 'stiffness', temp_vec, expected_size = 3)
    body%prm%K_lin = temp_vec

    call json_get(body_sub, 'damping', temp_vec, expected_size = 3)
    body%prm%C_lin = temp_vec

    call json_get(body_sub, 'rot_stiffness', temp_vec, expected_size = 3)
    body%prm%K_ang = temp_vec

    call json_get(body_sub, 'rot_damping', temp_vec, expected_size = 3)
    body%prm%C_ang = temp_vec

    call json_get(body_sub, 'pos_equilibrium', temp_vec, expected_size = 6)
    body%prm%pos_eq = temp_vec

    ! Fail on physically impossible inputs.
    call validate_body_params(body%prm, body%name)

    call json_get(body_sub, 'initial_velocity', temp_vec, expected_size = 6)
    body%initial_vel = temp_vec

    ! To be sure initial vel is 0 for inactive DOFs.
    do k = 1, 6
       if (body%active_dofs(k) .eq. 0) then
          body%initial_vel(k) = 0.0_rp
       end if
    end do

    body%disp_rel = 0.0_rp
    body%body_vel_lag = 0.0_rp
    body%body_vel = body%initial_vel

    ! I think it should be just 3. For the first two time
    ! lag(2) and lag(3) does not matter since we do not use them.
    ! but it's important to keep them in the lag array since
    ! we do not touch lag(2) and lag(3) in updating lag arrays
    ! in first two time steps, so we will lose this if we don't
    ! fill them now.
    do m = 1, 3
       body%body_vel_lag(:,m) = body%initial_vel
    end do
  end subroutine read_body

  !> Log the parameters of one body.
  subroutine log_body(body, com_position, cob_position, gravity_vec)
    type(fsi_body_t), intent(in) :: body
    real(kind=rp), intent(in) :: com_position(3), cob_position(3)
    real(kind=rp), intent(in) :: gravity_vec(3)
    real(kind=rp) :: I_pivot_log(3,3)
    character(len=200) :: log_buf

    call neko_log%message('Registered Body : ' // trim(body%name))

    write(log_buf, '(A,I0,A,I0)') &
         '   Zone ID       : ', body%zone_id, ' | ALE ID : ', body%ale_id
    call neko_log%message(log_buf)

    write(log_buf, '(A,6(I2,1X))') '   Active DOFs   : ', body%active_dofs
    call neko_log%message(log_buf)

    write(log_buf, '(A,ES18.11)') '   Mass          : ', body%prm%mass
    call neko_log%message(log_buf)

    if (body%prm%inertia_ref == FSI_INERTIA_ABOUT_COM) then
       call neko_log%message('   Inertia Tensor (Input at COM) :')
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_xx_xy_xz  : ', body%prm%inertia(1,:)
       call neko_log%message(trim(log_buf))
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_yx_yy_yz  : ', body%prm%inertia(2,:)
       call neko_log%message(trim(log_buf))
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_zx_zy_zz  : ', body%prm%inertia(3,:)
       call neko_log%message(trim(log_buf))

       I_pivot_log = params_inertia_about_pivot(body%prm)
       call neko_log%message('   Inertia Tensor (Derived at Pivot) :')
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_xx_xy_xz  : ', I_pivot_log(1,:)
       call neko_log%message(trim(log_buf))
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_yx_yy_yz  : ', I_pivot_log(2,:)
       call neko_log%message(trim(log_buf))
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_zx_zy_zz  : ', I_pivot_log(3,:)
       call neko_log%message(trim(log_buf))
    else
       call neko_log%message('   Inertia Tensor (Input at Pivot) :')
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_xx_xy_xz  : ', body%prm%inertia(1,:)
       call neko_log%message(trim(log_buf))
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_yx_yy_yz  : ', body%prm%inertia(2,:)
       call neko_log%message(trim(log_buf))
       write(log_buf, '(A,3(ES13.6,1X))') &
            '     I_zx_zy_zz  : ', body%prm%inertia(3,:)
       call neko_log%message(trim(log_buf))
    end if

    write(log_buf, '(A,3(ES13.6,1X))') '   COM Position  : ', com_position
    call neko_log%message(log_buf)
    write(log_buf, '(A,3(ES13.6,1X))') &
         '   COM Offset    : ', body%prm%offset_com
    call neko_log%message(log_buf)

    if (any(abs(gravity_vec) > 0.0_rp)) then
       write(log_buf, '(A,3(ES13.6,1X))') '   COB Position  : ', cob_position
       call neko_log%message(log_buf)
       write(log_buf, '(A,3(ES13.6,1X))') &
            '   COB Offset    : ', body%prm%offset_cob
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') '   Mass Disp.    : ', body%prm%mass_disp
       call neko_log%message(log_buf)
    end if

    write(log_buf, '(A,3(ES13.6,1X))') '   Lin Stiffness : ', body%prm%K_lin
    call neko_log%message(log_buf)
    write(log_buf, '(A,3(ES13.6,1X))') '   Lin Damping   : ', body%prm%C_lin
    call neko_log%message(log_buf)
    write(log_buf, '(A,3(ES13.6,1X))') '   Rot Stiffness : ', body%prm%K_ang
    call neko_log%message(log_buf)
    write(log_buf, '(A,3(ES13.6,1X))') '   Rot Damping   : ', body%prm%C_ang
    call neko_log%message(log_buf)

    write(log_buf, '(A,6(ES13.6,1X))') '   Pos Equilib   : ', body%prm%pos_eq
    call neko_log%message(log_buf)
    write(log_buf, '(A,6(ES13.6,1X))') '   Init Velocity : ', body%initial_vel
    call neko_log%message(log_buf)

    write(log_buf, '(A,6(ES13.6,1X))') &
         '   Prescribed F0 : ', body%prm%F_prescribed_pivot
    call neko_log%message(log_buf)

    call neko_log%message(' ')
  end subroutine log_body

  !> Solve a small dense linear system by Gaussian elimination with partial
  !! pivoting. CPU only.
  subroutine linsolve_dense(n, A_in, b_in, x_out)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: A_in(n,n), b_in(n)
    real(kind=rp), intent(out) :: x_out(n)

    real(kind=rp) :: A(n,n), b(n)
    integer :: i, j, k, p
    real(kind=rp) :: factor, temp, pmax

    A = A_in
    b = b_in

    do k = 1, n-1
       p = k
       pmax = abs(A(k,k))
       do i = k+1, n
          if (abs(A(i,k)) > pmax) then
             pmax = abs(A(i,k))
             p = i
          end if
       end do

       if (p /= k) then
          do j = k, n
             temp = A(k,j)
             A(k,j) = A(p,j)
             A(p,j) = temp
          end do
          temp = b(k)
          b(k) = b(p)
          b(p) = temp
       end if

       do i = k+1, n
          if (A(i,k) /= 0.0_rp) then
             factor = A(i,k) / A(k,k)
             do j = k, n
                A(i,j) = A(i,j) - factor * A(k,j)
             end do
             b(i) = b(i) - factor * b(k)
          end if
       end do
    end do

    x_out(n) = b(n) / A(n,n)
    do i = n-1, 1, -1
       do j = i+1, n
          b(i) = b(i) - A(i,j)*x_out(j)
       end do
       x_out(i) = b(i) / A(i,i)
    end do
  end subroutine linsolve_dense

  !> Solve the structural system M * X = B.
  !! @details The non-linear correction and the user terms depend on X, so M
  !! and B are rebuilt from their linear base and the system is solved again
  !! until X stops changing. When nothing depends on X, the second pass
  !! reproduces the first one and the loop exits with a zero residual.
  !! @param X_sol The initial guess on input, the solution on output.
  !! @param verbose If true, log the residual of every pass.
  subroutine fsi_solve_structure(nbodies_fsi, bodies, fsi_dof_map, n, &
       M_global, B_global, X_sol, rot_matrices, time, gravity_vec, &
       non_linear_correction_term, user_terms, verbose)
    integer, intent(in) :: nbodies_fsi
    type(fsi_body_t), intent(in) :: bodies(:)
    integer, intent(in) :: fsi_dof_map(:,:)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: M_global(n,n), B_global(n), X_sol(n)
    real(kind=rp), intent(in) :: rot_matrices(:,:,:)
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: gravity_vec(3)
    logical, intent(in) :: non_linear_correction_term
    procedure(user_fsi_structural_terms_intf), pointer, intent(in) :: &
         user_terms
    logical, intent(in) :: verbose
    real(kind=rp), allocatable :: M_linear(:,:), B_linear(:), X_prev(:)
    real(kind=rp) :: nr_residual
    integer :: nr_iter
    integer, parameter :: max_nr_iter = 20
    logical :: converged
    character(len=1000) :: msg

    allocate(M_linear(n, n))
    allocate(B_linear(n))
    allocate(X_prev(n))

    M_linear = M_global
    B_linear = B_global

    if (verbose) then
       call neko_log%message("  --- Fixed-point Iteration ---")
    end if

    converged = .false.
    do nr_iter = 1, max_nr_iter
       ! Reset to linear base
       M_global = M_linear
       B_global = B_linear

       if (non_linear_correction_term) then
          call add_fsi_non_linear_matrices(nbodies_fsi, bodies, fsi_dof_map, &
               M_global, X_sol, rot_matrices)
       end if

       call add_fsi_user_structural_terms(nbodies_fsi, bodies, fsi_dof_map, &
            M_global, B_global, X_sol, rot_matrices, time, gravity_vec, &
            user_terms)

       X_prev = X_sol

       call linsolve_dense(n, M_global, B_global, X_sol)

       nr_residual = maxval(abs(X_sol - X_prev))

       if (verbose) then
          write(msg, '(A, I2, A, ES13.6)') "    Iter: ", &
               nr_iter, " | Max Residual: ", nr_residual
          call neko_log%message(trim(msg))
       end if
       if (nr_residual .lt. 1.0e-14_rp) then
          converged = .true.
          exit
       end if
    end do

    if (.not. converged) then
       write(msg, '(A,I0,A,I0,A,ES13.6)') &
            "FSI structural loop did not converge at step ", time%tstep, &
            " after ", max_nr_iter, " passes. Max residual: ", nr_residual
       call neko_log%warning(trim(msg))
    end if

    deallocate(M_linear)
    deallocate(B_linear)
    deallocate(X_prev)
  end subroutine fsi_solve_structure

  !> Register the FSI rigid-body state as the "fsi" checkpoint payload.
  !! Every rank holds an identical copy of these arrays.
  subroutine fsi_register_checkpoint(chkp, global_disp_rel, &
       global_body_vel, global_body_vel_lag, &
       global_moving_frame_presc_vel, global_body_acc, global_frame_acc)
    type(chkp_t), intent(inout) :: chkp
    real(kind=rp), target, intent(inout) :: global_disp_rel(:)
    real(kind=rp), target, intent(inout) :: global_body_vel(:)
    real(kind=rp), target, intent(inout) :: global_body_vel_lag(:,:)
    real(kind=rp), target, intent(inout) :: global_moving_frame_presc_vel(:,:)
    !> Newmark previous-acceleration
    real(kind=rp), target, intent(inout), optional :: global_body_acc(:)
    !> Newmark prescribed-frame previous-acceleration
    real(kind=rp), target, intent(inout), optional :: global_frame_acc(:)
    type(checkpoint_payload_t), pointer :: payload

    payload => chkp%add_payload("fsi")
    call payload%add_array("disp_rel", global_disp_rel, replicated = .true.)
    call payload%add_array("body_vel", global_body_vel, replicated = .true.)
    call payload%add_array("body_vel_lag", global_body_vel_lag, &
         replicated = .true.)
    call payload%add_array("moving_frame_presc_vel", &
         global_moving_frame_presc_vel, replicated = .true.)
    if (present(global_body_acc)) then
       call payload%add_array("body_acc", global_body_acc, &
            replicated = .true.)
    end if
    if (present(global_frame_acc)) then
       call payload%add_array("frame_acc", global_frame_acc, &
            replicated = .true.)
    end if
  end subroutine fsi_register_checkpoint

  !> Flattens FSI body arrays into global 1D/2D arrays for checkpointing
  subroutine fsi_prep_checkpoint(nbodies_fsi, bodies, global_disp_rel, &
       global_body_vel, global_body_vel_lag, &
       global_moving_frame_presc_vel, global_body_acc, global_frame_acc)
    integer, intent(in) :: nbodies_fsi
    type(fsi_body_t), intent(in) :: bodies(:)
    real(kind=rp), intent(inout) :: global_disp_rel(:)
    real(kind=rp), intent(inout) :: global_body_vel(:)
    real(kind=rp), intent(inout) :: global_body_vel_lag(:,:)
    real(kind=rp), intent(inout) :: global_moving_frame_presc_vel(:,:)
    !> Newmark previous-acceleration
    real(kind=rp), intent(inout), optional :: global_body_acc(:)
    !> Newmark prescribed-frame previous-acceleration
    real(kind=rp), intent(inout), optional :: global_frame_acc(:)

    integer :: i, idx_base

    if (nbodies_fsi == 0) return

    do i = 1, nbodies_fsi
       idx_base = (i - 1) * 6

       global_disp_rel(idx_base + 1 : idx_base + 6) = &
            bodies(i)%disp_rel(1:6)

       global_body_vel(idx_base + 1 : idx_base + 6) = &
            bodies(i)%body_vel(1:6)

       global_body_vel_lag(idx_base + 1 : idx_base + 6, :) = &
            bodies(i)%body_vel_lag(1:6, :)

       global_moving_frame_presc_vel(idx_base + 1 : idx_base + 6, :) = &
            bodies(i)%moving_frame_presc_vel(1:6, :)

       if (present(global_body_acc)) then
          global_body_acc(idx_base + 1 : idx_base + 6) = &
               bodies(i)%body_acc(1:6)
       end if

       if (present(global_frame_acc)) then
          global_frame_acc(idx_base + 1 : idx_base + 6) = &
               bodies(i)%moving_frame_presc_acc_prev(1:6)
       end if
    end do
  end subroutine fsi_prep_checkpoint

  !> Restores FSI body arrays from global arrays after restart read
  subroutine fsi_restart_restore(nbodies_fsi, bodies, global_disp_rel, &
       global_body_vel, global_body_vel_lag, &
       global_moving_frame_presc_vel, global_body_acc, global_frame_acc)
    integer, intent(in) :: nbodies_fsi
    type(fsi_body_t), intent(inout) :: bodies(:)
    real(kind=rp), intent(in) :: global_disp_rel(:)
    real(kind=rp), intent(in) :: global_body_vel(:)
    real(kind=rp), intent(in) :: global_body_vel_lag(:,:)
    real(kind=rp), intent(in) :: global_moving_frame_presc_vel(:,:)
    !> Newmark previous-acceleration
    real(kind=rp), intent(in), optional :: global_body_acc(:)
    !> Newmark prescribed-frame previous-acceleration
    real(kind=rp), intent(in), optional :: global_frame_acc(:)

    integer :: i, idx_base

    if (nbodies_fsi == 0) return

    do i = 1, nbodies_fsi
       idx_base = (i - 1) * 6

       bodies(i)%disp_rel(1:6) = &
            global_disp_rel(idx_base + 1 : idx_base + 6)
       bodies(i)%body_vel(1:6) = &
            global_body_vel(idx_base + 1 : idx_base + 6)

       bodies(i)%body_vel_lag(1:6, :) = &
            global_body_vel_lag(idx_base + 1 : idx_base + 6, :)

       bodies(i)%moving_frame_presc_vel(1:6, :) = &
            global_moving_frame_presc_vel(idx_base + 1 : idx_base + 6, :)

       if (present(global_body_acc)) then
          bodies(i)%body_acc(1:6) = &
               global_body_acc(idx_base + 1 : idx_base + 6)
       end if

       if (present(global_frame_acc)) then
          bodies(i)%moving_frame_presc_acc_prev(1:6) = &
               global_frame_acc(idx_base + 1 : idx_base + 6)
       end if
    end do
  end subroutine fsi_restart_restore
end module fsi_manager
