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
module fluid_pnpn_fsi_greens
  use fsi_dynamics, only : fsi_body_t, assemble_structural_inertial_terms, &
       add_fsi_non_linear_matrices, add_fsi_user_structural_terms
  use fsi_manager, only: fsi_manager_init, linsolve_dense, &
       fsi_prep_checkpoint, fsi_restart_restore
  use fluid_pnpn, only : fluid_pnpn_t
  use force_torque, only : force_torque_t
  use field, only : field_t
  use field_math, only : field_add2, field_copy, field_cmult, field_rzero, &
       field_add2s2, field_cfill
  use num_types, only : rp, dp
  use mathops, only : opadd2cm
  use device_mathops, only : device_opadd2cm
  use neko_config, only : NEKO_BCKND_DEVICE
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use projection, only : projection_t
  use projection_vel, only : projection_vel_t
  use bc, only : bc_t
  use no_slip, only : no_slip_t
  use inflow, only : inflow_t
  use json_module, only : json_file
  use utils, only : neko_error
  use logger, only : neko_log, LOG_SIZE
  use mesh, only : mesh_t
  use user_intf, only : user_t, user_fsi_body_params_intf, &
       user_fsi_structural_terms_intf, dummy_fsi_structural_terms
  use checkpoint, only : chkp_t
  use mpi_f08, only: MPI_Wtime
  use ab_time_scheme, only : ab_time_scheme_t
  use math, only : rzero
  use fld_file, only : fld_file_t
  use file, only : file_t
  use dofmap, only : dofmap_t
  use ale_manager
  use mxm_wrapper, only : mxm
  use tensor, only : tnsr3d, trsp
  implicit none
  private

  ! Green's method of Fischer, P., Schmitt, M., & Tomboulides, A. (2017). 
  ! Recent developments in spectral element simulations of moving-domain problems
  ! Recent progress and modern challenges in 
  ! applied mathematics, modeling and computational science, 213-244.
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
     procedure(user_fsi_body_params_intf), nopass, pointer :: &
          user_fsi_body_params => null()
     !> User hook: extra structural equation terms.
     procedure(user_fsi_structural_terms_intf), nopass, pointer :: &
          user_fsi_structural_terms => null()
     !> True when the user registered structural terms.
     logical :: has_user_structural_terms = .false.
   contains
     procedure, pass(this) :: init => fluid_fsi_init
     procedure, pass(this) :: step => fluid_fsi_step
     procedure, pass(this) :: free => fluid_fsi_free
     !> Restart from a previous solution.
     procedure, pass(this) :: restart => fluid_fsi_restart
     procedure, pass(this) :: calc_fsi_terms => &
           assemble_fsi_structural_inertial_terms
     procedure, pass(this) :: log_fsi_results => fluid_fsi_log_results
     procedure, pass(this) :: query_frame_prescribed_motion => &
          fluid_fsi_query_frame_prescribed_motion
  end type fluid_pnpn_fsi_greens_t

contains

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
    call fsi_manager_init(params, msh, this%ale, this%c_Xh, this%dm_Xh, &
         this%if_fsi, this%nbodies_fsi, this%fsi_bodies, this%fsi_dof_map, &
         this%total_active_dofs, &
         this%M_global, this%B_global, this%X_sol, &
         this%u_g, this%v_g, this%w_g, this%p_g, &
         this%res_long_print, this%gravity_vec, this%proj_prs_green, &
         this%proj_vel_green, this%global_disp_rel, &
         this%global_body_vel, this%global_body_vel_lag, &
         this%global_moving_frame_presc_vel, this%skip_greens_solve, &
         this%non_linear_correction_term)

    call this%chkp%add_fsi(this%global_disp_rel, this%global_body_vel, &
         this%global_body_vel_lag, &
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
            override_ang = this%batch_ang,&
            out_prescribed_vels = this%temp_prescribed_vels, &
            mode = 0)
       do i = 1, this%nbodies_fsi
          this%fsi_bodies(i)%moving_frame_presc_vel(:, 0) = &
               this%temp_prescribed_vels(:, i)
       end do
    end if

  end subroutine fluid_fsi_init


  subroutine fluid_fsi_step(this, time, dt_controller)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    type(ab_time_scheme_t) :: ab_scheme_obj
    character(len=1000) :: msg

    real(kind=rp) :: ab_coeffs(4), dt_history(10)
    real(kind=rp) :: beta(0:3)
    integer :: nadv, n, i, j, k, row_g, col_g, k_row, idx_g
    real(kind=rp) :: F_fluid(6)
    real(kind=dp) :: start_time_s, end_time_s, step_time_s
    real(kind=dp) :: start_time_g, end_time_g, step_time_g
    real(kind=dp), save :: total_elapsed_s = 0.0_dp
    real(kind=dp), save :: total_elapsed_g = 0.0_dp

    ! Fixed-Point Iteration variables
    real(kind=rp), allocatable :: M_linear(:,:), B_linear(:), X_prev(:)
    real(kind=rp) :: nr_residual
    integer :: nr_iter, max_nr_iter
    logical :: iter_verbose, converged

    n = this%dm_Xh%size()
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

    !write(msg, '(A, 4(F10.5, 1X))') "DEBUG BETA (0-3):   ", &
    !      beta(0), beta(1), beta(2), beta(3)
    !call neko_log%message(trim(msg))

    !write(msg, '(A, 3(F10.5, 1X))') "DEBUG ALPHA_EXT (1-3):", &
    !      alpha_ext(1), alpha_ext(2), alpha_ext(3)
    !call neko_log%message(trim(msg))

    !write(msg, '(A, 4(F10.5, 1X))') "DEBUG AB_COEFFS (1-3):", &
    !      ab_coeffs(1), ab_coeffs(2), ab_coeffs(3), ab_coeffs(4)
    !call neko_log%message(trim(msg))


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
    ! I think I can move compute_rotationm matrics inside advance_mesh_explicit. 
    ! I need to remember if there was any reason that 
    ! I put it there at first place.
    ! Calling update_mesh_velocity here should be totally harmless.
    call this%step_ext(time, dt_controller, greens_function = .false., &
         skip_ale_msh_vel_update = .false.)

    ! Add FSI structural terms at current time step.
    ! Rotation matrix etc should be updated at this point.
    ! I moved this after the above step_ext, since the rotation matrix should for 
    ! the current time step. 
    ! Need to verify more.
    call this%calc_fsi_terms(time, beta, nadv)

    end_time_s = MPI_WTIME()
    step_time_s = end_time_s - start_time_s
    total_elapsed_s = total_elapsed_s + step_time_s

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

    ! Skip Green's function if skip_greens_solve is true
    if (.not. this%skip_greens_solve) then
       ! Green's Function Loop
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
                ! transaltional DOF
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
                  mode = 1)

             start_time_g = MPI_WTIME()
             call this%step_ext(time, dt_controller, greens_function = .true., &
                  skip_ale_msh_vel_update = .true., &
                  proj_prs_green = this%proj_prs_green(col_g), &
                  proj_vel_green = this%proj_vel_green(col_g))

             end_time_g = MPI_WTIME()
             step_time_g = end_time_g - start_time_g
             total_elapsed_g = total_elapsed_g + step_time_g

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
    else
       call neko_log%message("Weak coupling enabled: " // &
            "Skipping Green's function fluid feedback.")
       step_time_g = 0.0_dp
    end if

    call neko_log%message(' ')
    write(msg, '(A, E15.7, A, I0, A, E15.7)') &
         "Standard's step total elapsed time (s):  ", &
         total_elapsed_s, "  Step: ", time%tstep, "  time: ", time%t
    call neko_log%message(trim(msg))
    write(msg, '(A, E15.7, A, I0, A, E15.7)') &
         "Green's step total elapsed time (s):  ", &
         total_elapsed_g, "  Step: ", time%tstep, "  time: ", time%t
    call neko_log%message(trim(msg))
    call neko_log%message(' ')

    ! Calculate all FSI corrections: M_global * X_sol = B_global.
    ! With current algorithm, the correction is actually just
    ! the velocity change compared to the previous time step.
    ! X_sol enters holding the previous step's solution, which serves as
    ! the initial guess of the fixed-point loop below.
    if (this%total_active_dofs > 0) then
       allocate(M_linear(this%total_active_dofs, this%total_active_dofs))
       allocate(B_linear(this%total_active_dofs))
       allocate(X_prev(this%total_active_dofs))

       ! Save the linear base. Both M and B are rebuilt every pass: the
       ! built-in nonlinear correction touches M, user terms touch both.
       M_linear = this%M_global
       B_linear = this%B_global

       ! The loop always runs to convergence. When nothing depends on
       ! X_sol, pass 2 reproduces pass 1 (linsolve_dense is
       ! deterministic), the residual is exactly zero, and the result is
       ! identical to a single solve. X_sol enters with the previous
       ! step's solution as a warm start.
       max_nr_iter = 20
       iter_verbose = this%non_linear_correction_term .or. &
            this%has_user_structural_terms
       if (iter_verbose) then
          call neko_log%message("  --- Fixed-point Iteration ---")
       end if

       converged = .false.
       do nr_iter = 1, max_nr_iter
          ! Reset to linear base
          this%M_global = M_linear
          this%B_global = B_linear

          if (this%non_linear_correction_term) then
             call add_fsi_non_linear_matrices(this%nbodies_fsi, &
                  this%fsi_bodies, this%fsi_dof_map, this%M_global, &
                  this%X_sol, this%ale%body_rot_matrices)
          end if

          call add_fsi_user_structural_terms(this%nbodies_fsi, &
               this%fsi_bodies, this%fsi_dof_map, this%M_global, &
               this%B_global, this%X_sol, this%ale%body_rot_matrices, &
               time, this%gravity_vec, this%user_fsi_structural_terms)

          X_prev = this%X_sol

          call linsolve_dense(this%total_active_dofs, this%M_global, &
               this%B_global, this%X_sol)

          nr_residual = maxval(abs(this%X_sol - X_prev))

          if (iter_verbose) then
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
       if (iter_verbose) call neko_log%message(' ')

       deallocate(M_linear)
       deallocate(B_linear)
       deallocate(X_prev)
    end if

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

    ! Structure Update
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

    ! Calculate Final Velocity (FSI + Prescribed)
    ! This velocity will be used as the "guessed" velocity
    ! for the next time step, and also for the ALE mesh update.
    call this%ale%update_mesh_velocity(this%c_Xh, time, &
          override_ids = this%batch_ids, &
          override_trans = this%batch_trans, &
          override_ang = this%batch_ang,&
          out_prescribed_vels = this%temp_prescribed_vels, &
          mode = 0)
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

  subroutine fluid_fsi_free(this)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
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

    ! Free the base fluid_pnpn_t
    call this%fluid_pnpn_t%free()

  end subroutine fluid_fsi_free

  subroutine fluid_fsi_query_frame_prescribed_motion(this, time, beta, nadv)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    real(kind=rp) :: gamma
    integer, intent(in) :: nadv
    integer :: i, k

    gamma = beta(0) / time%dt

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
         mode = 2)

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
               (beta(k) * this%fsi_bodies(i)%moving_frame_presc_vel(:, k)) / time%dt
       end do
    end do
  end subroutine fluid_fsi_query_frame_prescribed_motion

  subroutine assemble_fsi_structural_inertial_terms(this, time, beta, nadv)
    class(fluid_pnpn_fsi_greens_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: beta(0:3)
    real(kind=rp) :: gamma
    integer, intent(in) :: nadv
    integer :: i, k

    gamma = beta(0) / time%dt

    ! Here we fill the M_global and B_global
    ! using the contributuon from strucutre and also
    ! bodies' inertial motion from previous time steps.
    call assemble_structural_inertial_terms(this%nbodies_fsi, &
         this%fsi_bodies, &
         this%fsi_dof_map, this%M_global, this%B_global, &
         this%ale%body_rot_matrices, &
         time, gamma, beta, nadv, this%gravity_vec, &
         this%user_fsi_body_params)

  end subroutine assemble_fsi_structural_inertial_terms

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

    call neko_log%message("variable, time step, time, body, x_val, y_val, z_val")

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

  subroutine fluid_fsi_restart(this, chkp)
    class(fluid_pnpn_fsi_greens_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    type(time_state_t) :: t_restart
    real(kind=rp) :: dtlag(10), tlag(10)
    integer :: i, n

    dtlag = chkp%dtlag
    tlag = chkp%tlag

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
       t_restart%dt = chkp%dtlag(1)

       do i = 1, this%nbodies_fsi
          this%batch_ids(i) = this%fsi_bodies(i)%ale_id
          this%batch_trans(:, i) = this%fsi_bodies(i)%body_vel(1:3)
          this%batch_ang(:, i) = this%fsi_bodies(i)%body_vel(4:6)
       end do

       call this%ale%update_mesh_velocity(this%c_Xh, t_restart, &
            override_ids = this%batch_ids, &
            override_trans = this%batch_trans, &
            override_ang = this%batch_ang, &
            mode = 0)
    end if
  end subroutine fluid_fsi_restart

end module fluid_pnpn_fsi_greens