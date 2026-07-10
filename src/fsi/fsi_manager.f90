module fsi_manager
  use fsi_dynamics, only : fsi_body_t, apply_parallel_axis_theorem
  use fluid_pnpn, only : fluid_pnpn_t
  use force_torque, only : force_torque_t
  use field, only : field_t
  use field_math, only : field_rzero
  use num_types, only : rp, dp
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_lookup, json_get_or_default, &
       json_extract_item
  use utils, only : neko_error
  use logger, only : neko_log, LOG_SIZE
  use mesh, only : mesh_t
  use user_intf, only : user_t
  use checkpoint, only : chkp_t
  use ale_manager, only : ale_manager_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use projection, only : projection_t
  use projection_vel, only : projection_vel_t
  implicit none
  private

  public :: fsi_manager_init
  public :: linsolve_dense
  public :: fsi_prep_checkpoint
  public :: fsi_restart_restore
contains

  subroutine fsi_manager_init(params, msh, ale, c_Xh, dm_Xh, &
        if_fsi, nbodies_fsi, bodies, fsi_dof_map, &
        total_active_dofs, M_global, B_global, X_sol, &
        u_g, v_g, w_g, p_g, res_long_print, gravity_vec, &
        proj_prs_green, proj_vel_green, global_disp_rel, &
        global_body_vel, global_body_vel_lag, &
        global_moving_frame_presc_vel, skip_greens_solve, &
        non_linear_correction_term, global_body_acc, global_frame_acc)

    ! Inputs needed for setup
    type(json_file), target, intent(inout) :: params
    type(mesh_t), target, intent(inout) :: msh
    type(ale_manager_t), target, intent(inout) :: ale
    type(coef_t), intent(inout) :: c_Xh
    type(dofmap_t), intent(in) :: dm_Xh

    ! FSI data structures to be initialized
    logical, intent(out) :: if_fsi
    integer, intent(out) :: nbodies_fsi
    type(fsi_body_t), allocatable, intent(inout) :: bodies(:)
    integer, allocatable, intent(inout) :: fsi_dof_map(:,:)
    integer, intent(out) :: total_active_dofs
    real(kind=rp), allocatable, intent(inout) :: M_global(:,:), B_global(:)
    real(kind=rp), allocatable, intent(inout) :: X_sol(:)
    type(field_t), allocatable, intent(inout) :: u_g(:), v_g(:), w_g(:)
    type(field_t), allocatable, intent(inout) :: p_g(:)
    type(projection_t), allocatable, intent(out) :: proj_prs_green(:)
    type(projection_vel_t), allocatable, intent(out) :: proj_vel_green(:)
    integer :: fsi_pr_projection_dim !< Steps to activate projection for ksp_vel
    integer :: fsi_vel_projection_dim !< Steps to activate projection for ksp_vel

    integer :: fsi_vel_projection_activ_step !< Steps to activate projection for ksp_pr
    integer :: fsi_pr_projection_activ_step !< Steps to activate projection for ksp_pr for FSI problem (green's solve)
    logical :: fsi_pr_projection_reorthogonalize_basis
    real(kind=rp), intent(out) :: gravity_vec(3)
    logical, intent(out) :: res_long_print
    logical, intent(out) :: skip_greens_solve
    logical, intent(out) :: non_linear_correction_term

    type(json_file) :: body_sub
    integer :: i, j, k, m, n_bodies
    character(len=200) :: name_buf, name_buf2, log_buf, field_name
    real(kind=rp) :: center_dummy(3), force_scale
    real(kind=rp) :: pivot_init(3)
    real(kind=rp) :: I_pivot_temp(3,3), I_com_temp(3,3)
    real(kind=rp), allocatable :: temp_vec(:)
    character(len=:), allocatable :: temp_str
    integer, allocatable :: temp_vec_int(:)
    character(len=20) :: coupling_mode
    integer :: temp_vec_int6(6)
    logical :: log_forces, long_print

    real(kind=rp), allocatable, intent(inout) :: global_disp_rel(:)
    real(kind=rp), allocatable, intent(inout) :: global_body_vel(:)
    real(kind=rp), allocatable, intent(inout) :: global_body_vel_lag(:,:)
    real(kind=rp), allocatable, intent(inout) :: global_moving_frame_presc_vel(:,:)
    !> Newmark previous-acceleration.
    real(kind=rp), allocatable, intent(inout), optional :: global_body_acc(:)
    !> Newmark prescribed-frame previous-acceleration.
    real(kind=rp), allocatable, intent(inout), optional :: global_frame_acc(:)

    center_dummy = 0.0_rp

    call neko_log%section("Fluid-Structure Interaction")

    if (params%valid_path('case.fluid.fsi')) then
       call json_get(params, 'case.fluid.fsi.enabled', if_fsi)
       if (.not. if_fsi) then
          call neko_error("Scheme pnpn_fsi: FSI block present, but 'fsi.enabled': false in case file")
       end if

       call json_get_or_default(params, 'case.fluid.fsi.log_forces', log_forces, .true.)
       call json_get_or_default(params, 'case.fluid.fsi.long_print', long_print, .false.)
       call json_get_or_default(params, 'case.fluid.fsi.results_long_print', res_long_print, .false.)
       call json_get_or_default(params, 'force_scale', force_scale, 1.0_rp)
       call json_get_or_default(params, 'case.fluid.fsi.skip_greens_solve', skip_greens_solve, .false.)
       call json_get_or_default(params, 'case.fluid.fsi.non_linear_correction_term', non_linear_correction_term, .false.)

       call json_get(params, 'case.fluid.fsi.coupling', temp_str)
       if (trim(temp_str) .eq. 'subiteration') then
         skip_greens_solve = .true.
         coupling_mode = "subiteration"
       elseif (trim(temp_str) .eq. 'greens') then
         coupling_mode = "greens"
       else
         call neko_error("FSI: Unknown coupling mode: " // trim(temp_str) // &
                ". Must be 'subiteration' or 'greens'.")
       end if
        
       if ( (.not. skip_greens_solve) .and. coupling_mode == "greens" ) then
          call neko_log%message(" ")
          call neko_log%message("-----------------------------------------------------")
          call neko_log%message("FSI Coupling: Greens + Strong Coupling (Implicit FSI)")
          call neko_log%message("-----------------------------------------------------")
          call neko_log%message(" ")
       else if ( (skip_greens_solve) .and. coupling_mode == "greens" ) then
          call neko_log%message(" ")
          call neko_log%message("---------------------------------------------------")
          call neko_log%message("FSI Coupling: Greens + Weak Coupling (Explicit FSI)")
          call neko_log%message("---------------------------------------------------")
          call neko_log%message(" ")
       else if (coupling_mode == "subiteration") then
          call neko_log%message(" ")
          call neko_log%message("-----------------------------------------------------")
          call neko_log%message("FSI Coupling: Subiteration")
          call neko_log%message("-----------------------------------------------------")
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
              fsi_pr_projection_reorthogonalize_basis, .false.)

       call json_get_or_default(params, &
              'case.fluid.fsi.velocity_solver.projection_space_size', &
              fsi_vel_projection_dim, 0)

       call json_get_or_default(params, &
              'case.fluid.fsi.velocity_solver.projection_hold_steps', &
              fsi_vel_projection_activ_step, 5)

       gravity_vec = 0.0_rp
       if (params%valid_path('case.fluid.fsi.gravity_vec')) then
          call json_get(params, 'case.fluid.fsi.gravity_vec', temp_vec, expected_size = 3)
          gravity_vec = temp_vec
       end if

       call params%info('case.fluid.fsi.bodies', n_children = n_bodies)
       nbodies_fsi = n_bodies
       if (nbodies_fsi == 0) then
          call neko_error("FSI: 'case.fluid.fsi.bodies' is empty (no FSI bodies defined)")
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

       if (allocated(bodies)) deallocate(bodies)
       allocate(bodies(nbodies_fsi))

       if (allocated(fsi_dof_map)) deallocate(fsi_dof_map)
       allocate(fsi_dof_map(nbodies_fsi, 6))

       if (allocated(global_disp_rel)) deallocate(global_disp_rel)
       allocate(global_disp_rel(6 * nbodies_fsi))

       if (allocated(global_body_vel)) deallocate(global_body_vel)
       allocate(global_body_vel(6 * nbodies_fsi))

       if (present(global_body_acc)) then
          if (allocated(global_body_acc)) deallocate(global_body_acc)
          allocate(global_body_acc(6 * nbodies_fsi))
          global_body_acc = 0.0_rp
       end if

       if (present(global_frame_acc)) then
          if (allocated(global_frame_acc)) deallocate(global_frame_acc)
          allocate(global_frame_acc(6 * nbodies_fsi))
          global_frame_acc = 0.0_rp
       end if

       if (allocated(global_body_vel_lag)) deallocate(global_body_vel_lag)
       allocate(global_body_vel_lag(6 * nbodies_fsi, &
              lbound(bodies(1)%body_vel_lag, 2) : ubound(bodies(1)%body_vel_lag, 2)))

       if (allocated(global_moving_frame_presc_vel)) deallocate(global_moving_frame_presc_vel)
       allocate(global_moving_frame_presc_vel(6 * nbodies_fsi, &
              lbound(bodies(1)%moving_frame_presc_vel, 2) : ubound(bodies(1)%moving_frame_presc_vel, 2)))

       global_disp_rel = 0.0_rp
       global_body_vel = 0.0_rp
       global_body_vel_lag = 0.0_rp
       global_moving_frame_presc_vel = 0.0_rp

       total_active_dofs = 0

       do i = 1, nbodies_fsi
          call json_extract_item(params, 'case.fluid.fsi.bodies', i, body_sub)

          if (body_sub%valid_path('name')) then
             call json_get(body_sub, 'name', temp_str)
             bodies(i)%name = temp_str
          else
             write(bodies(i)%name, '(A,I0)') 'fsi_body_', i
          endif

          call json_get(body_sub, 'zone_id', bodies(i)%zone_id)

          bodies(i)%ale_id = -1
          do j = 1, ale%config%nbodies
             if (any(ale%config%bodies(j)%zone_indices == bodies(i)%zone_id)) then
                bodies(i)%ale_id = j
                exit
             end if
          end do

          if (bodies(i)%ale_id == -1) then
             call neko_error("FSI: Body " // trim(bodies(i)%name) // " zone_id not found in ALE bodies!")
          end if

          call json_get(body_sub, 'mass', bodies(i)%mass)

          bodies(i)%I_body_tensor = 0.0_rp
          if (body_sub%valid_path('I_xx_xy_xz_p') .or. &
               body_sub%valid_path('I_yx_yy_yz_p') .or. &
               body_sub%valid_path('I_zx_zy_zz_p')) then

             call json_get(body_sub, 'I_xx_xy_xz_p', temp_vec, expected_size = 3)
             bodies(i)%I_body_tensor(1,:) = temp_vec

             call json_get(body_sub, 'I_yx_yy_yz_p', temp_vec, expected_size = 3)
             bodies(i)%I_body_tensor(2,:) = temp_vec

             call json_get(body_sub, 'I_zx_zy_zz_p', temp_vec, expected_size = 3)
             bodies(i)%I_body_tensor(3,:) = temp_vec

             bodies(i)%inertia_ref_frame = 'pivot'

          else if (body_sub%valid_path('I_xx_xy_xz_com') .or. &
               body_sub%valid_path('I_yx_yy_yz_com') .or. &
               body_sub%valid_path('I_zx_zy_zz_com')) then

             call json_get(body_sub, 'I_xx_xy_xz_com', temp_vec, expected_size = 3)
             bodies(i)%I_body_tensor(1,:) = temp_vec

             call json_get(body_sub, 'I_yx_yy_yz_com', temp_vec, expected_size = 3)
             bodies(i)%I_body_tensor(2,:) = temp_vec

             call json_get(body_sub, 'I_zx_zy_zz_com', temp_vec, expected_size = 3)
             bodies(i)%I_body_tensor(3,:) = temp_vec

             bodies(i)%inertia_ref_frame = 'com'

          else
             call neko_error("FSI: Body " // trim(bodies(i)%name) // &
                    " missing full inertia tensor inputs (I_xx_xy_xz_p OR I_xx_xy_xz_com)")
          end if

          call json_get(body_sub, 'active_dofs', temp_vec_int, expected_size = 6)
          bodies(i)%active_dofs = temp_vec_int


          call json_get(body_sub, 'center_of_mass', temp_vec, expected_size = 3)
          bodies(i)%center_of_mass = temp_vec

          pivot_init = ale%ale_pivot(bodies(i)%ale_id)%pos

          bodies(i)%local_offset_com = bodies(i)%center_of_mass - &
                 pivot_init

          if (trim(bodies(i)%inertia_ref_frame) == 'com') then
             I_com_temp = bodies(i)%I_body_tensor
             call apply_parallel_axis_theorem(bodies(i)%I_body_tensor,&
                    bodies(i)%mass, &
                    bodies(i)%local_offset_com, &
                    I_pivot_temp)
             bodies(i)%I_body_tensor = I_pivot_temp
             bodies(i)%inertia_ref_frame = 'pivot (shifted from com)'
          end if

          bodies(i)%center_of_buoyancy = 0.0_rp
          bodies(i)%local_offset_cob = 0.0_rp
          bodies(i)%mass_disp = 0.0_rp
          if (params%valid_path('case.fluid.fsi.gravity_vec')) then
             ! IF Gravity exists, buoyancy is mandatory;
             ! otherwise, the physics is wrong!
             call json_get(body_sub, 'center_of_buoyancy', temp_vec, expected_size = 3)
             bodies(i)%center_of_buoyancy = temp_vec

             bodies(i)%local_offset_cob = bodies(i)%center_of_buoyancy - &
                    pivot_init
             call json_get(body_sub, 'mass_disp', bodies(i)%mass_disp)
          end if

          ! Prescribed constant forcing at pivot location
          bodies(i)%F_prescribed_pivot = 0.0_rp
          if (body_sub%valid_path('F0_p')) then
             call json_get(body_sub, 'F0_p', temp_vec, &
                    expected_size = 6)
             bodies(i)%F_prescribed_pivot = temp_vec
          end if

          call json_get(body_sub, 'stiffness', temp_vec, expected_size = 3)
          bodies(i)%K_lin = temp_vec


          call json_get(body_sub, 'damping', temp_vec, expected_size = 3)
          bodies(i)%C_lin = temp_vec


          call json_get(body_sub, 'rot_stiffness', temp_vec, expected_size = 3)
          bodies(i)%K_ang = temp_vec


          call json_get(body_sub, 'rot_damping', temp_vec, expected_size = 3)
          bodies(i)%C_ang = temp_vec


          call json_get(body_sub, 'pos_equilibrium', temp_vec, expected_size = 6)
          bodies(i)%pos_eq = temp_vec

          call json_get(body_sub, 'initial_velocity', temp_vec, expected_size = 6)
          bodies(i)%initial_vel = temp_vec

          ! To be sure inital vel is 0 for inactive DOFs.
          do k = 1, 6
             if (bodies(i)%active_dofs(k) .eq. 0) then
                bodies(i)%initial_vel(k) = 0.0_rp
             end if
          end do


          bodies(i)%disp_rel = 0.0_rp
          bodies(i)%body_vel_lag = 0.0_rp
          bodies(i)%body_vel = bodies(i)%initial_vel

          ! I think it should be just 3. For the first two time
          ! lag(2) and lag(3) does not matter since we do not use them.
          ! but it's important to keep the in lag array since
          ! we do not touch lag(2) and lag(3) in updating lag arrays
          ! in first two time steps, so we will loose this if we don't
          ! fill them now.
          do m = 1, 3
             bodies(i)%body_vel_lag(:,m) = bodies(i)%initial_vel
          end do

          do k = 1, 6
             if (bodies(i)%active_dofs(k) == 1) then
                total_active_dofs = total_active_dofs + 1
                fsi_dof_map(i, k) = total_active_dofs
             else
                fsi_dof_map(i, k) = 0
             end if
          end do

          ! Comprehensive Per-Body Logging
          call neko_log%message('Registered Body : ' // trim(bodies(i)%name))

          write(log_buf, '(A,I0,A,I0)') '   Zone ID       : ', bodies(i)%zone_id, ' | ALE ID : ', bodies(i)%ale_id
          call neko_log%message(log_buf)

          write(log_buf, '(A,6(I2,1X))') '   Active DOFs   : ', bodies(i)%active_dofs
          call neko_log%message(log_buf)

          write(log_buf, '(A,ES18.11)') '   Mass          : ', bodies(i)%mass
          call neko_log%message(log_buf)

          if (trim(bodies(i)%inertia_ref_frame) == 'pivot (shifted from com)') then
             call neko_log%message('   Inertia Tensor (Input at COM) :')
             write(log_buf, '(A,3(ES13.6,1X))') '     I_xx_xy_xz  : ', I_com_temp(1,:)
             call neko_log%message(trim(log_buf))
             write(log_buf, '(A,3(ES13.6,1X))') '     I_yx_yy_yz  : ', I_com_temp(2,:)
             call neko_log%message(trim(log_buf))
             write(log_buf, '(A,3(ES13.6,1X))') '     I_zx_zy_zz  : ', I_com_temp(3,:)
             call neko_log%message(trim(log_buf))

             call neko_log%message('   Inertia Tensor (Shifted to Pivot) :')
             write(log_buf, '(A,3(ES13.6,1X))') '     I_xx_xy_xz  : ', bodies(i)%I_body_tensor(1,:)
             call neko_log%message(trim(log_buf))
             write(log_buf, '(A,3(ES13.6,1X))') '     I_yx_yy_yz  : ', bodies(i)%I_body_tensor(2,:)
             call neko_log%message(trim(log_buf))
             write(log_buf, '(A,3(ES13.6,1X))') '     I_zx_zy_zz  : ', bodies(i)%I_body_tensor(3,:)
             call neko_log%message(trim(log_buf))
          else
             call neko_log%message('   Inertia Tensor (Input at Pivot) :')
             write(log_buf, '(A,3(ES13.6,1X))') '     I_xx_xy_xz  : ', bodies(i)%I_body_tensor(1,:)
             call neko_log%message(trim(log_buf))
             write(log_buf, '(A,3(ES13.6,1X))') '     I_yx_yy_yz  : ', bodies(i)%I_body_tensor(2,:)
             call neko_log%message(trim(log_buf))
             write(log_buf, '(A,3(ES13.6,1X))') '     I_zx_zy_zz  : ', bodies(i)%I_body_tensor(3,:)
             call neko_log%message(trim(log_buf))
          end if

          write(log_buf, '(A,3(ES13.6,1X))') '   COM Position  : ', bodies(i)%center_of_mass
          call neko_log%message(log_buf)
          write(log_buf, '(A,3(ES13.6,1X))') '   COM Offset    : ', bodies(i)%local_offset_com
          call neko_log%message(log_buf)

          if (any(abs(gravity_vec) > 0.0_rp)) then
             write(log_buf, '(A,3(ES13.6,1X))') '   COB Position  : ', bodies(i)%center_of_buoyancy
             call neko_log%message(log_buf)
             write(log_buf, '(A,3(ES13.6,1X))') '   COB Offset    : ', bodies(i)%local_offset_cob
             call neko_log%message(log_buf)
             write(log_buf, '(A,ES13.6)') '   Mass Disp.    : ', bodies(i)%mass_disp
             call neko_log%message(log_buf)
          end if

          write(log_buf, '(A,3(ES13.6,1X))') '   Lin Stiffness : ', bodies(i)%K_lin
          call neko_log%message(log_buf)
          write(log_buf, '(A,3(ES13.6,1X))') '   Lin Damping   : ', bodies(i)%C_lin
          call neko_log%message(log_buf)
          write(log_buf, '(A,3(ES13.6,1X))') '   Rot Stiffness : ', bodies(i)%K_ang
          call neko_log%message(log_buf)
          write(log_buf, '(A,3(ES13.6,1X))') '   Rot Damping   : ', bodies(i)%C_ang
          call neko_log%message(log_buf)

          write(log_buf, '(A,6(ES13.6,1X))') '   Pos Equilib   : ', bodies(i)%pos_eq
          call neko_log%message(log_buf)
          write(log_buf, '(A,6(ES13.6,1X))') '   Init Velocity : ', bodies(i)%initial_vel
          call neko_log%message(log_buf)

          write(log_buf, '(A,6(ES13.6,1X))') '   Prescribed F0 : ', bodies(i)%F_prescribed_pivot
          call neko_log%message(log_buf)

          call neko_log%message(' ')

          write(name_buf, '(A,I0)') "fsi_force_monitor_", i
          write(name_buf2, '(A,A)') "fsi_body: ", trim(bodies(i)%name)

          call bodies(i)%force_monitor%init_common(name=trim(name_buf), &
                 fluid_name="fluid", &
                 zone_id=bodies(i)%zone_id, &
                 zone_name=name_buf2, &
                 center=center_dummy, &
                 scale=force_scale, &
                 coef=c_Xh, &
                 long_print=long_print, &
                 center_type='pivot', &
                 msh=msh, &
                 full_log=log_forces)
       end do

       if (total_active_dofs > 0) then
          if (allocated(M_global)) deallocate(M_global)
          allocate(M_global(total_active_dofs, total_active_dofs))

          if (allocated(B_global)) deallocate(B_global)
          allocate(B_global(total_active_dofs))

          if (allocated(X_sol)) deallocate(X_sol)
          allocate(X_sol(total_active_dofs))
          X_sol = 0.0_rp

          if (allocated(u_g)) deallocate(u_g)
          if (allocated(v_g)) deallocate(v_g)
          if (allocated(w_g)) deallocate(w_g)
          if (allocated(p_g)) deallocate(p_g)

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
                      fsi_pr_projection_activ_step, fsi_pr_projection_reorthogonalize_basis)

             call proj_vel_green(k)%init(dm_Xh%size(), fsi_vel_projection_dim, &
                      fsi_vel_projection_activ_step)
             call field_rzero(u_g(k))
             call field_rzero(v_g(k))
             call field_rzero(w_g(k))
             call field_rzero(p_g(k))
          end do

       else
          call neko_error("FSI: No active DOFs found!")
       end if

    else
       call neko_error("Missing 'case.fluid.fsi' block in case.json")
    end if

    call neko_log%end_section()

  end subroutine fsi_manager_init

  ! Gaussian elimination with partial pivoting
  ! CPU-only -> For small matrices
  subroutine linsolve_dense(n, A_in, b_in, x_out)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: A_in(n,n), b_in(n)
    real(kind=rp), intent(out) :: x_out(n)

    real(kind=rp) :: A(n,n), b(n)
    integer :: i, j, k, p
    real(kind=rp) :: factor, temp, pmax

    A = A_in; b = b_in

    do k = 1, n-1
       p = k
       pmax = abs(A(k,k))
       do i = k+1, n
          if (abs(A(i,k)) > pmax) then
             pmax = abs(A(i,k)); p = i
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
