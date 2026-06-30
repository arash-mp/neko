module fsi_dynamics
  use force_torque, only : force_torque_t
  use num_types, only : rp
  use time_state, only : time_state_t
  use logger, only : neko_log
  implicit none
  private

  public :: assemble_structural_inertial_terms
  public :: apply_parallel_axis_theorem
  public :: add_fsi_non_linear_matrices

  !> single FSI body properties
  type, public :: fsi_body_t
     character(len=256) :: name
     integer :: zone_id
     integer :: ale_id

     ! Active DOF Control
     ! 0: Deactive / 1: Active
     ! [u, v, w, omega_x, omega_y, omega_z]
     integer :: active_dofs(6)

     ! Structural Properties
     real(kind=rp) :: mass
     real(kind=rp) :: mass_disp
     real(kind=rp) :: I_body_diag(3) = 0.0_rp
     real(kind=rp) :: I_body_tensor(3,3) = 0.0_rp
     character(len=50) :: inertia_ref_frame = 'pivot'
     real(kind=rp) :: C_lin(3), C_ang(3)
     real(kind=rp) :: K_lin(3), K_ang(3)
     real(kind=rp) :: pos_eq(6)
     real(kind=rp) :: initial_vel(6) = 0.0_rp
     real(kind=rp) :: F_prescribed_pivot(6) = 0.0_rp

     ! Geometry
     real(kind=rp) :: local_offset_com(3)
     real(kind=rp) :: center_of_mass(3)
     real(kind=rp) :: local_offset_cob(3)
     real(kind=rp) :: center_of_buoyancy(3)

     ! State History
     real(kind=rp) :: disp_rel(6)
     real(kind=rp) :: body_vel(6)
     real(kind=rp) :: body_vel_guess(6)
     real(kind=rp) :: body_vel_lag(6, 3)

     ! Stores the prescribed velocity of the moving frame
     real(kind=rp) :: moving_frame_presc_acc(6) = 0.0_rp
     real(kind=rp) :: moving_frame_presc_vel(6, 0:3) = 0.0_rp

     type(force_torque_t) :: force_monitor
  end type fsi_body_t

contains

  subroutine assemble_structural_inertial_terms(nbodies_fsi, bodies, &
       fsi_dof_map, M_global, B_global, rot_matrices, &
       time, gamma, beta, nadv, gravity_vec)

    integer, intent(in) :: nbodies_fsi
    type(fsi_body_t), intent(inout) :: bodies(:)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: fsi_dof_map(:,:)
    real(kind=rp), intent(inout) :: M_global(:,:), B_global(:)
    real(kind=rp), intent(in) :: rot_matrices(:,:,:)
    real(kind=rp), intent(in) :: gravity_vec(3)
    real(kind=rp) :: dt, gamma, beta(0:3)
    integer, intent(in) :: nadv

    integer :: i, j, k, row_g, col_g
    real(kind=rp) :: m, m_disp
    real(kind=rp) :: R_mat(3,3), R_T(3,3), I_body(3,3), I_P(3,3)
    real(kind=rp) :: c(3), r_rel(3), r_cb(3)
    real(kind=rp) :: v_s(3), w_s(3)
    real(kind=rp) :: a_f(3), alpha_f(3), w_f(3)
    real(kind=rp) :: a_rel_s(3), alpha_rel_s(3), V_hist(6)

    real(kind=rp) :: M_local(6,6), B_local(6)
    real(kind=rp) :: C_skew(3,3), Wf_skew(3,3), Ws_skew(3,3), I3(3,3)
!    character(len=2048) :: msg

    dt = time%dt
    M_global = 0.0_rp
    B_global = 0.0_rp

    ! 3x3 Identity Matrix
    I3 = 0.0_rp
    I3(1,1) = 1.0_rp; I3(2,2) = 1.0_rp; I3(3,3) = 1.0_rp

    do i = 1, nbodies_fsi

       ! We use final corrected FSI (relative to the frame)
       ! as the gussed velocity for current time step.
       ! Contribution of frame movement is already included in
       ! in mesh velocity arrays.
       bodies(i)%body_vel_guess = bodies(i)%body_vel

       m = bodies(i)%mass
       m_disp = bodies(i)%mass_disp
       I_body = bodies(i)%I_body_tensor

       !> Frame States
       ! Translational acceleration
       a_f = bodies(i)%moving_frame_presc_acc(1:3)
       ! Angular acceleration
       alpha_f = bodies(i)%moving_frame_presc_acc(4:6)
       ! Angular velocity
       w_f = bodies(i)%moving_frame_presc_vel(4:6, 0)

       ! FSI (relative to the frame) displacements
       r_rel = bodies(i)%disp_rel(1:3)
       ! translational velocity
       v_s = bodies(i)%body_vel_guess(1:3)
       ! Angular velocity
       w_s = bodies(i)%body_vel_guess(4:6)

       ! Rotation matrix
       R_mat = rot_matrices(:,:, bodies(i)%ale_id)
       R_T = transpose(R_mat)

       ! Global vector from Pivot to Center of Mass
       c = matmul(R_mat, bodies(i)%local_offset_com)
       ! Global vector from Pivot to Center of Buoyancy
       r_cb = matmul(R_mat, bodies(i)%local_offset_cob)

       ! Rotate Inertia Tensor and Center of Mass (com) Offset to Global Frame
       I_P = matmul(R_mat, matmul(I_body, R_T))

!       write(msg, '(A, 3(ES23.15, 1X))') "DEBUG I_P(1,:): ", I_P(1,1), I_P(1,2), I_P(1,3)
!       call neko_log%message(trim(msg))
!       write(msg, '(A, 3(ES23.15, 1X))') "DEBUG I_P(2,:): ", I_P(2,1), I_P(2,2), I_P(2,3)
!       call neko_log%message(trim(msg))
!       write(msg, '(A, 3(ES23.15, 1X))') "DEBUG I_P(3,:): ", I_P(3,1), I_P(3,2), I_P(3,3)
!       call neko_log%message(trim(msg))


       ! skew-symmetric matrices for cross-products
       C_skew = skew_tensor(c)
       Wf_skew = skew_tensor(w_f)
       Ws_skew = skew_tensor(w_s)

       ! BDF-k History terms for acceleration
       V_hist = beta(0) * bodies(i)%body_vel_guess

!       k=0
!       write(msg, '(A, I0, A, ES23.15)') "beta(", k, ") = ", beta(k)
!       call neko_log%message(trim(msg))

       do k = 1, nadv
          V_hist = V_hist + beta(k)*bodies(i)%body_vel_lag(:, k)
 !         write(msg, '(A, I0, A, ES23.15)') "beta(", k, ") = ", beta(k)
 !         call neko_log%message(trim(msg))
       end do

       ! History acceleration (Evaluated with current guess)
       a_rel_s = V_hist(1:3) / dt !< Linear acceleration
       alpha_rel_s = V_hist(4:6) / dt !< Angular acceleration

!       write(msg, '(A, 3(ES23.15, 1X))') "DEBUG a_rel_s: ", a_rel_s(1), a_rel_s(2), a_rel_s(3)
!       call neko_log%message(trim(msg))
!       write(msg, '(A, 3(ES23.15, 1X))') "DEBUG alpha_rel_s: ", alpha_rel_s(1), alpha_rel_s(2), alpha_rel_s(3)
!       call neko_log%message(trim(msg))

       ! Initialize local matrices for this body
       M_local = 0.0_rp
       B_local = 0.0_rp

       ! ----------------------------------------------------------------------
       ! External forces
       ! ----------------------------------------------------------------------
       ! Gravity
       B_local(1:3) = B_local(1:3) + (m * gravity_vec)
       B_local(4:6) = B_local(4:6) + cross(c, m * gravity_vec)

       ! Bouyancy
       B_local(1:3) = B_local(1:3) - (m_disp * gravity_vec)
       B_local(4:6) = B_local(4:6) + cross(r_cb, -m_disp * gravity_vec)

       ! Structural Springs
       B_local(1:3) = B_local(1:3) - &
            bodies(i)%K_lin * (bodies(i)%disp_rel(1:3) - bodies(i)%pos_eq(1:3))
       B_local(4:6) = B_local(4:6) - &
            bodies(i)%K_ang * (bodies(i)%disp_rel(4:6) - bodies(i)%pos_eq(4:6))

       ! Structural Dampers
       B_local(1:3) = B_local(1:3) - bodies(i)%C_lin * v_s
       B_local(4:6) = B_local(4:6) - bodies(i)%C_ang * w_s
       do j = 1, 3
          M_local(j, j) = M_local(j, j) + bodies(i)%C_lin(j)
          M_local(j+3, j+3) = M_local(j+3, j+3) + bodies(i)%C_ang(j)
       end do

       ! Known external force applied at "pivot location"
       B_local = B_local + bodies(i)%F_prescribed_pivot


       ! ----------------------------------------------------------------------
       ! Rigid Body Inertial Forces
       ! ----------------------------------------------------------------------

       ! Term 1: Frame Linear Acceleration
       B_local(1:3) = B_local(1:3) - m * a_f

       ! Term 2: Local Linear Acceleration
       B_local(1:3) = B_local(1:3) - m * a_rel_s
       M_local(1:3, 1:3) = M_local(1:3, 1:3) + m * gamma * I3

       ! Term 3: Coriolis Acceleration
       B_local(1:3) = B_local(1:3) - 2.0_rp * m * cross(w_f, v_s)
       M_local(1:3, 1:3) = M_local(1:3, 1:3) + 2.0_rp * m * Wf_skew

       ! Term 4: Tangential due to Displacement
       B_local(1:3) = B_local(1:3) - m * cross(alpha_f, r_rel)

       ! Term 5: Centripetal due to Displacement
       B_local(1:3) = B_local(1:3) - m * cross(w_f, cross(w_f, r_rel))

       ! Term 6: Tangential Offset (Frame)
       B_local(1:3) = B_local(1:3) - m * cross(alpha_f, c)

       ! Term 7: Tangential Offset (Relative)
       B_local(1:3) = B_local(1:3) - m * cross(alpha_rel_s, c)
       M_local(1:3, 4:6) = M_local(1:3, 4:6) - m * gamma * C_skew

       ! Term 8: Gyroscopic Offset
       B_local(1:3) = B_local(1:3) - m * cross(cross(w_f, w_s), c)
       M_local(1:3, 4:6) = M_local(1:3, 4:6) - m * matmul(C_skew, Wf_skew)

       ! Term 9: Centripetal Offset (Frame)
       B_local(1:3) = B_local(1:3) - m * cross(w_f, cross(w_f, c))

       ! Term 10: Mixed Centripetal 1
       B_local(1:3) = B_local(1:3) - m * cross(w_f, cross(w_s, c))
       M_local(1:3, 4:6) = M_local(1:3, 4:6) - m * matmul(Wf_skew, C_skew)

       ! Term 11: Mixed Centripetal 2
       B_local(1:3) = B_local(1:3) - m * cross(w_s, cross(w_f, c))
       M_local(1:3, 4:6) = M_local(1:3, 4:6) - m * skew_tensor(cross(w_f, c))

       ! Term 12: Centripetal Offset
       ! (Relative - Non-Linear term dropped from LHS)
       B_local(1:3) = B_local(1:3) - m * cross(w_s, cross(w_s, c))
       M_local(1:3, 4:6) = M_local(1:3, 4:6) - m * matmul(Ws_skew, C_skew) &
            - m * skew_tensor(cross(w_s, c))

       ! ----------------------------------------------------------------------
       ! Rigid Body Inertial Torques
       ! ----------------------------------------------------------------------

       ! Torque Term 1: Frame Rotational Inertia
       B_local(4:6) = B_local(4:6) - matmul(I_P, alpha_f)

       ! Torque Term 2: Relative Rotational Inertia
       B_local(4:6) = B_local(4:6) - matmul(I_P, alpha_rel_s)
       M_local(4:6, 4:6) = M_local(4:6, 4:6) + gamma * I_P

!       write(msg, '(A, 3(ES23.15, 1X))') "I_p * alpha_ref: ", matmul(I_P, alpha_rel_s)
!       call neko_log%message(trim(msg))

       ! Torque Term 3: Gyroscopic Rotational Inertia
       B_local(4:6) = B_local(4:6) - matmul(I_P, cross(w_f, w_s))
       M_local(4:6, 4:6) = M_local(4:6, 4:6) + matmul(I_P, Wf_skew)

       ! Torque Term 4: Frame Gyroscopic
       B_local(4:6) = B_local(4:6) - cross(w_f, matmul(I_P, w_f))

       ! Torque Term 5: Mixed Gyroscopic 1
       B_local(4:6) = B_local(4:6) - cross(w_f, matmul(I_P, w_s))
       M_local(4:6, 4:6) = M_local(4:6, 4:6) + matmul(Wf_skew, I_P)

       ! Torque Term 6: Mixed Gyroscopic 2
       B_local(4:6) = B_local(4:6) - cross(w_s, matmul(I_P, w_f))
       M_local(4:6, 4:6) = M_local(4:6, 4:6) - skew_tensor(matmul(I_P, w_f))

       ! Torque Term 7: Relative Gyroscopic (Non-Linear dropped from LHS)
       B_local(4:6) = B_local(4:6) - cross(w_s, matmul(I_P, w_s))
       M_local(4:6, 4:6) = M_local(4:6, 4:6) + matmul(Ws_skew, I_P) - &
            skew_tensor(matmul(I_P, w_s))

       ! Torque Term 8: Frame Acceleration Torque
       B_local(4:6) = B_local(4:6) - m * cross(c, a_f)

       ! Torque Term 9: Local Acceleration Torque
       B_local(4:6) = B_local(4:6) - m * cross(c, a_rel_s)
       M_local(4:6, 1:3) = M_local(4:6, 1:3) + m * gamma * C_skew

       ! Torque Term 10: Coriolis Torque
       B_local(4:6) = B_local(4:6) - 2.0_rp * m * cross(c, cross(w_f, v_s))
       M_local(4:6, 1:3) = M_local(4:6, 1:3) + &
            2.0_rp * m * matmul(C_skew, Wf_skew)

       ! Torque Term 11: Displacement Tangential Torque
       B_local(4:6) = B_local(4:6) - m * cross(c, cross(alpha_f, r_rel))

       ! Torque Term 12: Displacement Centripetal Torque
       B_local(4:6) = B_local(4:6) - &
            m * cross(c, cross(w_f, cross(w_f, r_rel)))

       ! Map local matrices to global ones
       do j = 1, 6
          row_g = fsi_dof_map(i, j)
          if (row_g > 0) then
             B_global(row_g) = B_global(row_g) + B_local(j)
             do k = 1, 6
                col_g = fsi_dof_map(i, k)
                if (col_g > 0) then
                   M_global(row_g, col_g) = M_global(row_g, col_g) + &
                        M_local(j, k)
                end if
             end do
          end if
       end do

    end do

  end subroutine assemble_structural_inertial_terms
  
  !> Computes and maps the Non-Linear terms for FSI
  subroutine add_fsi_non_linear_matrices(nbodies_fsi, bodies, fsi_dof_map, &
       M_global, X_sol, rot_matrices)
    integer, intent(in) :: nbodies_fsi
    type(fsi_body_t), intent(in) :: bodies(:)
    integer, intent(in) :: fsi_dof_map(:,:)
    real(kind=rp), intent(inout) :: M_global(:,:)
    real(kind=rp), intent(in) :: X_sol(:)
    real(kind=rp), intent(in) :: rot_matrices(:,:,:)

    integer :: i, j, k, row_g, col_g
    real(kind=rp) :: R_mat(3,3), R_T(3,3), I_body(3,3), I_P(3,3)
    real(kind=rp) :: c(3), d_omega(3)
    real(kind=rp) :: C_skew(3,3), dW_skew(3,3)
    real(kind=rp) :: M_local(6,6)

    do i = 1, nbodies_fsi
       ! Extract d_omega from X_sol
       d_omega = 0.0_rp
       do k = 4, 6
          row_g = fsi_dof_map(i, k)
          if (row_g .gt. 0) d_omega(k-3) = X_sol(row_g)
       end do

       ! If no active rotational DOFs or correction is exactly zero, skip
       if (all(d_omega .eq. 0.0_rp)) cycle

       R_mat = rot_matrices(:,:, bodies(i)%ale_id)
       R_T = transpose(R_mat)
       c = matmul(R_mat, bodies(i)%local_offset_com)
       I_body = bodies(i)%I_body_tensor
       I_P = matmul(R_mat, matmul(I_body, R_T))

       dW_skew = skew_tensor(d_omega)
       C_skew = skew_tensor(c)

       M_local = 0.0_rp
       ! Term 12 Force Non-linear: -m * [dW_skew] * [C_skew]
       M_local(1:3, 4:6) = -bodies(i)%mass * matmul(dW_skew, C_skew)
       
       ! Term 4 Torque Non-linear: + [dW_skew] * I_P
       M_local(4:6, 4:6) = matmul(dW_skew, I_P)

       ! Map only the modified Non-Linear blocks to M_global
       do j = 1, 6
          row_g = fsi_dof_map(i, j)
          if (row_g .gt. 0) then
             do k = 4, 6 ! Only columns 4-6 are modified by non-linear rotation terms
                col_g = fsi_dof_map(i, k)
                if (col_g .gt. 0) then
                   M_global(row_g, col_g) = M_global(row_g, col_g) + M_local(j, k)
                end if
             end do
          end if
       end do
    end do
  end subroutine add_fsi_non_linear_matrices

  !> Computes the standard 3D cross product of two vectors
  pure function cross(a, b) result(c)
    real(kind=rp), intent(in) :: a(3), b(3)
    real(kind=rp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  !> Converts a 3D vector into its skew-symmetric tensor
  pure function skew_tensor(v) result(v_skew)
    real(kind=rp), intent(in) :: v(3)
    real(kind=rp) :: v_skew(3,3)

    v_skew(1,1) = 0.0_rp
    v_skew(2,2) = 0.0_rp
    v_skew(3,3) = 0.0_rp

    v_skew(1,2) = -v(3)
    v_skew(1,3) = v(2)

    v_skew(2,1) = v(3)
    v_skew(2,3) = -v(1)

    v_skew(3,1) = -v(2)
    v_skew(3,2) = v(1)
  end function skew_tensor

  !> Computes the Parallel Axis Theorem mapping for the inertia tensor
  subroutine apply_parallel_axis_theorem(I_in, mass, r, I_out)
    real(kind=rp), intent(in) :: I_in(3,3)
    real(kind=rp), intent(in) :: mass
    real(kind=rp), intent(in) :: r(3)
    real(kind=rp), intent(out) :: I_out(3,3)
    real(kind=rp) :: r_sq, J_steina(3,3)
    integer :: i, j

    r_sq = dot_product(r, r)
    J_steina = 0.0_rp
    do i = 1, 3
       do j = 1, 3
          J_steina(i,j) = -mass * r(i) * r(j)
       end do
       J_steina(i,i) = J_steina(i,i) + mass * r_sq
    end do
    I_out = I_in + J_steina
  end subroutine apply_parallel_axis_theorem

end module fsi_dynamics
