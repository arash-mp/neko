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
     !> Previous-step acceleration, used only by the Newmark structure
     !> integrator (constant average acceleration). Ignored for the BDF path.
     real(kind=rp) :: body_acc(6) = 0.0_rp

     ! The prescribed velocity of the moving frame
     real(kind=rp) :: moving_frame_presc_acc(6) = 0.0_rp
     !> Previous-step prescribed-frame acceleration, used only by the
     !> Newmark (CN-consistent) frame-acceleration recurrence. Not checkpointed
     real(kind=rp) :: moving_frame_presc_acc_prev(6) = 0.0_rp
     real(kind=rp) :: moving_frame_presc_vel(6, 0:3) = 0.0_rp

     type(force_torque_t) :: force_monitor
  end type fsi_body_t

contains

  subroutine assemble_structural_inertial_terms(nbodies_fsi, bodies, &
       fsi_dof_map, M_global, B_global, rot_matrices, &
       time, gamma, beta, nadv, gravity_vec, accel_hist)

    integer, intent(in) :: nbodies_fsi
    type(fsi_body_t), intent(inout) :: bodies(:)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: fsi_dof_map(:,:)
    real(kind=rp), intent(inout) :: M_global(:,:), B_global(:)
    real(kind=rp), intent(in) :: rot_matrices(:,:,:)
    real(kind=rp), intent(in) :: gravity_vec(3)
    real(kind=rp) :: dt, gamma, beta(0:3)
    integer, intent(in) :: nadv
    !> Per-body history acceleration a_hist for the
    !> Newmark integrator, such that a^{n+1} = gamma*v_guess + a_hist with
    !> gamma = 2/dt.
    real(kind=rp), intent(in), optional :: accel_hist(:, :)

    integer :: i, j, k, row_g, col_g
    real(kind=rp) :: m, m_disp
    real(kind=rp) :: R_mat(3,3), R_T(3,3), I_body(3,3), I_P(3,3)
    real(kind=rp) :: c(3), r_rel(3), r_cb(3)
    real(kind=rp) :: v_s(3), w_s(3)
    real(kind=rp) :: a_f(3), alpha_f(3), w_f(3)
    real(kind=rp) :: a_rel_s(3), alpha_rel_s(3), V_hist(6), a_full(6)

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

       ! Relative (body) acceleration a^{n+1} evaluated at the current velocity
       ! guess. In both schemes it has the form  a = gamma * v_guess + a_hist,
       ! where gamma = d(a)/d(v_guess) is passed in and used in the Jacobian M.
       if (present(accel_hist)) then
          ! Newmark: a_hist = -g^n
          ! gamma = 2/dt.
          a_full = gamma * bodies(i)%body_vel_guess + accel_hist(:, i)
       else
          ! BDF-k: V_hist = beta(0)*v_guess + sum_{k>=1} beta(k)*v_lag(k),
          ! a = V_hist/dt, gamma = beta(0)/dt.
          V_hist = beta(0) * bodies(i)%body_vel_guess
          do k = 1, nadv
             V_hist = V_hist + beta(k) * bodies(i)%body_vel_lag(:, k)
          end do
          a_full = V_hist / dt
       end if

       ! History acceleration (Evaluated with current guess)
       a_rel_s = a_full(1:3) !< Linear acceleration
       alpha_rel_s = a_full(4:6) !< Angular acceleration

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

       ! ----------------------------------------------------------------------
       ! REMOVED (single-frame kinematics): Force Terms 3, 4, 5.
       ! These are the "r_rel transport" terms (Coriolis 2*w_f x v_rel,
       ! tangential alpha_f x r_rel, centripetal w_f x (w_f x r_rel)).
       ! They only arise when the body sits at an offset r_rel from a SEPARATE
       ! frame origin that the frame rotates about. In this code the body
       ! rotates about its OWN pivot (see add_kinematics_to_mesh_velocity_cpu),
       ! so that offset is identically zero and these terms are non-physical.
       ! Keeping them injected a spurious time-varying stiffness (-m*w_f^2)
       ! into the relative DOFs. See derivation note (Sec. "Single-Frame
       ! Correction").
       ! ----------------------------------------------------------------------

       ! Term 6: Tangential Offset (Frame)
       B_local(1:3) = B_local(1:3) - m * cross(alpha_f, c)

       ! Term 7: Tangential Offset (Relative)
       B_local(1:3) = B_local(1:3) - m * cross(alpha_rel_s, c)
       M_local(1:3, 4:6) = M_local(1:3, 4:6) - m * gamma * C_skew

       ! ----------------------------------------------------------------------
       ! REMOVED (lab-frame omega): Force Term 8 = m * (w_f x w_rel) x c.
       ! This came from the rotating-frame identity
       !   alpha_tot = alpha_f + alpha_rel + w_f x w_rel.
       ! Here w_f and w_rel are both stored in LAB components and differentiated
       ! with plain BDF, so alpha_tot = alpha_f + alpha_rel EXACTLY, and the
       ! w_f x w_rel cross term must NOT be added again. See derivation note.
       ! ----------------------------------------------------------------------

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

       ! ----------------------------------------------------------------------
       ! REMOVED (lab-frame omega): Torque Term 3 = I_P * (w_f x w_rel).
       ! Same reason as Force Term 8: alpha_tot = alpha_f + alpha_rel exactly
       ! when both angular velocities live in lab components, so this extra
       ! I_P (w_f x w_rel) contribution is a double-count. See derivation note.
       ! ----------------------------------------------------------------------

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

       ! ----------------------------------------------------------------------
       ! REMOVED (single-frame kinematics): Torque Terms 10, 11, 12.
       ! These are the pivot-offset torques of the removed r_rel transport
       ! forces:  m*c x [ 2 w_f x v_rel ], m*c x [ alpha_f x r_rel ],
       ! m*c x [ w_f x (w_f x r_rel) ]. They vanish for the same reason as
       ! Force Terms 3, 4, 5: the body rotates about its own pivot, so there is
       ! no r_rel offset from a separate frame origin. See derivation note.
       ! ----------------------------------------------------------------------

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
