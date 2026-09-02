! Reproduces the fsi_cylinder example exactly, but supplies the entire
! structural equation
!
!     m Y" + c Y' + k (Y - Y_eq) = F_L
!
! through user%fsi_structural_terms instead of through "mass", "damping" and
! "stiffness" in the case file. Those three keys are set to zero in
! fsi_cylinder_user_terms.case, so every built-in structural term vanishes and
! the equation comes entirely from this file.
!
! -----------------------------------------------------------------------------
!
! Instructions should come here (to do)
!
! -----------------------------------------------------------------------------
!
! Not structural terms, so they cannot come from here:
!   "active_dofs"        builds the DOF map at init
!   "initial_velocity"   applied once at init, before any hook call
!   "pos_equilibrium"    still read into prm; mirrored below as y_eq because
!                        the built-in spring is switched off
!
! -----------------------------------------------------------------------------
module user
  use neko
  implicit none

  !> Must match "name" of the body in the case file.
  character(len=*), parameter :: fsi_body = 'cylFSI'

  !> The structural equation, in the cross-flow (y) direction only.
  real(kind=rp), parameter :: cyl_mass = 7.853982_rp
  real(kind=rp), parameter :: cyl_damp = 0.0_rp
  real(kind=rp), parameter :: cyl_stiff = 8.64734_rp

  !> Mirrors "pos_equilibrium"(2).
  real(kind=rp), parameter :: cyl_y_eq = 0.0_rp

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fsi_structural_terms => fsi_structural_terms
  end subroutine user_setup

  !> Adds  m Y" + c Y' + k (Y - Y_eq)  to the cross-flow equation.
  !!
  !! Called once per body on every pass of the structural fixed-point loop,
  !! with body_vel and body_acc at the current iterate. All three outputs are
  !! zeroed by the caller beforehand, so only the entries written here are
  !! non-zero.
  !!
  !! Sign convention: write the term as a generalised FORCE on the RIGHT-hand
  !! side of (inertia) = (forces). Moving  m Y" + c Y' + k dY  to the right
  !! flips its sign, hence the leading minus below.
  subroutine fsi_structural_terms(body_name, body_id, time, prm, rot_mat, &
       disp_rel, body_vel, body_acc, gravity_vec, force_pivot, &
       dforce_dvel, dforce_dacc)
    character(len=*), intent(in) :: body_name
    integer, intent(in) :: body_id
    type(time_state_t), intent(in) :: time
    type(fsi_body_params_t), intent(in) :: prm
    real(kind=rp), intent(in) :: rot_mat(3, 3)
    real(kind=rp), intent(in) :: disp_rel(6)
    real(kind=rp), intent(in) :: body_vel(6)
    real(kind=rp), intent(in) :: body_acc(6)
    real(kind=rp), intent(in) :: gravity_vec(3)
    real(kind=rp), intent(inout) :: force_pivot(6)
    real(kind=rp), intent(inout) :: dforce_dvel(6, 6)
    real(kind=rp), intent(inout) :: dforce_dacc(6, 6)

    ! This case has exactly one FSI body, so any other name is a mismatch
    ! between this file and the case file.
    if (trim(body_name) .ne. fsi_body) then
       call neko_error("fsi_cylinder_user_terms: expected body '" // &
            fsi_body // "' but the solver passed '" // trim(body_name) // &
            "'. Check the 'name' key in the case file.")
    end if

    ! ---- the term itself, on the right-hand side --------------------------
    force_pivot(2) = - cyl_mass * body_acc(2) &
         - cyl_damp * body_vel(2) &
         - cyl_stiff * (disp_rel(2) - cyl_y_eq)

    ! ---- its derivatives --------------------------------------------------
    ! d(force)/d(velocity) and d(force)/d(acceleration). The solver combines
    ! them into the Jacobian itself.
    ! dforce_dvel(i, j) = dforce_pivot(i)/dbody_vel(j)
    ! dforce_dacc(i, j) = dforce_pivot(i)/dbody_acc(j)

    ! The spring depends on disp_rel.
    dforce_dvel(2, 2) = - cyl_damp
    dforce_dacc(2, 2) = - cyl_mass

  end subroutine fsi_structural_terms

end module user
