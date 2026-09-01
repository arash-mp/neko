! Reproduces the fsi_cylinder example exactly, but supplies every structural
! parameter of the FSI body through the user interface
!
! Two keys are NOT part of prm and so CANNOT be set from here. They must stay
! correct in the case file:
!   "initial_velocity"  applied once at init, before any hook call
!   "active_dofs"       controls the DOF map, not a structural parameter
module user
  use neko
  implicit none

  ! Structural properties of the cylinder, matching fsi_cylinder.case.
  character(len=*), parameter :: fsi_body = 'cylFSI'

  real(kind=rp), parameter :: cyl_mass = 7.853982_rp
  real(kind=rp), parameter :: cyl_mass_disp = 0.0_rp
  real(kind=rp), parameter :: cyl_k_lin(3) = [0.0_rp, 8.64734_rp, 0.0_rp]

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fsi_structural_parameters => fsi_structural_parameters
  end subroutine user_setup

  !> Supplies the structural parameters of the FSI body.
  !!
  !! Called once per body inside every structural assembly.
  subroutine fsi_structural_parameters(body_name, body_id, time, rot_mat, &
       disp_rel, body_vel, prm)
    character(len=*), intent(in) :: body_name
    integer, intent(in) :: body_id
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: rot_mat(3, 3)
    real(kind=rp), intent(in) :: disp_rel(6)
    real(kind=rp), intent(in) :: body_vel(6)
    type(fsi_body_params_t), intent(inout) :: prm

    if (trim(body_name) .ne. fsi_body) return
    
    ! Mass and displaced (buoyancy) mass.
    prm%mass = cyl_mass
    prm%mass_disp = cyl_mass_disp

    ! Geometry, body frame, referenced to the pivot. The case file has the
    ! center of mass, the center of buoyancy and the ALE pivot all at the
    ! origin, so both offsets are zero.
    prm%offset_com = 0.0_rp
    prm%offset_cob = 0.0_rp

    ! Inertia.
    prm%inertia = 0.0_rp
    prm%inertia_ref = FSI_INERTIA_ABOUT_PIVOT

    ! Springs and dampers. Only the transverse (y) spring is non-zero.
    prm%K_lin = cyl_k_lin
    prm%K_ang = 0.0_rp
    prm%C_lin = 0.0_rp
    prm%C_ang = 0.0_rp

    ! Spring equilibrium and prescribed pivot forcing.
    prm%pos_eq = 0.0_rp
    prm%F_prescribed_pivot = 0.0_rp

  end subroutine fsi_structural_parameters

end module user
