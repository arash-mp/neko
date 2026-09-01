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
!> Structural parameters of an FSI body.
module fsi_body_params
  use num_types, only : rp
  use utils, only : neko_error
  use math, only : NEKO_EPS
  implicit none
  private

  !> Inertia tensor reference point identifiers.
  integer, public, parameter :: FSI_INERTIA_ABOUT_PIVOT = 1
  integer, public, parameter :: FSI_INERTIA_ABOUT_COM = 2

  !> Structural parameters of a single FSI body.
  type, public :: fsi_body_params_t
     !> Mass of the body.
     real(kind=rp) :: mass = 0.0_rp
     !> Mass of the displaced fluid (buoyancy).
     real(kind=rp) :: mass_disp = 0.0_rp
     !> Reference point of inertia (FSI_INERTIA_ABOUT_PIVOT/_COM).
     integer :: inertia_ref = FSI_INERTIA_ABOUT_PIVOT
     !> Mass moment of inertia tensor, body frame, 
     !! about the point declared by inertia_ref.
     real(kind=rp) :: inertia(3, 3) = 0.0_rp
     !> Pivot -> center of mass, body frame.
     real(kind=rp) :: offset_com(3) = 0.0_rp
     !> Pivot -> center of buoyancy, body frame.
     real(kind=rp) :: offset_cob(3) = 0.0_rp
     !> Linear/angular spring stiffness.
     real(kind=rp) :: K_lin(3) = 0.0_rp
     real(kind=rp) :: K_ang(3) = 0.0_rp
     !> Linear/angular damping.
     real(kind=rp) :: C_lin(3) = 0.0_rp
     real(kind=rp) :: C_ang(3) = 0.0_rp
     !> Equilibrium position of the springs.
     real(kind=rp) :: pos_eq(6) = 0.0_rp
     !> Prescribed constant external force/torque at the pivot.
     real(kind=rp) :: F_prescribed_pivot(6) = 0.0_rp
  end type fsi_body_params_t

  public :: params_inertia_about_pivot
  public :: apply_parallel_axis_theorem
  public :: validate_body_params

contains

  !> Inertia tensor about the pivot (body frame), derived from the declared
  !! reference point.
  function params_inertia_about_pivot(prm) result(I_pivot)
    type(fsi_body_params_t), intent(in) :: prm
    real(kind=rp) :: I_pivot(3, 3)

    select case (prm%inertia_ref)
    case (FSI_INERTIA_ABOUT_PIVOT)
       I_pivot = prm%inertia
    case (FSI_INERTIA_ABOUT_COM)
       call apply_parallel_axis_theorem(prm%inertia, prm%mass, &
            prm%offset_com, I_pivot)
    case default
       call neko_error("fsi_body_params: invalid inertia_ref (must be " // &
            "FSI_INERTIA_ABOUT_PIVOT or FSI_INERTIA_ABOUT_COM)")
    end select
  end function params_inertia_about_pivot

  !> Computes the Parallel Axis Theorem mapping for the inertia tensor
  subroutine apply_parallel_axis_theorem(I_in, mass, r, I_out)
    real(kind=rp), intent(in) :: I_in(3, 3)
    real(kind=rp), intent(in) :: mass
    real(kind=rp), intent(in) :: r(3)
    real(kind=rp), intent(out) :: I_out(3, 3)
    real(kind=rp) :: r_sq, J_steina(3, 3)
    integer :: i, j

    r_sq = dot_product(r, r)
    J_steina = 0.0_rp
    do i = 1, 3
       do j = 1, 3
          J_steina(i, j) = -mass * r(i) * r(j)
       end do
       J_steina(i, i) = J_steina(i, i) + mass * r_sq
    end do
    I_out = I_in + J_steina
  end subroutine apply_parallel_axis_theorem

  !> Sanity checks on body parameters. Rejects physically
  !! impossible inputs (negative mass, non-symmetric or negative-diagonal
  !! inertia, invalid inertia_ref).
  !! NOTE: zero mass and zero inertia entries are accepted.
  subroutine validate_body_params(prm, body_name)
    type(fsi_body_params_t), intent(in) :: prm
    character(len=*), intent(in) :: body_name
    real(kind=rp) :: tol, asym
    character(len=256) :: msg
    integer :: i, j

    if (prm%mass .lt. 0.0_rp) then
       call neko_error("FSI body '" // trim(body_name) // &
            "': mass must be non-negative")
    end if

    if (prm%mass_disp .lt. 0.0_rp) then
       call neko_error("FSI body '" // trim(body_name) // &
            "': mass_disp must be non-negative")
    end if

    if (prm%inertia_ref .ne. FSI_INERTIA_ABOUT_PIVOT .and. &
         prm%inertia_ref .ne. FSI_INERTIA_ABOUT_COM) then
       call neko_error("FSI body '" // trim(body_name) // &
            "': invalid inertia_ref")
    end if

    ! Symmetry tolerance, scaled by the magnitude of the tensor so it is
    ! meaningful for both small and large inertias.
    tol = NEKO_EPS * max(abs(prm%inertia(1, 1)), abs(prm%inertia(2, 2)), &
         abs(prm%inertia(3, 3)), 1.0e-30_rp)

    do i = 1, 3
       if (prm%inertia(i, i) .lt. 0.0_rp) then
          write(msg, '(A,I0,A,I0,A,ES16.9)') &
               "': inertia tensor diagonal must be non-negative: I(", i, &
               ",", i, ") = ", prm%inertia(i, i)
          call neko_error("FSI body '" // trim(body_name) // trim(msg))
       end if
       do j = i + 1, 3
          asym = abs(prm%inertia(i, j) - prm%inertia(j, i))
          if (asym .gt. tol) then
             write(msg, '(A,I0,A,I0,A,ES16.9,A,I0,A,I0,A,ES16.9,A,ES16.9,A,&
                  &ES16.9)') &
                  "': inertia tensor must be symmetric: I(", i, ",", j, &
                  ") = ", prm%inertia(i, j), ", I(", j, ",", i, ") = ", &
                  prm%inertia(j, i), ", difference = ", asym, &
                  ", tolerance = ", tol
             call neko_error("FSI body '" // trim(body_name) // trim(msg))
          end if
       end do
    end do
  end subroutine validate_body_params

end module fsi_body_params
