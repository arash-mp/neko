! Copyright (c) 2026 The Neko Authors
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

!> Adams-Bashforth concrete of `ale_scheme_t`.
module ale_scheme_ab
  use num_types, only : rp
  use coefs, only : coef_t
  use field, only : field_t
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  use ale_scheme, only : ale_scheme_t
  use ale_rigid_kinematics, only : ab_integrate_point_pos
  use ale_routines_cpu, only : update_ale_mesh_ab_cpu
  use ale_routines_device, only : update_ale_mesh_ab_device
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  implicit none
  private

  type, extends(ale_scheme_t), public :: ale_scheme_ab_t
   contains
     procedure, pass(this) :: reposition_mesh => ab_reposition_mesh
     procedure, pass(this) :: integrate_point => ab_integrate_point
     procedure, pass(this) :: integrate_6dof => ab_integrate_6dof
     procedure, pass(this) :: commit_history => ab_commit_history
     procedure, pass(this) :: is_implicit => ab_is_implicit
  end type ale_scheme_ab_t

contains

  !> AB mesh reposition:  x^{n+1} = x^n + dt * sum_j ab_coeffs(j) * wm^{n+1-j}.
  subroutine ab_reposition_mesh(this, c_Xh, wm_x, wm_y, wm_z, time, nadv, &
       wm_x_lag, wm_y_lag, wm_z_lag, mesh_x_lag, mesh_y_lag, mesh_z_lag, beta, &
       wm_x_prev, wm_y_prev, wm_z_prev)
    class(ale_scheme_ab_t), intent(inout) :: this
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: wm_x, wm_y, wm_z
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    type(field_series_t), intent(in), optional :: wm_x_lag, wm_y_lag, wm_z_lag
    type(field_t), intent(in), optional :: mesh_x_lag(:), mesh_y_lag(:), &
         mesh_z_lag(:)
    real(kind=rp), intent(in), optional :: beta(0:3)
    type(field_t), intent(in), optional :: wm_x_prev, wm_y_prev, wm_z_prev

    if (.not. (present(wm_x_lag) .and. present(wm_y_lag) .and. &
         present(wm_z_lag))) &
         call neko_error("ale_scheme_ab: reposition_mesh requires wm_*_lag")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call update_ale_mesh_ab_device(c_Xh, wm_x, wm_y, wm_z, &
            wm_x_lag, wm_y_lag, wm_z_lag, time, nadv)
    else
       call update_ale_mesh_ab_cpu(c_Xh, wm_x, wm_y, wm_z, &
            wm_x_lag, wm_y_lag, wm_z_lag, time, nadv)
    end if
  end subroutine ab_reposition_mesh

  !> AB point integration
  subroutine ab_integrate_point(this, pos, vel, time, nadv, vel_lag, beta, &
       hist, vel_prev)
    class(ale_scheme_ab_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: pos(3)
    real(kind=rp), intent(in) :: vel(3)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    real(kind=rp), intent(inout), optional :: vel_lag(3, 3)
    real(kind=rp), intent(in), optional :: beta(0:3)
    real(kind=rp), intent(in), optional :: hist(:, :)
    real(kind=rp), intent(in), optional :: vel_prev(3)

    if (.not. present(vel_lag)) &
         call neko_error("ale_scheme_ab: integrate_point requires vel_lag")

    call ab_integrate_point_pos(pos, vel_lag, vel, time, nadv)
  end subroutine ab_integrate_point

  !> x6 <- x6 + dt*ab_coeffs(1)*v6 + sum_{j=2}^{nadv} dt*ab_coeffs(j)*v6_lag(:,j)
  subroutine ab_integrate_6dof(this, x6, v6, time, nadv, v6_lag, ab_coeffs, &
       beta, hist6, v6_prev)
    class(ale_scheme_ab_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: x6(6)
    real(kind=rp), intent(in) :: v6(6)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    real(kind=rp), intent(in), optional :: v6_lag(:, :)
    real(kind=rp), intent(in), optional :: ab_coeffs(:)
    real(kind=rp), intent(in), optional :: beta(0:3)
    real(kind=rp), intent(in), optional :: hist6(:, :)
    real(kind=rp), intent(in), optional :: v6_prev(6)
    integer :: j
    real(kind=rp) :: dt

    if (.not. (present(v6_lag) .and. present(ab_coeffs))) &
         call neko_error("ale_scheme_ab: integrate_6dof requires v6_lag " // &
         "and ab_coeffs")

    dt = time%dt
    x6 = x6 + dt * ab_coeffs(1) * v6
    do j = 2, nadv
       x6 = x6 + dt * ab_coeffs(j) * v6_lag(:, j)
    end do
  end subroutine ab_integrate_6dof

  !> Updates mesh-velocity lag history.
  subroutine ab_commit_history(this, wm_x_lag, wm_y_lag, wm_z_lag)
    class(ale_scheme_ab_t), intent(inout) :: this
    type(field_series_t), intent(inout) :: wm_x_lag, wm_y_lag, wm_z_lag

    call wm_x_lag%update()
    call wm_y_lag%update()
    call wm_z_lag%update()
  end subroutine ab_commit_history

  pure function ab_is_implicit(this) result(implicit_scheme)
    class(ale_scheme_ab_t), intent(in) :: this
    logical :: implicit_scheme
    implicit_scheme = .false.
  end function ab_is_implicit

end module ale_scheme_ab