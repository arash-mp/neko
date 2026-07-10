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
!> Implicit Crank-Nicolson (trapezoidal) concrete of `ale_scheme_t`.
!>
!> This is the mesh/geometry integrator that matches a Newmark constant
!> average-acceleration structure integrator (Xu & Peet, J. Comput. Phys. X
!> 10 (2021) 100084). The structure displacement update for Newmark
!> (beta, gamma) = (1/4, 1/2) is trapezoidal (their Eq. 33),
!>   d^{n+1} = d^n + (dt/2)(d'^{n+1} + d'^n),
!> and the matching "CN-ALE" mesh map update (their Eq. 36) is
!>   x^{n+1} = x^n + (dt/2)(w^{n+1} + w^n).
!> Both the mesh reposition and the rigid-body point / 6-DOF integrators below
!> apply this same trapezoidal rule, so the mesh, the pivot / ghost trackers
!> and the logged relative displacement all move consistently.
module ale_scheme_cn
  use num_types, only : rp
  use coefs, only : coef_t
  use field, only : field_t
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  use ale_scheme, only : ale_scheme_t
  use ale_routines_cpu, only : update_ale_mesh_cn_cpu
  use ale_routines_device, only : update_ale_mesh_cn_device
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  implicit none
  private

  type, extends(ale_scheme_t), public :: ale_scheme_cn_t
   contains
     procedure, pass(this) :: reposition_mesh => cn_reposition_mesh
     procedure, pass(this) :: integrate_point => cn_integrate_point
     procedure, pass(this) :: integrate_6dof => cn_integrate_6dof
     procedure, pass(this) :: commit_history => cn_commit_history
     procedure, pass(this) :: is_implicit => cn_is_implicit
  end type ale_scheme_cn_t

contains

  !> CN mesh reposition:
  !> x^{n+1} = x^n + (dt/2)(wm^{n+1} + wm^n),
  !> with x^n = mesh_*_lag(1) and wm^n = wm_*_prev.
  subroutine cn_reposition_mesh(this, c_Xh, wm_x, wm_y, wm_z, time, nadv, &
       wm_x_lag, wm_y_lag, wm_z_lag, mesh_x_lag, mesh_y_lag, mesh_z_lag, beta, &
       wm_x_prev, wm_y_prev, wm_z_prev)
    class(ale_scheme_cn_t), intent(inout) :: this
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: wm_x, wm_y, wm_z
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    type(field_series_t), intent(in), optional :: wm_x_lag, wm_y_lag, wm_z_lag
    type(field_t), intent(in), optional :: mesh_x_lag(:), mesh_y_lag(:), &
         mesh_z_lag(:)
    real(kind=rp), intent(in), optional :: beta(0:3)
    type(field_t), intent(in), optional :: wm_x_prev, wm_y_prev, wm_z_prev

    if (.not. (present(mesh_x_lag) .and. present(mesh_y_lag) .and. &
         present(mesh_z_lag))) &
         call neko_error("ale_scheme_cn: reposition_mesh requires mesh_*_lag")
    if (.not. (present(wm_x_prev) .and. present(wm_y_prev) .and. &
         present(wm_z_prev))) &
         call neko_error("ale_scheme_cn: reposition_mesh requires wm_*_prev")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call update_ale_mesh_cn_device(c_Xh, wm_x, wm_y, wm_z, &
            wm_x_prev, wm_y_prev, wm_z_prev, &
            mesh_x_lag(1), mesh_y_lag(1), mesh_z_lag(1), time%dt)
    else
       call update_ale_mesh_cn_cpu(c_Xh, wm_x, wm_y, wm_z, &
            wm_x_prev, wm_y_prev, wm_z_prev, &
            mesh_x_lag(1), mesh_y_lag(1), mesh_z_lag(1), time%dt)
    end if
  end subroutine cn_reposition_mesh

  !> CN point integration:
  !> pos = pos^n + (dt/2)(vel + vel^n),
  !> with pos^n = hist(:,1) and vel^n = vel_prev. Overwrites pos.
  subroutine cn_integrate_point(this, pos, vel, time, nadv, vel_lag, beta, &
       hist, vel_prev)
    class(ale_scheme_cn_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: pos(3)
    real(kind=rp), intent(in) :: vel(3)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    real(kind=rp), intent(inout), optional :: vel_lag(3, 3)
    real(kind=rp), intent(in), optional :: beta(0:3)
    real(kind=rp), intent(in), optional :: hist(:, :)
    real(kind=rp), intent(in), optional :: vel_prev(3)

    if (.not. (present(hist) .and. present(vel_prev))) &
         call neko_error("ale_scheme_cn: integrate_point requires hist " // &
         "and vel_prev")

    pos = hist(:, 1) + 0.5_rp * time%dt * (vel + vel_prev)
  end subroutine cn_integrate_point

  !> CN 6-DOF displacement:
  !> x6 = x6^n + (dt/2)(v6 + v6^n),
  !> with x6^n = hist6(:,1) and v6^n = v6_prev. Overwrites x6.
  !> Matches the `pos6` combine in `subiter_update_kinematics`.
  subroutine cn_integrate_6dof(this, x6, v6, time, nadv, v6_lag, ab_coeffs, &
       beta, hist6, v6_prev)
    class(ale_scheme_cn_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: x6(6)
    real(kind=rp), intent(in) :: v6(6)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    real(kind=rp), intent(in), optional :: v6_lag(:, :)
    real(kind=rp), intent(in), optional :: ab_coeffs(:)
    real(kind=rp), intent(in), optional :: beta(0:3)
    real(kind=rp), intent(in), optional :: hist6(:, :)
    real(kind=rp), intent(in), optional :: v6_prev(6)

    if (.not. (present(hist6) .and. present(v6_prev))) &
         call neko_error("ale_scheme_cn: integrate_6dof requires hist6 " // &
         "and v6_prev")

    x6 = hist6(:, 1) + 0.5_rp * time%dt * (v6 + v6_prev)
  end subroutine cn_integrate_6dof

  !> CN commit is a no-op (implicit scheme; mesh-velocity history managed by
  !> the sub-iteration driver via the frozen wm^n snapshot).
  subroutine cn_commit_history(this, wm_x_lag, wm_y_lag, wm_z_lag)
    class(ale_scheme_cn_t), intent(inout) :: this
    type(field_series_t), intent(inout) :: wm_x_lag, wm_y_lag, wm_z_lag
    ! No op!
  end subroutine cn_commit_history

  pure function cn_is_implicit(this) result(implicit_scheme)
    class(ale_scheme_cn_t), intent(in) :: this
    logical :: implicit_scheme
    implicit_scheme = .true.
  end function cn_is_implicit

end module ale_scheme_cn