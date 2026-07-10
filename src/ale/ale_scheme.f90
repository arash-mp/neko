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
! > Abstract interface for ALE mesh-position integration schemes.
module ale_scheme
  use num_types, only : rp
  use coefs, only : coef_t
  use field, only : field_t
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  implicit none
  private

  !> Mesh position-integration strategy.
  type, abstract, public :: ale_scheme_t
   contains
     !> Recompute the mesh coordinates x^{n+1}.
     procedure(ale_reposition_intf), pass(this), deferred :: reposition_mesh
     !> Integrate a rigid point in time.
     procedure(ale_point_intf), pass(this), deferred :: integrate_point
     !> Integrate the 6-DOF displacement vector in time.
     procedure(ale_6dof_intf), pass(this), deferred :: integrate_6dof
     !> Update the mesh-velocity lag history.
     procedure(ale_commit_intf), pass(this), deferred :: commit_history
     !> .true. for BDF (implicit), .false. for AB (explicit).
     procedure(ale_flag_intf), pass(this), deferred :: is_implicit
  end type ale_scheme_t

  abstract interface

     !> Finds the new mesh position.
     !> Explicit (AB): uses `wm_*_lag` + `time`.
     !> Implicit (BDF): uses `mesh_*_lag` + `beta` + `time%dt`.
     !> Implicit (CN): uses `mesh_*_lag(1)` (= x^n) + `wm_*_prev` (= wm^n)
     !>   + `time%dt`, i.e. x^{n+1} = x^n + (dt/2)(wm^{n+1} + wm^n).
     subroutine ale_reposition_intf(this, c_Xh, wm_x, wm_y, wm_z, time, nadv, &
          wm_x_lag, wm_y_lag, wm_z_lag, mesh_x_lag, mesh_y_lag, mesh_z_lag, beta, &
          wm_x_prev, wm_y_prev, wm_z_prev)
       import :: ale_scheme_t, coef_t, field_t, field_series_t, time_state_t, rp
       class(ale_scheme_t), intent(inout) :: this
       type(coef_t), intent(inout) :: c_Xh
       type(field_t), intent(in) :: wm_x, wm_y, wm_z
       type(time_state_t), intent(in) :: time
       integer, intent(in) :: nadv
       !> Mesh-velocity lags for the AB scheme.
       type(field_series_t), intent(in), optional :: wm_x_lag, wm_y_lag, wm_z_lag
       !> Mesh-coordinate lags for the BDF/CN schemes (CN uses lag 1 only).
       type(field_t), intent(in), optional :: mesh_x_lag(:), mesh_y_lag(:), &
            mesh_z_lag(:)
       !> BDF coefficients beta(0:3) for the BDF scheme.
       real(kind=rp), intent(in), optional :: beta(0:3)
       !> Previous-step mesh velocity wm^n for the CN scheme.
       type(field_t), intent(in), optional :: wm_x_prev, wm_y_prev, wm_z_prev
     end subroutine ale_reposition_intf

     !> integrate_point: advances a point's position from its velocity.
     !> Explicit (AB) path uses `vel_lag`.
     !> Implicit (BDF) path uses `beta` + `hist`.
     !>   pos_bdf = (dt*vel - sum_j beta(j)*hist(:,j)) / beta(0).
     !> Implicit (CN) path uses `hist(:,1)` (= pos^n) + `vel_prev` (= vel^n):
     !>   pos_cn = pos^n + (dt/2)(vel + vel^n).
     subroutine ale_point_intf(this, pos, vel, time, nadv, vel_lag, beta, hist, &
          vel_prev)
       import :: ale_scheme_t, time_state_t, rp
       class(ale_scheme_t), intent(inout) :: this
       real(kind=rp), intent(inout) :: pos(3)
       real(kind=rp), intent(in) :: vel(3)
       type(time_state_t), intent(in) :: time
       integer, intent(in) :: nadv
       real(kind=rp), intent(inout), optional :: vel_lag(3, 3)
       real(kind=rp), intent(in), optional :: beta(0:3)
       real(kind=rp), intent(in), optional :: hist(:, :)
       !> Previous-step point velocity vel^n for the CN scheme.
       real(kind=rp), intent(in), optional :: vel_prev(3)
     end subroutine ale_point_intf

     !> integrate_6dof: advance the 6-DOF relative displacement `x6`.
     !> Explicit (AB): x6 <- x6 + dt*ab_coeffs(1)*v6
     !> + sum_{j=2}^{nadv} dt*ab_coeffs(j)*v6_lag(:,j)
     !> Implicit (BDF): x6 = (dt*v6 - sum_{j=1}^{nadv} beta(j)*hist6(:,j))
     !> / beta(0).
     !> Implicit (CN): x6 = hist6(:,1) + (dt/2)(v6 + v6_prev), with
     !>   hist6(:,1) = x6^n and v6_prev = v6^n.
     subroutine ale_6dof_intf(this, x6, v6, time, nadv, v6_lag, ab_coeffs, &
          beta, hist6, v6_prev)
       import :: ale_scheme_t, time_state_t, rp
       class(ale_scheme_t), intent(inout) :: this
       real(kind=rp), intent(inout) :: x6(6)
       real(kind=rp), intent(in) :: v6(6)
       type(time_state_t), intent(in) :: time
       integer, intent(in) :: nadv
       real(kind=rp), intent(in), optional :: v6_lag(:, :)
       real(kind=rp), intent(in), optional :: ab_coeffs(:)
       real(kind=rp), intent(in), optional :: beta(0:3)
       real(kind=rp), intent(in), optional :: hist6(:, :)
       !> Previous-step 6-DOF velocity v6^n for the CN scheme.
       real(kind=rp), intent(in), optional :: v6_prev(6)
     end subroutine ale_6dof_intf

     !> commit_history: update the mesh-velocity lag history.
     !> Explicit: wm_*_lag%update().
     !> Implicit: no-op.
     subroutine ale_commit_intf(this, wm_x_lag, wm_y_lag, wm_z_lag)
       import :: ale_scheme_t, field_series_t
       class(ale_scheme_t), intent(inout) :: this
       type(field_series_t), intent(inout) :: wm_x_lag, wm_y_lag, wm_z_lag
     end subroutine ale_commit_intf

     !> is_implicit:  .true. for BDF, .false. for AB.
     pure function ale_flag_intf(this) result(implicit_scheme)
       import :: ale_scheme_t
       class(ale_scheme_t), intent(in) :: this
       logical :: implicit_scheme
     end function ale_flag_intf

  end interface

end module ale_scheme