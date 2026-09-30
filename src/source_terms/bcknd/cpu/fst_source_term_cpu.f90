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
!> Implements the cpu kernels for the `fst_source_term_t` type.
module fst_source_term_cpu
  use num_types, only : rp
  implicit none
  private

  public :: fst_source_term_compute_cpu, fst_source_term_preview_cpu, &
       fst_fringe

contains

  !> Computes the FST fringe forcing on the cpu: adds
  !! coeff * lambda(x) * (u_bf + u' - u) at the masked points.
  !! Coordinates and velocities are read through explicit-shape dummies so
  !! that the (lx, ly, lz, nelv) field arrays can be indexed linearly in a
  !! standard-conforming way. Coordinates are the CURRENT ones (ALE-safe).
  subroutine fst_source_term_compute_cpu(n, n_mask, mask, xc, yc, zc, u, v, w, &
       fu, fv, fw, u_bf, v_bf, w_bf, k_length, kx, ky, kz, ax, ay, az, &
       phase, shift, coeff, fringe_smooth, fringe_start, &
       fringe_end, fringe_rise, fringe_fall)
    integer, intent(in) :: n, n_mask, k_length
    integer, intent(in) :: mask(n_mask)
    real(kind=rp), intent(in) :: xc(n), yc(n), zc(n)
    real(kind=rp), intent(in) :: u(n), v(n), w(n)
    real(kind=rp), intent(inout) :: fu(n), fv(n), fw(n)
    real(kind=rp), intent(in) :: u_bf(n_mask), v_bf(n_mask), w_bf(n_mask)
    real(kind=rp), intent(in) :: kx(k_length), ky(k_length), kz(k_length)
    real(kind=rp), intent(in) :: ax(k_length), ay(k_length), az(k_length)
    real(kind=rp), intent(in) :: phase(k_length)
    real(kind=rp), intent(in) :: shift(3), coeff
    logical, intent(in) :: fringe_smooth(3)
    real(kind=rp), intent(in) :: fringe_start(3), fringe_end(3), &
         fringe_rise(3), fringe_fall(3)

    integer :: idx, i
    real(kind=rp) :: lam, c, rv(3)

    do idx = 1, n_mask
       i = mask(idx)

       lam = fst_fringe(xc(i), yc(i), zc(i), fringe_smooth, fringe_start, &
            fringe_end, fringe_rise, fringe_fall)
       if (lam .le. 0.0_rp) cycle

       call fst_mode_sum(xc(i), yc(i), zc(i), shift, k_length, kx, ky, kz, &
            ax, ay, az, phase, rv)

       c = coeff*lam
       fu(i) = fu(i) + c*(u_bf(idx) + rv(1) - u(i))
       fv(i) = fv(i) + c*(v_bf(idx) + rv(2) - v(i))
       fw(i) = fw(i) + c*(w_bf(idx) + rv(3) - w(i))
    end do

  end subroutine fst_source_term_compute_cpu

  !> Sum of all Fourier modes at a point:
  !! u'_j = sum_m a_j(m) * sin(k(m).(x - U_c t) + phase(m))
  !! with shift = U_c*t (frozen turbulence at the constant velocity
  !! vector U_c, oblique/cross flow supported) and a_j the
  !! amplitude-weighted direction vectors built by `pack_modes`.
  pure subroutine fst_mode_sum(x, y, z, shift, k_length, kx, ky, kz, &
       ax, ay, az, phase, rv)
    real(kind=rp), intent(in) :: x, y, z, shift(3)
    integer, intent(in) :: k_length
    real(kind=rp), intent(in) :: kx(k_length), ky(k_length), kz(k_length)
    real(kind=rp), intent(in) :: ax(k_length), ay(k_length), az(k_length)
    real(kind=rp), intent(in) :: phase(k_length)
    real(kind=rp), intent(out) :: rv(3)

    integer :: m
    real(kind=rp) :: xs, ys, zs, sn

    xs = x - shift(1)
    ys = y - shift(2)
    zs = z - shift(3)

    rv = 0.0_rp
    do m = 1, k_length
       sn = sin(kx(m)*xs + ky(m)*ys + kz(m)*zs + phase(m))
       rv(1) = rv(1) + ax(m)*sn
       rv(2) = rv(2) + ay(m)*sn
       rv(3) = rv(3) + az(m)*sn
    end do

  end subroutine fst_mode_sum

  !> Product of the per-direction fringes at a point; flat directions
  !! contribute 1. Public so the main module can evaluate it for
  !! diagnostics.
  pure function fst_fringe(x, y, z, smooth, fstart, fend, frise, ffall) &
       result(lam)
    real(kind=rp), intent(in) :: x, y, z
    logical, intent(in) :: smooth(3)
    real(kind=rp), intent(in) :: fstart(3), fend(3), frise(3), ffall(3)
    real(kind=rp) :: lam

    real(kind=rp) :: c(3)
    integer :: d

    c(1) = x
    c(2) = y
    c(3) = z

    lam = 1.0_rp
    do d = 1, 3
       if (smooth(d)) then
          lam = lam*fringe_1d(c(d), fstart(d), fend(d), frise(d), ffall(d))
       end if
    end do

  end function fst_fringe

  !> SIMSON fringe in one direction: 0 for x <= start, ramps to 1 over
  !! rise, ramps back to 0 over fall, exactly 0 for x >= end
  !! (compact support [start, end]).
  pure function fringe_1d(x, xstart, xend, rise, fall) result(f)
    real(kind=rp), intent(in) :: x, xstart, xend, rise, fall
    real(kind=rp) :: f

    f = smooth_step((x - xstart)/rise) &
         - smooth_step((x - xend)/fall + 1.0_rp)

  end function fringe_1d

  !> Smooth step: 0 for x <= 0, 1 for x >= 1,
  !! 1/(1 + exp(1/(x-1) + 1/x)) in between (infinitely differentiable).
  pure function smooth_step(x) result(y)
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: y

    if (x .le. 0.0_rp) then
       y = 0.0_rp
    else if (x .ge. 1.0_rp) then
       y = 1.0_rp
    else
       y = 1.0_rp/(1.0_rp + exp(1.0_rp/(x - 1.0_rp) + 1.0_rp/x))
    end if

  end function smooth_step

  !> Fills the preview fields on the cpu at the masked points: the fringe
  !! lambda and the raw perturbation u' (no gain, ramp or velocity
  !! difference), for the validation dump.
  subroutine fst_source_term_preview_cpu(n, n_mask, mask, xc, yc, zc, lam_f, up, vp, wp, &
       k_length, kx, ky, kz, ax, ay, az, phase, shift, &
       fringe_smooth, fringe_start, fringe_end, fringe_rise, fringe_fall)
    integer, intent(in) :: n, n_mask, k_length
    integer, intent(in) :: mask(n_mask)
    real(kind=rp), intent(in) :: xc(n), yc(n), zc(n)
    real(kind=rp), intent(inout) :: lam_f(n), up(n), vp(n), wp(n)
    real(kind=rp), intent(in) :: kx(k_length), ky(k_length), kz(k_length)
    real(kind=rp), intent(in) :: ax(k_length), ay(k_length), az(k_length)
    real(kind=rp), intent(in) :: phase(k_length)
    real(kind=rp), intent(in) :: shift(3)
    logical, intent(in) :: fringe_smooth(3)
    real(kind=rp), intent(in) :: fringe_start(3), fringe_end(3), &
         fringe_rise(3), fringe_fall(3)

    integer :: idx, i
    real(kind=rp) :: rv(3)

    do idx = 1, n_mask
       i = mask(idx)

       lam_f(i) = fst_fringe(xc(i), yc(i), zc(i), fringe_smooth, &
            fringe_start, fringe_end, fringe_rise, fringe_fall)

       call fst_mode_sum(xc(i), yc(i), zc(i), shift, k_length, kx, ky, kz, &
            ax, ay, az, phase, rv)
       up(i) = rv(1)
       vp(i) = rv(2)
       wp(i) = rv(3)
    end do

  end subroutine fst_source_term_preview_cpu

end module fst_source_term_cpu
