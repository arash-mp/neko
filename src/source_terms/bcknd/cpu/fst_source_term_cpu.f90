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
!> CPU kernels for `fst_source_term_t`.
module fst_source_term_cpu
  use num_types, only : rp
  use math, only : math_stepf
  implicit none
  private

  public :: fst_source_term_compute_cpu, fst_source_term_fringe_cpu

contains

  !> f += coeff * lambda * (u_bf + u' - u) at the zone points.
  !! @param n Number of local dofs.
  !! @param mask Zone points (local linear indices).
  !! @param xc, yc, zc Current coordinates.
  !! @param shift Frozen-turbulence shift, U_c * t.
  !! @param coeff gain * ramp(t).
  subroutine fst_source_term_compute_cpu(n, n_mask, mask, xc, yc, zc, u, v, &
       w, fu, fv, fw, u_bf, v_bf, w_bf, k_length, kx, ky, kz, ax, ay, az, &
       phase, shift, coeff, fringe_smooth, fringe_start, fringe_end, &
       fringe_rise, fringe_fall)
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

    !$omp parallel do private(i, lam, c, rv)
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
    !$omp end parallel do

  end subroutine fst_source_term_compute_cpu

  !> u'_j = sum_m a_j(m) sin(k(m) . (x - shift) + phase(m)) at one point.
  pure subroutine fst_mode_sum(x, y, z, shift, k_length, kx, ky, kz, &
       ax, ay, az, phase, rv)
    real(kind=rp), intent(in) :: x, y, z, shift(3)
    integer, intent(in) :: k_length
    real(kind=rp), intent(in) :: kx(k_length), ky(k_length), kz(k_length)
    real(kind=rp), intent(in) :: ax(k_length), ay(k_length), az(k_length)
    real(kind=rp), intent(in) :: phase(k_length)
    real(kind=rp), intent(out) :: rv(3)

    integer :: m
    real(kind=rp) :: xs, ys, zs, sn, rx, ry, rz

    xs = x - shift(1)
    ys = y - shift(2)
    zs = z - shift(3)

    rx = 0.0_rp
    ry = 0.0_rp
    rz = 0.0_rp
    !$omp simd private(sn) reduction(+:rx, ry, rz)
    do m = 1, k_length
       sn = sin(kx(m)*xs + ky(m)*ys + kz(m)*zs + phase(m))
       rx = rx + ax(m)*sn
       ry = ry + ay(m)*sn
       rz = rz + az(m)*sn
    end do

    rv(1) = rx
    rv(2) = ry
    rv(3) = rz

  end subroutine fst_mode_sum

  !> Product of the smooth fringes at a point; flat directions give 1.
  function fst_fringe(x, y, z, smooth, fstart, fend, frise, ffall) &
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

  !> Rises from 0 at start over `rise`, falls back to 0 at end over `fall`.
  function fringe_1d(x, xstart, xend, rise, fall) result(f)
    real(kind=rp), intent(in) :: x, xstart, xend, rise, fall
    real(kind=rp) :: f

    f = math_stepf((x - xstart)/rise) - math_stepf((x - xend)/fall + 1.0_rp)

  end function fringe_1d

  !> Fringe lambda at the zone points, for the dump.
  subroutine fst_source_term_fringe_cpu(n, n_mask, mask, xc, yc, zc, lam, &
       fringe_smooth, fringe_start, fringe_end, fringe_rise, fringe_fall)
    integer, intent(in) :: n, n_mask
    integer, intent(in) :: mask(n_mask)
    real(kind=rp), intent(in) :: xc(n), yc(n), zc(n)
    real(kind=rp), intent(inout) :: lam(n)
    logical, intent(in) :: fringe_smooth(3)
    real(kind=rp), intent(in) :: fringe_start(3), fringe_end(3), &
         fringe_rise(3), fringe_fall(3)

    integer :: idx, i

    do idx = 1, n_mask
       i = mask(idx)
       lam(i) = fst_fringe(xc(i), yc(i), zc(i), fringe_smooth, &
            fringe_start, fringe_end, fringe_rise, fringe_fall)
    end do

  end subroutine fst_source_term_fringe_cpu

end module fst_source_term_cpu
