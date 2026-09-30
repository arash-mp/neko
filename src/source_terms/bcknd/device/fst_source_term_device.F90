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
!> Device backend for the FST fringe forcing (`fst_source_term_t`).
module fst_source_term_device
  use num_types, only : rp, c_rp
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  public :: device_fst_apply

#ifdef HAVE_HIP
  interface
     subroutine hip_fst_apply(n_mask, mask_d, xc_d, yc_d, zc_d, &
          u_d, v_d, w_d, fu_d, fv_d, fw_d, u_bf_d, v_bf_d, w_bf_d, &
          k_length, kx_d, ky_d, kz_d, ax_d, ay_d, az_d, phase_d, &
          shift, coeff, fringe_smooth, fringe_start, fringe_end, &
          fringe_rise, fringe_fall) &
          bind(c, name = 'hip_fst_apply')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       implicit none
       integer(c_int) :: n_mask, k_length
       type(c_ptr), value :: mask_d, xc_d, yc_d, zc_d
       type(c_ptr), value :: u_d, v_d, w_d, fu_d, fv_d, fw_d
       type(c_ptr), value :: u_bf_d, v_bf_d, w_bf_d
       type(c_ptr), value :: kx_d, ky_d, kz_d, ax_d, ay_d, az_d, phase_d
       real(c_rp) :: shift(3), coeff
       integer(c_int) :: fringe_smooth(3)
       real(c_rp) :: fringe_start(3), fringe_end(3), fringe_rise(3), &
            fringe_fall(3)
     end subroutine hip_fst_apply
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_fst_apply(n_mask, mask_d, xc_d, yc_d, zc_d, &
          u_d, v_d, w_d, fu_d, fv_d, fw_d, u_bf_d, v_bf_d, w_bf_d, &
          k_length, kx_d, ky_d, kz_d, ax_d, ay_d, az_d, phase_d, &
          shift, coeff, fringe_smooth, fringe_start, fringe_end, &
          fringe_rise, fringe_fall) &
          bind(c, name = 'cuda_fst_apply')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       implicit none
       integer(c_int) :: n_mask, k_length
       type(c_ptr), value :: mask_d, xc_d, yc_d, zc_d
       type(c_ptr), value :: u_d, v_d, w_d, fu_d, fv_d, fw_d
       type(c_ptr), value :: u_bf_d, v_bf_d, w_bf_d
       type(c_ptr), value :: kx_d, ky_d, kz_d, ax_d, ay_d, az_d, phase_d
       real(c_rp) :: shift(3), coeff
       integer(c_int) :: fringe_smooth(3)
       real(c_rp) :: fringe_start(3), fringe_end(3), fringe_rise(3), &
            fringe_fall(3)
     end subroutine cuda_fst_apply
  end interface
#endif

contains

  !> Adds the FST fringe forcing to the right-hand-side fields on device.
  !! All array arguments are device pointers. Coordinates are the current
  !! ones, so the kernel is ALE-safe.
  !! @param n_mask Number of points in the zone (local).
  !! @param shift The frozen-turbulence shift U_c * t.
  !! @param coeff gain * ramp(t).
  !! @param fringe_smooth 1 for a smooth fringe in that direction, 0 flat.
  subroutine device_fst_apply(n_mask, mask_d, xc_d, yc_d, zc_d, &
       u_d, v_d, w_d, fu_d, fv_d, fw_d, u_bf_d, v_bf_d, w_bf_d, &
       k_length, kx_d, ky_d, kz_d, ax_d, ay_d, az_d, phase_d, &
       shift, coeff, fringe_smooth, fringe_start, fringe_end, &
       fringe_rise, fringe_fall)
    integer, intent(in) :: n_mask, k_length
    type(c_ptr) :: mask_d, xc_d, yc_d, zc_d
    type(c_ptr) :: u_d, v_d, w_d, fu_d, fv_d, fw_d
    type(c_ptr) :: u_bf_d, v_bf_d, w_bf_d
    type(c_ptr) :: kx_d, ky_d, kz_d, ax_d, ay_d, az_d, phase_d
    real(kind=rp), intent(in) :: shift(3), coeff
    integer, intent(in) :: fringe_smooth(3)
    real(kind=rp), intent(in) :: fringe_start(3), fringe_end(3), &
         fringe_rise(3), fringe_fall(3)

#ifdef HAVE_HIP
    call hip_fst_apply(n_mask, mask_d, xc_d, yc_d, zc_d, &
         u_d, v_d, w_d, fu_d, fv_d, fw_d, u_bf_d, v_bf_d, w_bf_d, &
         k_length, kx_d, ky_d, kz_d, ax_d, ay_d, az_d, phase_d, &
         shift, coeff, fringe_smooth, fringe_start, fringe_end, &
         fringe_rise, fringe_fall)
#elif HAVE_CUDA
    call cuda_fst_apply(n_mask, mask_d, xc_d, yc_d, zc_d, &
         u_d, v_d, w_d, fu_d, fv_d, fw_d, u_bf_d, v_bf_d, w_bf_d, &
         k_length, kx_d, ky_d, kz_d, ax_d, ay_d, az_d, phase_d, &
         shift, coeff, fringe_smooth, fringe_start, fringe_end, &
         fringe_rise, fringe_fall)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for the FST source term')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_fst_apply

end module fst_source_term_device