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
!> Device backend for `fst_inflow_t`.
module fst_inflow_device
  use num_types, only : rp, c_rp
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  public :: fst_inflow_update_device, fst_inflow_add_device

#ifdef HAVE_HIP
  interface
     subroutine hip_fst_inflow_update(m, msk, xc, yc, zc, gx, gy, gz, &
          k_length, kx, ky, kz, ax, ay, az, phase, shift, coeff, &
          fringe_smooth, fringe_start, fringe_end, fringe_rise, &
          fringe_fall, strm) bind(c, name = 'hip_fst_inflow_update')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       implicit none
       integer(c_int) :: m, k_length
       type(c_ptr), value :: msk, xc, yc, zc, gx, gy, gz
       type(c_ptr), value :: kx, ky, kz, ax, ay, az, phase
       real(c_rp) :: shift(3), coeff
       integer(c_int) :: fringe_smooth(3)
       real(c_rp) :: fringe_start(3), fringe_end(3), fringe_rise(3), &
            fringe_fall(3)
       type(c_ptr), value :: strm
     end subroutine hip_fst_inflow_update
  end interface

  interface
     subroutine hip_fst_inflow_add(m, msk, x, y, z, gx, gy, gz, strm) &
          bind(c, name = 'hip_fst_inflow_add')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z, gx, gy, gz, strm
     end subroutine hip_fst_inflow_add
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_fst_inflow_update(m, msk, xc, yc, zc, gx, gy, gz, &
          k_length, kx, ky, kz, ax, ay, az, phase, shift, coeff, &
          fringe_smooth, fringe_start, fringe_end, fringe_rise, &
          fringe_fall, strm) bind(c, name = 'cuda_fst_inflow_update')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       implicit none
       integer(c_int) :: m, k_length
       type(c_ptr), value :: msk, xc, yc, zc, gx, gy, gz
       type(c_ptr), value :: kx, ky, kz, ax, ay, az, phase
       real(c_rp) :: shift(3), coeff
       integer(c_int) :: fringe_smooth(3)
       real(c_rp) :: fringe_start(3), fringe_end(3), fringe_rise(3), &
            fringe_fall(3)
       type(c_ptr), value :: strm
     end subroutine cuda_fst_inflow_update
  end interface

  interface
     subroutine cuda_fst_inflow_add(m, msk, x, y, z, gx, gy, gz, strm) &
          bind(c, name = 'cuda_fst_inflow_add')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z, gx, gy, gz, strm
     end subroutine cuda_fst_inflow_add
  end interface
#endif

contains

  !> g = coeff * lambda * u' at the m points of a bc_t mask, on device.
  subroutine fst_inflow_update_device(m, msk, xc, yc, zc, gx, gy, gz, &
       k_length, kx, ky, kz, ax, ay, az, phase, shift, coeff, &
       fringe_smooth, fringe_start, fringe_end, fringe_rise, fringe_fall, &
       strm)
    integer, intent(in) :: m, k_length
    type(c_ptr) :: msk, xc, yc, zc, gx, gy, gz
    type(c_ptr) :: kx, ky, kz, ax, ay, az, phase, strm
    real(kind=rp), intent(in) :: shift(3), coeff
    integer, intent(in) :: fringe_smooth(3)
    real(kind=rp), intent(in) :: fringe_start(3), fringe_end(3), &
         fringe_rise(3), fringe_fall(3)

    if (m .lt. 1) return
#ifdef HAVE_HIP
    call hip_fst_inflow_update(m, msk, xc, yc, zc, gx, gy, gz, &
         k_length, kx, ky, kz, ax, ay, az, phase, shift, coeff, &
         fringe_smooth, fringe_start, fringe_end, fringe_rise, &
         fringe_fall, strm)
#elif HAVE_CUDA
    call cuda_fst_inflow_update(m, msk, xc, yc, zc, gx, gy, gz, &
         k_length, kx, ky, kz, ax, ay, az, phase, shift, coeff, &
         fringe_smooth, fringe_start, fringe_end, fringe_rise, &
         fringe_fall, strm)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for the FST inflow')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine fst_inflow_update_device

  !> x += g at the m points of a bc_t mask, on device.
  subroutine fst_inflow_add_device(m, msk, x, y, z, gx, gy, gz, strm)
    integer, intent(in) :: m
    type(c_ptr) :: msk, x, y, z, gx, gy, gz, strm

    if (m .lt. 1) return
#ifdef HAVE_HIP
    call hip_fst_inflow_add(m, msk, x, y, z, gx, gy, gz, strm)
#elif HAVE_CUDA
    call cuda_fst_inflow_add(m, msk, x, y, z, gx, gy, gz, strm)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for the FST inflow')
#else
    call neko_error('No device backend configured')
#endif

  end subroutine fst_inflow_add_device

end module fst_inflow_device
