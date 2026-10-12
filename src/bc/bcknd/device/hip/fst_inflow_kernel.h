#ifndef __BC_FST_INFLOW_KERNEL_H__
#define __BC_FST_INFLOW_KERNEL_H__
/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

#include <fst/bcknd/device/hip/fst_device.h>

/**
 * Free-stream turbulence at the points of a boundary: g = coeff * lambda * u'.
 * The mask is a bc_t mask, 1-based with the count in element 0.
 */
template< typename T >
__global__ void fst_inflow_update_kernel(const int * __restrict__ msk,
                                         const int m,
                                         const T * __restrict__ xc,
                                         const T * __restrict__ yc,
                                         const T * __restrict__ zc,
                                         T * __restrict__ gx,
                                         T * __restrict__ gy,
                                         T * __restrict__ gz,
                                         const int k_length,
                                         const T * __restrict__ kx,
                                         const T * __restrict__ ky,
                                         const T * __restrict__ kz,
                                         const T * __restrict__ ax,
                                         const T * __restrict__ ay,
                                         const T * __restrict__ az,
                                         const T * __restrict__ phase,
                                         const T sx, const T sy, const T sz,
                                         const T coeff,
                                         const fst_fringe_t<T> fringe) {
  const T zero = 0.0;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < m; i += str) {
    const int k = msk[i + 1] - 1;
    const T lam = fst_fringe<T>(xc[k], yc[k], zc[k], fringe);
    if (lam <= zero) {
      gx[i] = zero;
      gy[i] = zero;
      gz[i] = zero;
      continue;
    }
    T rx, ry, rz;
    fst_mode_sum<T>(xc[k] - sx, yc[k] - sy, zc[k] - sz, k_length,
                    kx, ky, kz, ax, ay, az, phase, rx, ry, rz);
    const T c = coeff * lam;
    gx[i] = c * rx;
    gy[i] = c * ry;
    gz[i] = c * rz;
  }
}

/** x += g at the points of a boundary. */
template< typename T >
__global__ void fst_inflow_add_kernel(const int * __restrict__ msk,
                                      const int m,
                                      T * __restrict__ x,
                                      T * __restrict__ y,
                                      T * __restrict__ z,
                                      const T * __restrict__ gx,
                                      const T * __restrict__ gy,
                                      const T * __restrict__ gz) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < m; i += str) {
    const int k = msk[i + 1] - 1;
    x[k] += gx[i];
    y[k] += gy[i];
    z[k] += gz[i];
  }
}

#endif // __BC_FST_INFLOW_KERNEL_H__
