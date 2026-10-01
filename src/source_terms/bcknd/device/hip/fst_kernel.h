#ifndef __SOURCE_TERMS_FST_KERNEL_H__
#define __SOURCE_TERMS_FST_KERNEL_H__
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
 * FST fringe forcing, one thread per zone point. All threads walk the mode
 * list in lockstep, so mode data is broadcast from cache.
 */
/**
 * f += coeff * lambda * (u_bf + u' - u) at the zone points, with
 * u'_j = sum_m a_j(m) sin(k(m) . (x - shift) + phase(m)).
 */
template< typename T >
__global__ void fst_apply_kernel(const int n_mask,
                                 const int * __restrict__ mask,
                                 const T * __restrict__ xc,
                                 const T * __restrict__ yc,
                                 const T * __restrict__ zc,
                                 const T * __restrict__ u,
                                 const T * __restrict__ v,
                                 const T * __restrict__ w,
                                 T * __restrict__ fu,
                                 T * __restrict__ fv,
                                 T * __restrict__ fw,
                                 const T * __restrict__ u_bf,
                                 const T * __restrict__ v_bf,
                                 const T * __restrict__ w_bf,
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
  if (idx >= n_mask) return;

  /* The device copy of a mask_t is already 0-based */
  const int i = mask[idx];

  const T x = xc[i];
  const T y = yc[i];
  const T z = zc[i];

  /* No barriers in this kernel, so returning early is safe */
  const T lam = fst_fringe<T>(x, y, z, fringe);
  if (lam <= zero) return;

  const T xs = x - sx;
  const T ys = y - sy;
  const T zs = z - sz;

  T rx, ry, rz;
  fst_mode_sum<T>(xs, ys, zs, k_length, kx, ky, kz, ax, ay, az, phase,
                  rx, ry, rz);

  const T c = coeff * lam;
  fu[i] += c * (u_bf[idx] + rx - u[i]);
  fv[i] += c * (v_bf[idx] + ry - v[i]);
  fw[i] += c * (w_bf[idx] + rz - w[i]);
}

#endif // __SOURCE_TERMS_FST_KERNEL_H__
