#ifndef __FST_FST_DEVICE_H__
#define __FST_FST_DEVICE_H__
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

/**
 * Device functions shared by the free-stream turbulence kernels. Constants
 * are typed as T to avoid double promotion in single precision builds.
 */

/** Fringe description per direction, passed by value. */
template< typename T >
struct fst_fringe_t {
  int smooth[3];
  T start[3];
  T end[3];
  T rise[3];
  T fall[3];
};

/** Smooth step, same function and bounds as math_stepf on the host. */
template< typename T >
__device__ __forceinline__ T fst_stepf(const T x) {
  const T zero = 0.0;
  const T one = 1.0;
  const T xdmin = 0.0001;
  const T xdmax = 0.9999;

  if (x <= xdmin) return zero;
  if (x >= xdmax) return one;
  return one / (one + exp(one / (x - one) + one / x));
}

/** Product of the smooth fringes; flat directions contribute 1. */
template< typename T >
__device__ __forceinline__ T fst_fringe(const T x, const T y, const T z,
                                        const fst_fringe_t<T> f) {
  const T c[3] = {x, y, z};
  const T one = 1.0;
  T lam = one;

#pragma unroll
  for (int d = 0; d < 3; d++) {
    if (f.smooth[d]) {
      lam *= fst_stepf<T>((c[d] - f.start[d]) / f.rise[d])
           - fst_stepf<T>((c[d] - f.end[d]) / f.fall[d] + one);
    }
  }
  return lam;
}

/** u'_j = sum_m a_j(m) sin(k(m) . (x - shift) + phase(m)) at one point. */
template< typename T >
__device__ __forceinline__ void fst_mode_sum(const T xs, const T ys,
                                             const T zs, const int k_length,
                                             const T * __restrict__ kx,
                                             const T * __restrict__ ky,
                                             const T * __restrict__ kz,
                                             const T * __restrict__ ax,
                                             const T * __restrict__ ay,
                                             const T * __restrict__ az,
                                             const T * __restrict__ phase,
                                             T &rx, T &ry, T &rz) {
  const T zero = 0.0;
  rx = zero;
  ry = zero;
  rz = zero;

  for (int m = 0; m < k_length; m++) {
    const T s = sin(kx[m] * xs + ky[m] * ys + kz[m] * zs + phase[m]);
    rx += ax[m] * s;
    ry += ay[m] * s;
    rz += az[m] * s;
  }
}

#endif // __FST_FST_DEVICE_H__
