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
 * Device kernel for the FST fringe forcing.
 *
 * The kernels are templated on the working precision T (instantiated with
 * `real` from device_config.h, which follows the single/double precision
 * build). All literal constants are declared as typed `const T` values so
 * that no expression is silently promoted to double in a single precision
 * build; the transcendental calls take T-typed arguments and resolve to the
 * matching device overload.
 *
 * One thread per masked (zone) point. Every thread iterates over the same
 * mode index in lockstep, so all threads in a warp read identical mode
 * addresses and the values are broadcast from cache; no shared-memory
 * staging is needed. Threads outside the fringe support (lambda = 0)
 * return early, which is safe because the kernel contains no barriers.
 */

/** Per-direction fringe description, passed by value. */
template< typename T >
struct fst_fringe_t {
  int smooth[3];
  T start[3];
  T end[3];
  T rise[3];
  T fall[3];
};

/** Smooth step: 0 for x <= 0, 1 for x >= 1, C-infinity in between. */
template< typename T >
__device__ __forceinline__ T fst_smooth_step(const T x) {
  const T zero = 0.0;
  const T one = 1.0;

  if (x <= zero) return zero;
  if (x >= one) return one;
  return one / (one + exp(one / (x - one) + one / x));
}

/** Product of the per-direction fringes; flat directions contribute 1. */
template< typename T >
__device__ __forceinline__ T fst_fringe3(const T x, const T y, const T z,
                                         const fst_fringe_t<T> f) {
  const T c[3] = {x, y, z};
  const T one = 1.0;
  T lam = one;

#pragma unroll
  for (int d = 0; d < 3; d++) {
    if (f.smooth[d]) {
      lam *= fst_smooth_step<T>((c[d] - f.start[d]) / f.rise[d])
           - fst_smooth_step<T>((c[d] - f.end[d]) / f.fall[d] + one);
    }
  }
  return lam;
}

/**
 * Adds  f_i += coeff * lambda(x) * (u_bf_i + u'_i(x,t) - u_i)
 * at the masked points, with
 *   u'_j = sum_m a_j(m) * sin(k(m).(x - shift) + phase(m)).
 * Coordinates are the current ones, so the kernel is ALE-safe.
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

  /* Fortran mask holds 1-based linear indices. */
  const int i = mask[idx] - 1;

  const T x = xc[i];
  const T y = yc[i];
  const T z = zc[i];

  const T lam = fst_fringe3<T>(x, y, z, fringe);
  if (lam <= zero) return;

  /* Shifted coordinates: k.(x - U_c t) */
  const T xs = x - sx;
  const T ys = y - sy;
  const T zs = z - sz;

  T rx = zero;
  T ry = zero;
  T rz = zero;

  for (int m = 0; m < k_length; m++) {
    const T s = sin(kx[m] * xs + ky[m] * ys + kz[m] * zs + phase[m]);
    rx += ax[m] * s;
    ry += ay[m] * s;
    rz += az[m] * s;
  }

  const T c = coeff * lam;
  fu[i] += c * (u_bf[idx] + rx - u[i]);
  fv[i] += c * (v_bf[idx] + ry - v[i]);
  fw[i] += c * (w_bf[idx] + rz - w[i]);
}