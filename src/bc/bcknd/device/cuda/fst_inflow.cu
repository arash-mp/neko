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

#include <device/device_config.h>
#include <device/cuda/check.h>
#include "fst_inflow_kernel.h"

extern "C" {

  /** g = coeff * lambda * u' at the boundary points. */
  void cuda_fst_inflow_update(int *m, void *msk,
                            void *xc, void *yc, void *zc,
                            void *gx, void *gy, void *gz,
                            int *k_length,
                            void *kx, void *ky, void *kz,
                            void *ax, void *ay, void *az, void *phase,
                            real *shift, real *coeff,
                            int *fringe_smooth, real *fringe_start,
                            real *fringe_end, real *fringe_rise,
                            real *fringe_fall,
                            cudaStream_t strm) {
    if (*m == 0) return;
    fst_fringe_t<real> fringe;
    for (int d = 0; d < 3; d++) {
      fringe.smooth[d] = fringe_smooth[d];
      fringe.start[d] = fringe_start[d];
      fringe.end[d] = fringe_end[d];
      fringe.rise[d] = fringe_rise[d];
      fringe.fall[d] = fringe_fall[d];
    }

    /* Heavy per thread: a smaller block keeps register use in bounds */
    const dim3 nthrds(256, 1, 1);
    const dim3 nblcks(((*m) + 256 - 1) / 256, 1, 1);

    fst_inflow_update_kernel<real>
    <<<nblcks, nthrds, 0, strm>>>((int *) msk, *m,
       (real *) xc, (real *) yc, (real *) zc,
       (real *) gx, (real *) gy, (real *) gz,
       *k_length,
       (real *) kx, (real *) ky, (real *) kz,
       (real *) ax, (real *) ay, (real *) az, (real *) phase,
       shift[0], shift[1], shift[2], *coeff, fringe);
    CUDA_CHECK(cudaGetLastError());
  }

  /** x += g at the boundary points. */
  void cuda_fst_inflow_add(int *m, void *msk, void *x, void *y, void *z,
                         void *gx, void *gy, void *gz, cudaStream_t strm) {
    if (*m == 0) return;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1) / 1024, 1, 1);

    fst_inflow_add_kernel<real>
    <<<nblcks, nthrds, 0, strm>>>((int *) msk, *m,
       (real *) x, (real *) y, (real *) z,
       (real *) gx, (real *) gy, (real *) gz);
    CUDA_CHECK(cudaGetLastError());
  }

}
