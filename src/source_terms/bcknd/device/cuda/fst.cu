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
#include "fst_kernel.h"

extern "C" {

  /** Launches the FST fringe forcing kernel. */
  void cuda_fst_apply(int *n_mask, void *mask,
                      void *xc, void *yc, void *zc,
                      void *u, void *v, void *w,
                      void *fu, void *fv, void *fw,
                      void *u_bf, void *v_bf, void *w_bf,
                      int *k_length,
                      void *kx, void *ky, void *kz,
                      void *ax, void *ay, void *az, void *phase,
                      real *shift, real *coeff,
                      int *fringe_smooth, real *fringe_start,
                      real *fringe_end, real *fringe_rise,
                      real *fringe_fall) {

    if (*n_mask == 0) return;

    fst_fringe_t<real> fringe;
    for (int d = 0; d < 3; d++) {
      fringe.smooth[d] = fringe_smooth[d];
      fringe.start[d] = fringe_start[d];
      fringe.end[d] = fringe_end[d];
      fringe.rise[d] = fringe_rise[d];
      fringe.fall[d] = fringe_fall[d];
    }

    const int nthrds = 256;
    const int nblcks = ((*n_mask) + nthrds - 1) / nthrds;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    fst_apply_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>(*n_mask, (int *) mask,
                                      (real *) xc, (real *) yc, (real *) zc,
                                      (real *) u, (real *) v, (real *) w,
                                      (real *) fu, (real *) fv, (real *) fw,
                                      (real *) u_bf, (real *) v_bf,
                                      (real *) w_bf,
                                      *k_length,
                                      (real *) kx, (real *) ky, (real *) kz,
                                      (real *) ax, (real *) ay, (real *) az,
                                      (real *) phase,
                                      shift[0], shift[1], shift[2],
                                      *coeff, fringe);
    CUDA_CHECK(cudaGetLastError());
  }

}
