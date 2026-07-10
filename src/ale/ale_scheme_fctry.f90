! Copyright (c) 2026 The Neko Authors
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

!> Factory for the ALE mesh position-integration scheme.
!> 'ab' -> ale_scheme_ab_t (Adams-Bashforth; prescribed ALE / Green's)
!> 'bdf' -> ale_scheme_bdf_t (implicit BDF-k; FSI sub-iteration)
!> 'cn' -> ale_scheme_cn_t (implicit Crank-Nicolson; Newmark FSI sub-iteration)
module ale_scheme_fctry
  use ale_scheme, only : ale_scheme_t
  use ale_scheme_ab, only : ale_scheme_ab_t
  use ale_scheme_bdf, only : ale_scheme_bdf_t
  use ale_scheme_cn, only : ale_scheme_cn_t
  use utils, only : neko_error
  implicit none
  private

  public :: ale_scheme_factory

contains

  subroutine ale_scheme_factory(scheme, mode)
    class(ale_scheme_t), allocatable, intent(inout) :: scheme
    character(len=*), intent(in) :: mode

    if (allocated(scheme)) deallocate(scheme)

    select case (trim(mode))
    case ('ab')
       allocate(ale_scheme_ab_t :: scheme)
    case ('bdf')
       allocate(ale_scheme_bdf_t :: scheme)
    case ('cn')
       allocate(ale_scheme_cn_t :: scheme)
    case default
       call neko_error("ale_scheme_factory: unknown ALE scheme '" // &
            trim(mode) // "' (expected 'ab', 'bdf' or 'cn')")
    end select
  end subroutine ale_scheme_factory

end module ale_scheme_fctry