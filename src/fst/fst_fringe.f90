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
!> Fringe for free-stream turbulence: in space, a product of one smooth step
!! per direction (a direction without an entry is flat), and in time, an
!! active window with a linear ramp.
module fst_fringe
  use num_types, only : rp
  use json_module, only : json_file, json_string
  use json_utils, only : json_get, json_get_or_lookup, &
       json_get_or_lookup_or_default
  use math, only : math_stepf
  use utils, only : neko_error, neko_warning
  use comm, only : pe_rank
  implicit none
  private

  public :: fst_fringe_value

  type, public :: fst_fringe_t
     !> Directions with a smooth fringe; the others are flat.
     logical :: smooth(3) = .false.
     !> Support [start, end] and ramp lengths per direction.
     real(kind=rp) :: start(3) = 0.0_rp
     real(kind=rp) :: end(3) = 0.0_rp
     real(kind=rp) :: rise(3) = 0.0_rp
     real(kind=rp) :: fall(3) = 0.0_rp
   contains
     procedure, pass(this) :: init => fst_fringe_init
     procedure, pass(this) :: fill => fst_fringe_fill
  end type fst_fringe_t

  !> Activation in time: active between start and end, ramped over ramp.
  type, public :: fst_ramp_t
     real(kind=rp) :: start = 0.0_rp
     real(kind=rp) :: end = huge(0.0_rp)
     real(kind=rp) :: ramp = 0.0_rp
   contains
     procedure, pass(this) :: init => fst_ramp_init
     procedure, pass(this) :: value => fst_ramp_value
  end type fst_ramp_t

contains

  !> Read the fringe under `key`: "none", or an object with optional x, y, z
  !! entries, each with start, end, rise and fall.
  !! @param json The json object holding the key.
  !! @param key The name of the fringe entry.
  !! @param required Stop if the key is missing; otherwise missing is flat.
  subroutine fst_fringe_init(this, json, key, required)
    class(fst_fringe_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: key
    logical, intent(in) :: required
    character(len=1), parameter :: dir_char(3) = ['x', 'y', 'z']
    character(len=:), allocatable :: path, str
    integer :: d, var_type
    logical :: found

    this%smooth = .false.

    call json%info(key, found = found, var_type = var_type)
    if (.not. found) then
       if (required) then
          call neko_error("(FST) '" // key // "' is required: set it " // &
               "to ""none"" or give a fringe.")
       end if
       return
    end if

    if (var_type .eq. json_string) then
       call json_get(json, key, str)
       if (trim(str) .ne. "none") then
          call neko_error("(FST) '" // key // "' must be ""none"" " // &
               "or an object")
       end if
       return
    end if

    do d = 1, 3
       path = key // "." // dir_char(d)
       if (.not. json%valid_path(path)) cycle
       this%smooth(d) = .true.
       call json_get_or_lookup(json, path // ".start", this%start(d))
       call json_get_or_lookup(json, path // ".end", this%end(d))
       call json_get_or_lookup(json, path // ".rise", this%rise(d))
       call json_get_or_lookup(json, path // ".fall", this%fall(d))

       if (this%end(d) .le. this%start(d)) then
          call neko_error("(FST) " // path // ": end must be > start")
       end if
       if (this%rise(d) .le. 0.0_rp .or. this%fall(d) .le. 0.0_rp) then
          call neko_error("(FST) " // path // ": rise and fall must be > 0")
       end if
       if (this%end(d) - this%start(d) .lt. &
            this%rise(d) + this%fall(d) .and. pe_rank .eq. 0) then
          call neko_warning("(FST) " // path // ": rise + fall exceeds" // &
               " end - start." // new_line('A') // "      The ramps " // &
               "overlap and the fringe never reaches 1.")
       end if
    end do

  end subroutine fst_fringe_init

  !> Fringe value at a point; flat directions give 1.
  function fst_fringe_value(x, y, z, smooth, fstart, fend, frise, ffall) &
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
          lam = lam*(math_stepf((c(d) - fstart(d))/frise(d)) &
               - math_stepf((c(d) - fend(d))/ffall(d) + 1.0_rp))
       end if
    end do

  end function fst_fringe_value

  !> Fringe at a list of points, lam(mask(i)) for i = 1, ..., n_mask.
  !! @param n Size of the coordinate and output arrays.
  !! @param n_mask Number of points.
  !! @param mask The points (1-based).
  !! @param xc, yc, zc Coordinates.
  !! @param lam Output.
  subroutine fst_fringe_fill(this, n, n_mask, mask, xc, yc, zc, lam)
    class(fst_fringe_t), intent(in) :: this
    integer, intent(in) :: n, n_mask
    integer, intent(in) :: mask(n_mask)
    real(kind=rp), intent(in) :: xc(n), yc(n), zc(n)
    real(kind=rp), intent(inout) :: lam(n)
    integer :: idx, i

    do idx = 1, n_mask
       i = mask(idx)
       lam(i) = fst_fringe_value(xc(i), yc(i), zc(i), this%smooth, &
            this%start, this%end, this%rise, this%fall)
    end do

  end subroutine fst_fringe_fill

  !> Read start_time, end_time and ramp_time from `json`.
  subroutine fst_ramp_init(this, json)
    class(fst_ramp_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    call json_get_or_lookup_or_default(json, "start_time", this%start, &
         0.0_rp)
    call json_get_or_lookup_or_default(json, "end_time", this%end, &
         huge(0.0_rp))
    call json_get_or_lookup_or_default(json, "ramp_time", this%ramp, 0.0_rp)
    if (this%ramp .lt. 0.0_rp) then
       call neko_error("(FST) ramp_time must be >= 0")
    end if

  end subroutine fst_ramp_init

  !> 0 before start and after end, rising linearly to 1 over the ramp.
  pure function fst_ramp_value(this, t) result(r)
    class(fst_ramp_t), intent(in) :: this
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: r

    if (t .le. this%start .or. t .gt. this%end) then
       r = 0.0_rp
    else if (this%ramp .le. 0.0_rp) then
       r = 1.0_rp
    else
       r = min(1.0_rp, (t - this%start)/this%ramp)
    end if

  end function fst_ramp_value

end module fst_fringe
