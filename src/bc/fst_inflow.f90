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
!> Free-stream turbulence added to a velocity inflow condition. The base
!! condition sets the mean, and ramp(t) * lambda(x) * u'(x, t) is added at
!! the boundary points, once per time step, at the current coordinates.
module fst_inflow
  use num_types, only : rp, sp
  use bc, only : bc_t, BC_DIRICHLET
  use coefs, only : coef_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, &
       json_get_or_lookup, json_get_or_lookup_or_default
  use time_state, only : time_state_t
  use field, only : field_t
  use field_series, only : field_series_t
  use fst_modes, only : fst_modes_t, fst_mode_sum, fst_free_mapped
  use fst_fringe, only : fst_fringe_t, fst_fringe_value
  use fst_inflow_device, only : fst_inflow_update_device, &
       fst_inflow_add_device
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, HOST_TO_DEVICE
  use scratch_registry, only : neko_scratch_registry
  use fld_file_output, only : fld_file_output_t
  use logger, only : neko_log
  use utils, only : neko_error, NEKO_FNAME_LEN
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  public :: fst_inflow_wrap

  type, public, extends(bc_t) :: fst_inflow_t
     !> The inflow condition the turbulence is added to.
     class(bc_t), pointer :: base => null()
     !> Mode set and fringe.
     type(fst_modes_t) :: modes
     type(fst_fringe_t) :: fringe
     !> Active between start_time and end_time, ramped over ramp_time.
     real(kind=rp) :: start_time = 0.0_rp
     real(kind=rp) :: end_time = huge(0.0_rp)
     real(kind=rp) :: ramp_time = 0.0_rp
     !> ramp * lambda * u' at the boundary points for the current step.
     logical :: active = .false.
     real(kind=rp), allocatable :: gx(:), gy(:), gz(:)
     type(c_ptr) :: gx_d = C_NULL_PTR
     type(c_ptr) :: gy_d = C_NULL_PTR
     type(c_ptr) :: gz_d = C_NULL_PTR
     !> Dump lambda and u' at the boundary at the first step.
     logical :: dump = .false.
     character(len=NEKO_FNAME_LEN) :: dump_fname
   contains
     procedure, pass(this) :: init => fst_inflow_init
     procedure, pass(this) :: free => fst_inflow_free
     procedure, pass(this) :: finalize => fst_inflow_finalize
     procedure, pass(this) :: apply_scalar => fst_inflow_apply_scalar
     procedure, pass(this) :: apply_vector => fst_inflow_apply_vector
     procedure, pass(this) :: apply_scalar_dev => &
          fst_inflow_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => &
          fst_inflow_apply_vector_dev
     procedure, pass(this) :: restart_vector => fst_inflow_restart_vector
  end type fst_inflow_t

contains

  !> Turn `object` into an fst inflow whose base is the inflow condition
  !! `object` pointed to.
  !! @param object An allocated inflow condition; on return the fst inflow.
  subroutine fst_inflow_wrap(object)
    class(bc_t), pointer, intent(inout) :: object
    type(fst_inflow_t), pointer :: fst_bc

    allocate(fst_bc)
    fst_bc%base => object
    object => fst_bc

  end subroutine fst_inflow_wrap

  !> Constructor. The base condition is set by the factory before this.
  !! @param coef SEM coefficients.
  !! @param json The boundary condition block, with the `fst` object.
  subroutine fst_inflow_init(this, coef, json)
    class(fst_inflow_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    type(json_file) :: fst, spectrum
    real(kind=rp), allocatable :: conv_vel(:)
    character(len=:), allocatable :: dump_name

    if (.not. associated(this%base)) then
       call neko_error("(FST) The inflow has no base boundary condition")
    end if
    call this%init_base(coef)
    call this%base%init(coef, json)
    this%bc_type = BC_DIRICHLET

    call neko_log%section("FST INFLOW")
    call json_get(json, "fst", fst)
    call json_get_or_lookup(fst, "convection_velocity", conv_vel)
    if (size(conv_vel) .ne. 3) then
       call neko_error("(FST) convection_velocity must have 3 elements")
    end if
    call json_get(fst, "spectrum", spectrum)
    call this%modes%init(spectrum, conv_vel, coef)
    call this%fringe%init(fst, "fringe", .true.)

    call json_get_or_lookup_or_default(fst, "start_time", this%start_time, &
         0.0_rp)
    call json_get_or_lookup_or_default(fst, "end_time", this%end_time, &
         huge(0.0_rp))
    call json_get_or_lookup_or_default(fst, "ramp_time", this%ramp_time, &
         0.0_rp)
    if (this%ramp_time .lt. 0.0_rp) then
       call neko_error("(FST) ramp_time must be >= 0")
    end if
    call json_get_or_default(fst, "dump_fields", this%dump, .false.)
    call json_get_or_default(fst, "dump_file_name", dump_name, &
         "fst_inflow_fields")
    this%dump_fname = trim(dump_name)
    call neko_log%end_section()

  end subroutine fst_inflow_init

  !> Build the masks, here and in the base, which gets the same facets.
  subroutine fst_inflow_finalize(this)
    class(fst_inflow_t), target, intent(inout) :: this
    integer :: m

    call this%base%mark_facets(this%marked_facet)
    call this%base%finalize()
    call this%finalize_base()

    m = this%msk(0)
    allocate(this%gx(m), this%gy(m), this%gz(m))
    this%gx = 0.0_rp
    this%gy = 0.0_rp
    this%gz = 0.0_rp
    if (NEKO_BCKND_DEVICE .eq. 1 .and. m .gt. 0) then
       call device_map(this%gx, this%gx_d, m)
       call device_map(this%gy, this%gy_d, m)
       call device_map(this%gz, this%gz_d, m)
    end if

  end subroutine fst_inflow_finalize

  !> Destructor.
  subroutine fst_inflow_free(this)
    class(fst_inflow_t), target, intent(inout) :: this

    if (associated(this%base)) then
       call this%base%free()
       deallocate(this%base)
    end if
    call this%modes%free()
    call fst_free_mapped(this%gx, this%gx_d)
    call fst_free_mapped(this%gy, this%gy_d)
    call fst_free_mapped(this%gz, this%gz_d)
    call this%free_base()

  end subroutine fst_inflow_free

  !> Apply the base, then add the turbulence (strong application only).
  subroutine fst_inflow_apply_vector(this, x, y, z, n, time, strong)
    class(fst_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x, y, z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_
    integer :: i, k

    strong_ = .true.
    if (present(strong)) strong_ = strong

    ! The solver resets only our flag each step, so pass it on to the base
    if (strong_ .and. .not. this%updated) this%base%updated = .false.
    call this%base%apply_vector(x, y, z, n, time, strong)
    if (.not. strong_) return

    ! On every rank, also those without points: the dump is collective
    call fst_inflow_update(this, time)
    if (.not. this%active .or. this%msk(0) .eq. 0) return
    do i = 1, this%msk(0)
       k = this%msk(i)
       x(k) = x(k) + this%gx(i)
       y(k) = y(k) + this%gy(i)
       z(k) = z(k) + this%gz(i)
    end do

  end subroutine fst_inflow_apply_vector

  !> Apply the base, then add the turbulence (device version).
  subroutine fst_inflow_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(fst_inflow_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d, y_d, z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    strong_ = .true.
    if (present(strong)) strong_ = strong

    if (strong_ .and. .not. this%updated) this%base%updated = .false.
    call this%base%apply_vector_dev(x_d, y_d, z_d, time, strong, strm)
    if (.not. strong_) return

    ! On every rank, also those without points: the dump is collective
    call fst_inflow_update(this, time, strm)
    if (.not. this%active .or. this%msk(0) .eq. 0) return
    call fst_inflow_add_device(this%msk(0), this%msk_d, x_d, y_d, z_d, &
         this%gx_d, this%gy_d, this%gz_d, strm)

  end subroutine fst_inflow_apply_vector_dev

  !> Scalar application is the base's.
  subroutine fst_inflow_apply_scalar(this, x, n, time, strong)
    class(fst_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    call this%base%apply_scalar(x, n, time, strong)

  end subroutine fst_inflow_apply_scalar

  !> Scalar application is the base's (device version).
  subroutine fst_inflow_apply_scalar_dev(this, x_d, time, strong, strm)
    class(fst_inflow_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    call this%base%apply_scalar_dev(x_d, time, strong, strm)

  end subroutine fst_inflow_apply_scalar_dev

  !> Restart is the base's.
  subroutine fst_inflow_restart_vector(this, u, v, w, ulag, vlag, wlag)
    class(fst_inflow_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w
    type(field_series_t), intent(in) :: ulag, vlag, wlag

    call this%base%restart_vector(u, v, w, ulag, vlag, wlag)

  end subroutine fst_inflow_restart_vector

  !> Compute ramp * lambda * u' at the boundary points, once per step.
  subroutine fst_inflow_update(this, time, strm)
    class(fst_inflow_t), intent(inout) :: this
    type(time_state_t), intent(in), optional :: time
    type(c_ptr), intent(inout), optional :: strm
    real(kind=rp) :: t, coeff

    if (this%updated) return
    if (.not. present(time)) then
       call neko_error("(FST) The inflow needs the time state")
    end if

    t = real(time%t, kind=rp)
    if (t .le. this%start_time .or. t .gt. this%end_time) then
       coeff = 0.0_rp
    else if (this%ramp_time .gt. 0.0_rp) then
       coeff = min(1.0_rp, (t - this%start_time)/this%ramp_time)
    else
       coeff = 1.0_rp
    end if
    this%active = coeff .gt. 0.0_rp

    if (this%dump) then
       call fst_inflow_dump(this, time)
       this%dump = .false.
    end if

    if (this%active) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          if (.not. present(strm)) then
             call neko_error("(FST) The inflow needs a device stream")
          end if
          call fst_inflow_update_device(this%msk(0), this%msk_d, &
               this%coef%dof%x%x_d, this%coef%dof%y%x_d, &
               this%coef%dof%z%x_d, this%gx_d, this%gy_d, this%gz_d, &
               this%modes%k_length, this%modes%kx_d, this%modes%ky_d, &
               this%modes%kz_d, this%modes%ax_d, this%modes%ay_d, &
               this%modes%az_d, this%modes%phase_d, &
               this%modes%conv_vel*t, coeff, merge(1, 0, &
               this%fringe%smooth), this%fringe%start, this%fringe%end, &
               this%fringe%rise, this%fringe%fall, strm)
       else
          call fst_inflow_update_cpu(this%coef%dof%size(), this%msk(0), &
               this%msk, this%coef%dof%x%x, this%coef%dof%y%x, &
               this%coef%dof%z%x, this%gx, this%gy, this%gz, &
               this%modes, this%fringe, this%modes%conv_vel*t, coeff)
       end if
    end if

    this%updated = .true.

  end subroutine fst_inflow_update

  !> g = coeff * lambda * u' at the points of a bc_t mask.
  subroutine fst_inflow_update_cpu(n, m, msk, xc, yc, zc, gx, gy, gz, &
       modes, fringe, shift, coeff)
    integer, intent(in) :: n, m
    integer, intent(in) :: msk(0:m)
    real(kind=rp), intent(in) :: xc(n), yc(n), zc(n)
    real(kind=rp), intent(inout) :: gx(m), gy(m), gz(m)
    type(fst_modes_t), intent(in) :: modes
    type(fst_fringe_t), intent(in) :: fringe
    real(kind=rp), intent(in) :: shift(3), coeff
    real(kind=rp) :: lam, rv(3)
    integer :: i, k

    !$omp parallel do private(k, lam, rv)
    do i = 1, m
       k = msk(i)
       lam = fst_fringe_value(xc(k), yc(k), zc(k), fringe%smooth, &
            fringe%start, fringe%end, fringe%rise, fringe%fall)
       if (lam .le. 0.0_rp) then
          gx(i) = 0.0_rp
          gy(i) = 0.0_rp
          gz(i) = 0.0_rp
          cycle
       end if
       call fst_mode_sum(xc(k), yc(k), zc(k), shift, modes%k_length, &
            modes%kx, modes%ky, modes%kz, modes%ax, modes%ay, modes%az, &
            modes%phase, rv)
       gx(i) = coeff*lam*rv(1)
       gy(i) = coeff*lam*rv(2)
       gz(i) = coeff*lam*rv(3)
    end do
    !$omp end parallel do

  end subroutine fst_inflow_update_cpu

  !> Write lambda, u', v', w' at the boundary points as fields 1-4. Built on
  !! the host, at the first step, before the mesh is moved.
  subroutine fst_inflow_dump(this, time)
    class(fst_inflow_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(fld_file_output_t) :: fout
    type(field_t), pointer :: f_lam, f_up, f_vp, f_wp
    integer :: i1, i2, i3, i4

    call neko_log%message("[FST] Writing inflow fields 1-4 " // &
         "(fringe, u', v', w') to '" // trim(this%dump_fname) // "'")

    call neko_scratch_registry%request_field(f_lam, i1, .false.)
    call neko_scratch_registry%request_field(f_up, i2, .false.)
    call neko_scratch_registry%request_field(f_vp, i3, .false.)
    call neko_scratch_registry%request_field(f_wp, i4, .false.)
    f_lam%x = 0.0_rp
    f_up%x = 0.0_rp
    f_vp%x = 0.0_rp
    f_wp%x = 0.0_rp

    call dump_values(f_lam%size(), this%msk(0), this%msk, &
         this%coef%dof%x%x, this%coef%dof%y%x, this%coef%dof%z%x, &
         f_lam%x, f_up%x, f_vp%x, f_wp%x, this%modes, this%fringe, &
         this%modes%conv_vel*real(time%t, kind=rp))

    call f_lam%copy_from(HOST_TO_DEVICE, .false.)
    call f_up%copy_from(HOST_TO_DEVICE, .false.)
    call f_vp%copy_from(HOST_TO_DEVICE, .false.)
    call f_wp%copy_from(HOST_TO_DEVICE, .true.)

    call fout%init(sp, trim(this%dump_fname), 4)
    call fout%fields%assign_to_ptr(1, f_lam)
    call fout%fields%assign_to_ptr(2, f_up)
    call fout%fields%assign_to_ptr(3, f_vp)
    call fout%fields%assign_to_ptr(4, f_wp)
    call fout%sample(time%t)
    call fout%free()

    call neko_scratch_registry%relinquish_field(i1)
    call neko_scratch_registry%relinquish_field(i2)
    call neko_scratch_registry%relinquish_field(i3)
    call neko_scratch_registry%relinquish_field(i4)

  end subroutine fst_inflow_dump

  !> Fringe and raw u' (no ramp) at the points of a bc_t mask.
  subroutine dump_values(n, m, msk, xc, yc, zc, lam, up, vp, wp, modes, &
       fringe, shift)
    integer, intent(in) :: n, m
    integer, intent(in) :: msk(0:m)
    real(kind=rp), intent(in) :: xc(n), yc(n), zc(n)
    real(kind=rp), intent(inout) :: lam(n), up(n), vp(n), wp(n)
    type(fst_modes_t), intent(in) :: modes
    type(fst_fringe_t), intent(in) :: fringe
    real(kind=rp), intent(in) :: shift(3)
    real(kind=rp) :: rv(3)
    integer :: i, k

    do i = 1, m
       k = msk(i)
       lam(k) = fst_fringe_value(xc(k), yc(k), zc(k), fringe%smooth, &
            fringe%start, fringe%end, fringe%rise, fringe%fall)
       call fst_mode_sum(xc(k), yc(k), zc(k), shift, modes%k_length, &
            modes%kx, modes%ky, modes%kz, modes%ax, modes%ay, modes%az, &
            modes%phase, rv)
       up(k) = rv(1)
       vp(k) = rv(2)
       wp(k) = rv(3)
    end do

  end subroutine dump_values

end module fst_inflow
