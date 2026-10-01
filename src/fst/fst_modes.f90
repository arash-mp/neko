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
!> The mode set used by the free-stream turbulence kernels: generated (and
!! saved) or read from files, packed as a_j = u_hat_j * amplitude and copied
!! to the device.
module fst_modes
  use num_types, only : rp
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, &
       json_get_or_lookup, json_get_or_lookup_or_default
  use fst_spectrum, only : fst_spectrum_t
  use coefs, only : coef_t
  use math, only : glmax, glmin
  use utils, only : neko_error, neko_warning, mkdir
  use logger, only : LOG_SIZE
  use comm, only : pe_rank
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, device_unmap, device_memcpy, &
       HOST_TO_DEVICE
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  implicit none
  private

  public :: fst_mode_sum, fst_free_mapped

  type, public :: fst_modes_t
     !> Convection velocity of the frozen turbulence and its magnitude.
     real(kind=rp) :: conv_vel(3) = 0.0_rp
     real(kind=rp) :: u_ref = 0.0_rp
     !> Largest wavenumber, for the resolution check.
     real(kind=rp) :: k_max = 0.0_rp
     !> Energy of each velocity component.
     real(kind=rp) :: energy(3) = 0.0_rp
     !> Wavenumbers, a_j = u_hat_j * amplitude and phases of the modes.
     integer :: k_length = 0
     real(kind=rp), allocatable :: kx(:), ky(:), kz(:)
     real(kind=rp), allocatable :: ax(:), ay(:), az(:)
     real(kind=rp), allocatable :: phase(:)
     !> Device copies.
     type(c_ptr) :: kx_d = C_NULL_PTR
     type(c_ptr) :: ky_d = C_NULL_PTR
     type(c_ptr) :: kz_d = C_NULL_PTR
     type(c_ptr) :: ax_d = C_NULL_PTR
     type(c_ptr) :: ay_d = C_NULL_PTR
     type(c_ptr) :: az_d = C_NULL_PTR
     type(c_ptr) :: phase_d = C_NULL_PTR
   contains
     procedure, pass(this) :: init => fst_modes_init
     procedure, pass(this) :: free => fst_modes_free
  end type fst_modes_t

contains

  !> Build the mode set.
  !! @param json The object with `convection_velocity` and `spectrum`.
  !! @param coef SEM coefficients, for the domain lengths.
  subroutine fst_modes_init(this, json, coef)
    class(fst_modes_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(json_file) :: spec_json
    type(fst_spectrum_t) :: spectrum
    real(kind=rp), allocatable :: conv_vel(:)
    character(len=*), parameter :: gen_keys(13) = [character(len=21) :: &
         "turbulence_intensity", "integral_length_scale", "n_shells", &
         "modes_per_shell", "k_min", "k_max", "periodic", "seed", &
         "output_path", "config_file", "modes_file", "phases_file", &
         "sphere_file"]
    character(len=:), allocatable :: path, config_file, modes_file, &
         phases_file, sphere_file
    logical, allocatable :: periodic_json(:)
    logical :: periodic(3)
    real(kind=rp) :: lx, ly, lz, ti, il, k_min, k_max
    integer :: n, n_shells, modes_per_shell, seed, i
    character(len=LOG_SIZE) :: log_buf

    call this%free()

    call json_get_or_lookup(json, "convection_velocity", conv_vel)
    if (size(conv_vel) .ne. 3) then
       call neko_error("(FST) convection_velocity must have 3 elements")
    end if
    call json_get(json, "spectrum", spec_json)
    this%conv_vel = conv_vel
    this%u_ref = norm2(conv_vel)
    if (this%u_ref .le. 0.0_rp) then
       call neko_error("(FST) convection_velocity must not be zero")
    end if

    n = coef%dof%size()
    lx = glmax(coef%dof%x%x, n) - glmin(coef%dof%x%x, n)
    ly = glmax(coef%dof%y%x, n) - glmin(coef%dof%y%x, n)
    lz = glmax(coef%dof%z%x, n) - glmin(coef%dof%z%x, n)

    if (spec_json%valid_path("read_from_file")) then
       do i = 1, size(gen_keys)
          if (spec_json%valid_path(trim(gen_keys(i)))) then
             call neko_error("(FST) spectrum: give either read_from_file " // &
                  "or the spectrum parameters, not both")
          end if
       end do
       call json_get(spec_json, "read_from_file.config_file", config_file)
       call json_get(spec_json, "read_from_file.modes_file", modes_file)
       call json_get(spec_json, "read_from_file.phases_file", phases_file)
       call spectrum%load(config_file, modes_file, phases_file, lx, ly, lz)

       if (abs(spectrum%u_ref - this%u_ref) .gt. &
            1.0e3_rp*epsilon(1.0_rp)*this%u_ref .and. pe_rank .eq. 0) then
          write(log_buf, '(A,E13.5,A,E13.5)') "(FST) |U_c| = ", &
               this%u_ref, ", spectrum made with U = ", spectrum%u_ref
          call neko_warning(trim(log_buf) // new_line('A') // &
               "      The turbulence intensity refers to the second.")
       end if
    else
       call json_get_or_lookup(spec_json, "turbulence_intensity", ti)
       call json_get_or_lookup(spec_json, "integral_length_scale", il)
       call json_get_or_lookup(spec_json, "n_shells", n_shells)
       call json_get_or_lookup(spec_json, "modes_per_shell", modes_per_shell)
       call json_get_or_lookup(spec_json, "k_min", k_min)
       call json_get_or_lookup(spec_json, "k_max", k_max)
       periodic = .false.
       if (spec_json%valid_path("periodic")) then
          call json_get(spec_json, "periodic", periodic_json)
          if (size(periodic_json) .ne. 3) then
             call neko_error("(FST) periodic must have 3 elements")
          end if
          periodic = periodic_json
       end if
       call json_get_or_lookup_or_default(spec_json, "seed", seed, -143)
       call json_get_or_default(spec_json, "output_path", path, "./fst_files")
       call json_get_or_default(spec_json, "config_file", config_file, &
            "fst.config")
       call json_get_or_default(spec_json, "modes_file", modes_file, &
            "fst_spectrum.csv")
       call json_get_or_default(spec_json, "phases_file", phases_file, &
            "fst_phases.csv")
       call json_get_or_default(spec_json, "sphere_file", sphere_file, &
            "sphere.dat")

       if (pe_rank .eq. 0) call mkdir(trim(path))
       sphere_file = trim(path) // "/" // sphere_file

       call spectrum%init(n_shells, modes_per_shell, k_min, k_max, ti, il, &
            this%u_ref, periodic, seed, sphere_file)
       call spectrum%generate(lx, ly, lz)
       call spectrum%save(trim(path) // "/" // config_file, &
            trim(path) // "/" // modes_file, trim(path) // "/" // phases_file)
    end if

    this%k_max = spectrum%k_end
    call pack_modes(this, spectrum)
    call spectrum%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call map_and_copy(this%kx, this%kx_d)
       call map_and_copy(this%ky, this%ky_d)
       call map_and_copy(this%kz, this%kz_d)
       call map_and_copy(this%ax, this%ax_d)
       call map_and_copy(this%ay, this%ay_d)
       call map_and_copy(this%az, this%az_d)
       call map_and_copy(this%phase, this%phase_d)
    end if

  end subroutine fst_modes_init

  !> Copy the modes into flat arrays and fold the shell amplitude into a_j.
  subroutine pack_modes(this, spectrum)
    class(fst_modes_t), intent(inout) :: this
    type(fst_spectrum_t), intent(in) :: spectrum
    integer :: m
    real(kind=rp) :: amp

    this%k_length = spectrum%k_length
    allocate(this%kx(this%k_length), this%ky(this%k_length), &
         this%kz(this%k_length))
    allocate(this%ax(this%k_length), this%ay(this%k_length), &
         this%az(this%k_length))
    allocate(this%phase(this%k_length))

    do m = 1, this%k_length
       amp = spectrum%shell_amp(spectrum%shell(m))
       this%kx(m) = spectrum%k_num(m, 1)
       this%ky(m) = spectrum%k_num(m, 2)
       this%kz(m) = spectrum%k_num(m, 3)
       this%ax(m) = spectrum%u_hat(m, 1)*amp
       this%ay(m) = spectrum%u_hat(m, 2)*amp
       this%az(m) = spectrum%u_hat(m, 3)*amp
       this%phase(m) = spectrum%phase(m)
    end do

    this%energy(1) = sum(this%ax**2)/2.0_rp
    this%energy(2) = sum(this%ay**2)/2.0_rp
    this%energy(3) = sum(this%az**2)/2.0_rp

  end subroutine pack_modes

  !> Map an array to the device and copy it there.
  subroutine map_and_copy(x, x_d)
    real(kind=rp), intent(inout) :: x(:)
    type(c_ptr), intent(inout) :: x_d

    call device_map(x, x_d, size(x))
    call device_memcpy(x, x_d, size(x), HOST_TO_DEVICE, sync = .true.)

  end subroutine map_and_copy

  !> Destructor.
  subroutine fst_modes_free(this)
    class(fst_modes_t), intent(inout) :: this

    call fst_free_mapped(this%kx, this%kx_d)
    call fst_free_mapped(this%ky, this%ky_d)
    call fst_free_mapped(this%kz, this%kz_d)
    call fst_free_mapped(this%ax, this%ax_d)
    call fst_free_mapped(this%ay, this%ay_d)
    call fst_free_mapped(this%az, this%az_d)
    call fst_free_mapped(this%phase, this%phase_d)
    this%k_length = 0

  end subroutine fst_modes_free

  !> Unmap (on device) and deallocate an array created with device_map.
  subroutine fst_free_mapped(x, x_d)
    real(kind=rp), allocatable, intent(inout) :: x(:)
    type(c_ptr), intent(inout) :: x_d

    if (allocated(x)) then
       if (c_associated(x_d)) call device_unmap(x, x_d)
       deallocate(x)
    end if

  end subroutine fst_free_mapped

  !> u'_j = sum_m a_j(m) sin(k(m) . (x - shift) + phase(m)) at one point.
  pure subroutine fst_mode_sum(x, y, z, shift, k_length, kx, ky, kz, &
       ax, ay, az, phase, rv)
    real(kind=rp), intent(in) :: x, y, z, shift(3)
    integer, intent(in) :: k_length
    real(kind=rp), intent(in) :: kx(k_length), ky(k_length), kz(k_length)
    real(kind=rp), intent(in) :: ax(k_length), ay(k_length), az(k_length)
    real(kind=rp), intent(in) :: phase(k_length)
    real(kind=rp), intent(out) :: rv(3)
    integer :: m
    real(kind=rp) :: xs, ys, zs, sn, rx, ry, rz

    xs = x - shift(1)
    ys = y - shift(2)
    zs = z - shift(3)

    rx = 0.0_rp
    ry = 0.0_rp
    rz = 0.0_rp
    !$omp simd private(sn) reduction(+:rx, ry, rz)
    do m = 1, k_length
       sn = sin(kx(m)*xs + ky(m)*ys + kz(m)*zs + phase(m))
       rx = rx + ax(m)*sn
       ry = ry + ay(m)*sn
       rz = rz + az(m)*sn
    end do

    rv(1) = rx
    rv(2) = ry
    rv(3) = rz

  end subroutine fst_mode_sum

end module fst_modes
