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
  use math, only : glmax, glmin, pi
  use utils, only : neko_error, neko_warning, mkdir
  use logger, only : neko_log, LOG_SIZE
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_MIN, MPI_MAX, &
       MPI_SUM, MPI_INTEGER
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
     procedure, pass(this) :: check_resolution => &
          fst_modes_check_resolution
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
    real(kind=rp) :: iso(3)
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

    iso = 3.0_rp*this%energy/sum(this%energy)
    write(log_buf, '(A,3F7.3)') &
         "[FST] component energy ratio E_u:E_v:E_w (x3/sum): ", iso
    call neko_log%message(log_buf)
    if (maxval(abs(iso - 1.0_rp)) .gt. 0.2_rp .and. pe_rank .eq. 0) then
       call neko_warning("(FST) Component energies deviate more " // &
            "than 20% from isotropy." // new_line('A') // &
            "      Increase spectrum.modes_per_shell / n_shells, or " // &
            "raise k_min relative" // new_line('A') // &
            "      to 2*pi/L in periodic directions.")
    end if

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

  !> Log how well the smallest wavelength is resolved in the elements that
  !! contain the given points, and warn below 4 points per wavelength.
  !! Collective: every rank calls it, also with no points.
  !! @param coef SEM coefficients.
  !! @param mask The points (1-based).
  !! @param n_mask Number of points on this rank.
  subroutine fst_modes_check_resolution(this, coef, mask, n_mask)
    class(fst_modes_t), intent(in) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n_mask
    integer, intent(in) :: mask(n_mask)
    real(kind=rp) :: h_min, h_max, h_avg_max, h_avg, gap_sum, lambda_min, ppw
    integer :: gap_count, ierr
    character(len=LOG_SIZE) :: log_buf

    ! The smallest wavelength has to be resolved where the grid is coarsest.
    ! Spectral elements need about pi points per wavelength of the average
    ! spacing, so the check uses the coarsest line average (4 for margin).
    call gll_spacing(coef, mask, n_mask, h_min, h_max, &
         h_avg_max, gap_sum, gap_count)
    call MPI_Allreduce(MPI_IN_PLACE, h_min, 1, MPI_REAL_PRECISION, &
         MPI_MIN, NEKO_COMM, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, h_max, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, h_avg_max, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, gap_sum, 1, MPI_REAL_PRECISION, &
         MPI_SUM, NEKO_COMM, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, gap_count, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    h_avg = gap_sum/real(gap_count, kind=rp)

    lambda_min = 2.0_rp*pi/this%k_max
    ppw = lambda_min/h_avg_max

    write(log_buf, '(A,E13.5,A,E13.5)') "[FST] GLL gap: min ", &
         h_min, " mean ", h_avg
    call neko_log%message(log_buf)
    write(log_buf, '(A,E13.5)') "[FST] GLL gap: max ", h_max
    call neko_log%message(log_buf)
    write(log_buf, '(A,E13.5)') &
         "[FST] coarsest line-average spacing (criterion basis): ", &
         h_avg_max
    call neko_log%message(log_buf)
    write(log_buf, '(A,E13.5)') "[FST] smallest FST wavelength: ", &
         lambda_min
    call neko_log%message(log_buf)
    write(log_buf, '(A,F8.2,A,F8.2)') "[FST] points per wavelength: best ", &
         lambda_min/h_min, " mean ", lambda_min/h_avg
    call neko_log%message(log_buf)
    write(log_buf, '(A,F8.2)') "[FST] points per wavelength: worst gap ", &
         lambda_min/h_max
    call neko_log%message(log_buf)
    write(log_buf, '(A,F8.2)') &
         "[FST] points per wavelength (criterion): ", ppw
    call neko_log%message(log_buf)
    if (ppw .lt. 4.0_rp .and. pe_rank .eq. 0) then
       call neko_warning("(FST) Fewer than 4 points per smallest FST" // &
            new_line('A') // "      wavelength at the coarsest average " // &
            "spacing (spectral" // new_line('A') // &
            "      criterion ~pi with margin): the smallest scales may" // &
            new_line('A') // "      be under-resolved there. " // &
            "Reduce spectrum.k_max or refine.")
    end if

  end subroutine fst_modes_check_resolution

  !> GLL spacing over the elements that contain the points: smallest and largest
  !! gap, the sum and count of all gaps, and the largest line average.
  !! Local values; reduce over ranks outside.
  subroutine gll_spacing(coef, mask, n_mask, h_min, h_max, &
       h_avg_max, gap_sum, gap_count)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n_mask
    integer, intent(in) :: mask(n_mask)
    real(kind=rp), intent(out) :: h_min, h_max, h_avg_max
    real(kind=rp), intent(out) :: gap_sum
    integer, intent(out) :: gap_count

    logical, allocatable :: in_zone(:)
    integer :: lx, ly, lz, nelv, npts, idx, e, i, j, k
    real(kind=rp) :: d, line_len

    lx = coef%Xh%lx
    ly = coef%Xh%ly
    lz = coef%Xh%lz
    nelv = coef%msh%nelv
    npts = lx*ly*lz

    allocate(in_zone(nelv))
    in_zone = .false.
    do idx = 1, n_mask
       e = (mask(idx) - 1)/npts + 1
       in_zone(e) = .true.
    end do

    h_min = huge(0.0_rp)
    h_max = 0.0_rp
    h_avg_max = 0.0_rp
    gap_sum = 0.0_rp
    gap_count = 0
    do e = 1, nelv
       if (.not. in_zone(e)) cycle

       do k = 1, lz
          do j = 1, ly
             line_len = 0.0_rp
             do i = 2, lx
                d = sqrt((coef%dof%x%x(i,j,k,e) - coef%dof%x%x(i-1,j,k,e))**2 &
                     + (coef%dof%y%x(i,j,k,e) - coef%dof%y%x(i-1,j,k,e))**2 &
                     + (coef%dof%z%x(i,j,k,e) - coef%dof%z%x(i-1,j,k,e))**2)
                h_min = min(h_min, d)
                h_max = max(h_max, d)
                gap_sum = gap_sum + d
                gap_count = gap_count + 1
                line_len = line_len + d
             end do
             h_avg_max = max(h_avg_max, line_len/real(lx - 1, kind=rp))
          end do
       end do

       do k = 1, lz
          do i = 1, lx
             line_len = 0.0_rp
             do j = 2, ly
                d = sqrt((coef%dof%x%x(i,j,k,e) - coef%dof%x%x(i,j-1,k,e))**2 &
                     + (coef%dof%y%x(i,j,k,e) - coef%dof%y%x(i,j-1,k,e))**2 &
                     + (coef%dof%z%x(i,j,k,e) - coef%dof%z%x(i,j-1,k,e))**2)
                h_min = min(h_min, d)
                h_max = max(h_max, d)
                gap_sum = gap_sum + d
                gap_count = gap_count + 1
                line_len = line_len + d
             end do
             h_avg_max = max(h_avg_max, line_len/real(ly - 1, kind=rp))
          end do
       end do

       do j = 1, ly
          do i = 1, lx
             line_len = 0.0_rp
             do k = 2, lz
                d = sqrt((coef%dof%x%x(i,j,k,e) - coef%dof%x%x(i,j,k-1,e))**2 &
                     + (coef%dof%y%x(i,j,k,e) - coef%dof%y%x(i,j,k-1,e))**2 &
                     + (coef%dof%z%x(i,j,k,e) - coef%dof%z%x(i,j,k-1,e))**2)
                h_min = min(h_min, d)
                h_max = max(h_max, d)
                gap_sum = gap_sum + d
                gap_count = gap_count + 1
                line_len = line_len + d
             end do
             h_avg_max = max(h_avg_max, line_len/real(lz - 1, kind=rp))
          end do
       end do
    end do

    deallocate(in_zone)

  end subroutine gll_spacing

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
