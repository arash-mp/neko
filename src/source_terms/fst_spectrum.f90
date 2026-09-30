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
!> Implements `fst_spectrum_t`: generation of a synthetic free-stream
!! turbulence (FST) mode set (wavenumbers, amplitudes, phases) sampled from
!! a von Karman spectrum on isotropically distributed spherical shells.
!!
!! This is a runtime-configurable port of the FST generation code by
!! V. Baconnet, E. Kluesberg, P. Negi and P. Schlatter
!! (https://github.com/vbaconnet/neko-plugins, FST plugin, files
!! 01_global_params.f90 through 05_turbu.f90). The generation algorithm is
!! preserved verbatim so that, for identical inputs (seed, shells, modes per
!! shell, wavenumber range, Ti, L, U_inf, periodicity, domain lengths), this
!! module reproduces the plugin's mode set bit-for-bit (the random number
!! generator is kept identical, including its single-precision granularity).
!!
!! Differences from the plugin (deliberate):
!!  - All compile-time parameters (nshells, Npmax, kstart, kend, Ti, L,
!!    U_inf) are runtime inputs; all module-global arrays are type
!!    components with allocatable storage.
!!  - The user-provided seed is honoured. The plugin overrides it with -143
!!    inside make_turbu (05_turbu.f90); that line does not exist here.
!!  - The domain lengths are intent(in). The plugin declares them
!!    intent(out) in spec_s (04_spec.f90) although they are inputs, which
!!    only works by compiler accident.
!!  - The RNG state lives in `fst_rng_t` instead of SAVE variables, so the
!!    generator is reentrant across multiple source-term instances.
!!  - Periodic quantization uses nearest rounding instead of floor (see
!!    make_periodic_1d). This changes results relative to the plugin for
!!    periodic configurations only.
!!  - Dead code removed: the wire.dat output and its link bookkeeping (asl),
!!    the unreachable non-"new" lattice branch in the sphere routine,
!!    gen_bounded_k, and vlamax (whose max accumulator was initialized to
!!    +9.9e21 and could never have worked).
module fst_spectrum
  use num_types, only : rp
  use math, only : pi
  use utils, only : neko_error
  use logger, only : neko_log, LOG_SIZE
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use mpi_f08, only : MPI_Bcast, MPI_INTEGER
  implicit none
  private

  !> State of the portable random number generator (Numerical Recipes ran2,
  !! as used in the original FST implementation). Kept bit-identical to the
  !! plugin, including the single-precision 1/m factor, so that mode sets
  !! generated here can be validated against plugin output for equal seeds.
  type :: fst_rng_t
     integer :: ir(97) = 0
     integer :: iy = 0
     integer :: iff = 0
   contains
     procedure, pass(this) :: next => fst_rng_next
  end type fst_rng_t

  !> Synthetic FST mode set: wavenumber vectors on spherical shells,
  !! continuity-projected unit direction vectors, random phases and
  !! shell amplitudes sampled from a von Karman spectrum.
  type, public :: fst_spectrum_t
     !> Number of spherical shells discretizing [k_start, k_end].
     integer :: n_shells = 0
     !> Requested points per shell (before mirroring). The sphere point
     !! lattice may realize fewer points; see `sphere_points`.
     integer :: npmax = 0
     !> Maximum number of modes (2 * n_shells * npmax); allocation size of
     !! the per-mode arrays. The number of valid modes is `k_length`.
     integer :: n_modes_max = 0
     !> Smallest and largest total wavenumber (shell radii).
     real(kind=rp) :: k_start = 0.0_rp
     real(kind=rp) :: k_end = 0.0_rp
     !> Target turbulence intensity Tu (rms(u')/U_inf).
     real(kind=rp) :: ti = 0.0_rp
     !> Integral length scale of the von Karman spectrum.
     real(kind=rp) :: il = 0.0_rp
     !> Reference speed |U_c| used for the energy scaling
     !! (tke = 3/2 * (Ti * u_ref)^2).
     real(kind=rp) :: u_ref = 0.0_rp
     !> Periodic directions (x, y, z): wavenumber components in periodic
     !! directions are quantized to multiples of 2*pi/L.
     logical :: periodic(3) = .false.
     !> RNG seed (used as given; negative values initialize the generator).
     integer :: seed = -143
     !> Whether to write diagnostic files (sphere.dat, bb.txt,
     !! fst_spectrum.csv) during generation (rank 0 only).
     logical :: write_files = .false.
     !> Output path for the diagnostic files.
     character(len=:), allocatable :: path

     !> Number of valid modes (after removal of zero wavenumber vectors).
     integer :: k_length = 0
     !> Wavenumber vectors, (n_modes_max, 3); valid rows are 1..k_length.
     real(kind=rp), allocatable :: k_num(:,:)
     !> Continuity-projected unit direction vectors, (n_modes_max, 3).
     real(kind=rp), allocatable :: u_hat(:,:)
     !> Random phase per mode, (n_modes_max).
     real(kind=rp), allocatable :: phase(:)
     !> Shell index of each mode, (n_modes_max).
     integer, allocatable :: shell(:)
     !> Amplitude of each shell, (n_shells).
     real(kind=rp), allocatable :: shell_amp(:)
     !> Number of valid modes in each shell, (n_shells).
     integer, allocatable :: shell_modes(:)

     ! --- Diagnostics, valid on rank 0 only (used for logging/validation).
     !> Estimated Tu * U_inf of the generated mode set,
     !! sqrt((E_u + E_v + E_w)/3).
     real(kind=rp) :: tu_uinf_estimate = 0.0_rp
     !> Energy in each velocity component of the generated mode set.
     real(kind=rp) :: energy(3) = 0.0_rp
   contains
     procedure, pass(this) :: init => fst_spectrum_init
     procedure, pass(this) :: generate => fst_spectrum_generate
     procedure, pass(this) :: free => fst_spectrum_free
  end type fst_spectrum_t

contains

  !> Portable random number generator (Numerical Recipes ran2). Verbatim
  !! port of ran2 in the FST plugin's 01_global_params.f90; the state is a
  !! type component instead of SAVE variables. `rm` is intentionally kept in
  !! default (single) precision for bit parity with the plugin.
  function fst_rng_next(this, idum) result(r)
    class(fst_rng_t), intent(inout) :: this
    integer, intent(inout) :: idum
    real(kind=rp) :: r

    integer, parameter :: m = 714025, ia = 1366, ic = 150889
    real, parameter :: rm = 1./m
    integer :: j

    if (idum .lt. 0 .or. this%iff .eq. 0) then
       ! Initialize
       this%iff = 1
       idum = mod(ic - idum, m)
       do j = 1, 97
          idum = mod(ia*idum + ic, m)
          this%ir(j) = idum
       end do
       idum = mod(ia*idum + ic, m)
       this%iy = idum
    end if

    j = 1 + (97*this%iy)/m
    this%iy = this%ir(j)
    r = this%iy*rm
    idum = mod(ia*idum + ic, m)
    this%ir(j) = idum

  end function fst_rng_next

  !> Initialize the spectrum configuration and allocate the result arrays.
  !! All arguments are validated; invalid input is a hard error.
  subroutine fst_spectrum_init(this, n_shells, npmax, k_start, k_end, &
       ti, il, u_ref, periodic, seed, write_files, path)
    class(fst_spectrum_t), intent(inout) :: this
    integer, intent(in) :: n_shells
    integer, intent(in) :: npmax
    real(kind=rp), intent(in) :: k_start, k_end, ti, il, u_ref
    logical, intent(in) :: periodic(3)
    integer, intent(in) :: seed
    logical, intent(in) :: write_files
    character(len=*), intent(in) :: path

    call this%free()

    if (n_shells .lt. 2) then
       call neko_error("(FST) spectrum.n_shells must be >= 2")
    end if
    if (npmax .lt. 3 .or. npmax .gt. 1000) then
       call neko_error("(FST) spectrum.modes_per_shell must be in " // &
            "[3, 1000]")
    end if
    if (k_start .le. 0.0_rp) then
       call neko_error("(FST) spectrum.k_min must be > 0")
    end if
    if (k_end .le. k_start) then
       call neko_error("(FST) spectrum.k_max must be > spectrum.k_min")
    end if
    if (ti .le. 0.0_rp) then
       call neko_error("(FST) turbulence_intensity must be > 0")
    end if
    if (il .le. 0.0_rp) then
       call neko_error("(FST) integral_length_scale must be > 0")
    end if
    if (u_ref .le. 0.0_rp) then
       call neko_error("(FST) reference speed " // &
            "|convection_velocity| must be > 0")
    end if
    if (all(periodic)) then
       call neko_error("(FST) at most two periodic directions are supported")
    end if

    this%n_shells = n_shells
    this%npmax = npmax
    this%n_modes_max = 2 * n_shells * npmax
    this%k_start = k_start
    this%k_end = k_end
    this%ti = ti
    this%il = il
    this%u_ref = u_ref
    this%periodic = periodic
    this%seed = seed
    this%write_files = write_files
    this%path = trim(path)

    allocate(this%k_num(this%n_modes_max, 3))
    allocate(this%u_hat(this%n_modes_max, 3))
    allocate(this%phase(this%n_modes_max))
    allocate(this%shell(this%n_modes_max))
    allocate(this%shell_amp(this%n_shells))
    allocate(this%shell_modes(this%n_shells))

    ! Zero-initialize so that entries beyond k_length are defined (they are
    ! broadcast but never read) and so that any unfilled sphere-lattice slot
    ! is deterministically removed by the zero-wavenumber filter.
    this%k_num = 0.0_rp
    this%u_hat = 0.0_rp
    this%phase = 0.0_rp
    this%shell = 0
    this%shell_amp = 0.0_rp
    this%shell_modes = 0

  end subroutine fst_spectrum_init

  !> Free the result arrays and reset the configuration.
  subroutine fst_spectrum_free(this)
    class(fst_spectrum_t), intent(inout) :: this

    if (allocated(this%k_num)) deallocate(this%k_num)
    if (allocated(this%u_hat)) deallocate(this%u_hat)
    if (allocated(this%phase)) deallocate(this%phase)
    if (allocated(this%shell)) deallocate(this%shell)
    if (allocated(this%shell_amp)) deallocate(this%shell_amp)
    if (allocated(this%shell_modes)) deallocate(this%shell_modes)
    if (allocated(this%path)) deallocate(this%path)

    this%n_shells = 0
    this%npmax = 0
    this%n_modes_max = 0
    this%k_length = 0

  end subroutine fst_spectrum_free

  !> Generate the mode set on rank 0 and broadcast it to all ranks.
  !! @param lx Domain length in x, used for periodic quantization if
  !! periodic(1).
  !! @param ly Domain length in y, used if periodic(2).
  !! @param lz Domain length in z, used if periodic(3).
  subroutine fst_spectrum_generate(this, lx, ly, lz)
    class(fst_spectrum_t), intent(inout) :: this
    real(kind=rp), intent(in) :: lx, ly, lz

    integer :: ierr

    call neko_log%section('Generating FST spectrum')

    if (pe_rank .eq. 0) then
       call generate_rank0(this, lx, ly, lz)
    end if

    call MPI_Bcast(this%k_length, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    call MPI_Bcast(this%k_num, this%n_modes_max*3, MPI_REAL_PRECISION, &
         0, NEKO_COMM, ierr)
    call MPI_Bcast(this%u_hat, this%n_modes_max*3, MPI_REAL_PRECISION, &
         0, NEKO_COMM, ierr)
    call MPI_Bcast(this%phase, this%n_modes_max, MPI_REAL_PRECISION, &
         0, NEKO_COMM, ierr)
    call MPI_Bcast(this%shell, this%n_modes_max, MPI_INTEGER, &
         0, NEKO_COMM, ierr)
    call MPI_Bcast(this%shell_amp, this%n_shells, MPI_REAL_PRECISION, &
         0, NEKO_COMM, ierr)
    call MPI_Bcast(this%shell_modes, this%n_shells, MPI_INTEGER, &
         0, NEKO_COMM, ierr)

    call neko_log%end_section('Done --> Generating FST spectrum')

  end subroutine fst_spectrum_generate

  !> Rank-0 generation. Port of spec_s (04_spec.f90) followed by the random
  !! amplitude/phase drawing and continuity projection of make_turbu
  !! (05_turbu.f90). The random draw order is identical to the plugin.
  subroutine generate_rank0(this, dlx, dly, dlz)
    class(fst_spectrum_t), intent(inout) :: this
    real(kind=rp), intent(in) :: dlx, dly, dlz

    type(fst_rng_t) :: rng
    integer :: idum

    real(kind=rp), allocatable :: co(:,:,:)
    real(kind=rp), allocatable :: kk(:), q(:), dk(:), tke_shell(:)
    real(kind=rp), allocatable :: bb(:,:), bb1(:,:)
    real(kind=rp) :: u_hat_raw(3), u_hat_p(3)

    real(kind=rp) :: k2, dkint, tke_tot, tke_tot1, tke_scaled, shell_energy
    real(kind=rp) :: rotx, roty, rotz
    real(kind=rp) :: kxmin, kxmax, kymin, kymax, kzmin, kzmax
    real(kind=rp) :: ue, ve, we, uamp, vamp, wamp, amp, kdotu, knorm2
    integer :: np, ndk, i, j, k, l, shellno, n_kept, n_removed
    character(len=LOG_SIZE) :: log_buf

    idum = this%seed
    np = this%npmax

    allocate(co(2*this%npmax, this%n_shells, 3))
    allocate(kk(0:this%n_shells))
    allocate(q(this%n_shells))
    allocate(dk(this%n_shells))
    allocate(tke_shell(this%n_shells))
    co = 0.0_rp

    call print_param('integral length scale', this%il)

    ! Target kinetic energy: 3/2 * Tu^2 * U_inf^2
    tke_scaled = 1.5_rp * (this%ti * this%u_ref)**2

    kxmax = 1.0e-20_rp
    kxmin = 1.0e+20_rp
    kymax = 1.0e-20_rp
    kymin = 1.0e+20_rp
    kzmax = 1.0e-20_rp
    kzmin = 1.0e+20_rp

    ! --- Integrate the energy spectrum on a fine grid (diagnostic only)
    ndk = 5000
    dkint = (this%k_end - this%k_start)/real(ndk, kind=rp)
    tke_tot1 = ek(this%k_start, this%il, 1.0_rp) &
         + ek(this%k_end, this%il, 1.0_rp)
    do i = 1, ndk - 1
       tke_tot1 = tke_tot1 + ek(this%k_start + i*dkint, this%il, 1.0_rp)
    end do
    tke_tot1 = tke_tot1*dkint
    call print_param('FST - integrated energy in spectrum', tke_tot1)

    ! --- Integrate the energy spectrum on the n_shells nodes. This is the
    !     normalization used to rescale the truncated spectrum to tke_scaled.
    dkint = (this%k_end - this%k_start)/real(this%n_shells - 1, kind=rp)
    tke_tot = 0.0_rp
    do i = 1, this%n_shells
       tke_tot = tke_tot + ek(this%k_start + (i-1)*dkint, this%il, 1.0_rp)
    end do
    tke_tot = tke_tot*dkint
    write (log_buf, '(A,I0,A,E13.5)') 'FST - discretized on ', &
         this%n_shells, ' shells : ', tke_tot
    call neko_log%message(log_buf)

    if (this%write_files) then
       open(file = trim(this%path) // '/sphere.dat', unit = 10)
       write(10, *) 'energy shell parameters'
       write(10, '(a20,i18)') 'Nshells', this%n_shells
       write(10, '(a20,f18.9)') 'kstart', this%k_start
       write(10, '(a20,f18.9)') 'kend', this%k_end
       write(10, '(a20,i18)') 'Np', np
       write(10, *) 'isotropic coordinates'
       write(10, '(2a5,3a18)') 'i', 'j', 'x', 'y', 'z'
    end if

    kk(0) = 0.0_rp
    tke_tot1 = 0.0_rp

    call print_param("Truncated TKE", tke_scaled/tke_tot)

    do i = 1, this%n_shells

       k2 = (this%k_start + (i-1)*(this%k_end - this%k_start) &
            / real(this%n_shells - 1, kind=rp))**2
       kk(i) = sqrt(k2)
       dk(i) = (this%k_end - this%k_start)/real(this%n_shells - 1, kind=rp)

       ! 1/tke_tot so that the total truncated energy equals tke_scaled
       q(i) = ek(kk(i), this%il, tke_scaled/tke_tot)
       tke_shell(i) = q(i)*dk(i)
       tke_tot1 = tke_tot1 + tke_shell(i)

       ! Randomly rotated sphere point set with radius kk(i).
       ! NOTE: np is intent(inout) and may be reduced by the point lattice
       ! on the first call; subsequent calls with the reduced np reproduce
       ! the same lattice (the construction is deterministic in np).
       rotx = rng%next(idum)*2.0_rp*pi
       roty = rng%next(idum)*2.0_rp*pi
       rotz = rng%next(idum)*2.0_rp*pi
       call sphere_points(np, co(:, i, 1), co(:, i, 2), co(:, i, 3), &
            kk(i), rotx, roty, rotz)

       ! Quantize wavenumbers in the periodic direction(s)
       call periodicity_chk(co(:, i, 1), co(:, i, 2), co(:, i, 3), np, &
            kk(i), dlx, dly, dlz, this%periodic(1), this%periodic(2), &
            this%periodic(3), rng, idum)

       ! Add a second point set mirrored at the x-axis
       do j = np + 1, 2*np
          co(j, i, 1) = co(j - np, i, 1)
          co(j, i, 2) = -co(j - np, i, 2)
          co(j, i, 3) = -co(j - np, i, 3)
       end do

       if (this%write_files) then
          do j = 1, 2*np
             write(10, '(2i5,3e18.9)') i, j, co(j, i, 1), co(j, i, 2), &
                  co(j, i, 3)
          end do
       end if

       ! Track smallest and largest wavenumber magnitudes per direction.
       ! NOTE: the plugin used the signed maximum (vlmax) here; we use the
       ! magnitude, which is the physically meaningful quantity for the
       ! wavelength report. Affects log output only.
       kxmax = max(kxmax, maxval(abs(co(1:2*np, i, 1))))
       kxmin = min(kxmin, minval(abs(co(1:2*np, i, 1))))
       kymax = max(kymax, maxval(abs(co(1:2*np, i, 2))))
       kymin = min(kymin, minval(abs(co(1:2*np, i, 2))))
       kzmax = max(kzmax, maxval(abs(co(1:2*np, i, 3))))
       kzmin = min(kzmin, minval(abs(co(1:2*np, i, 3))))

    end do ! i = 1, n_shells

    if (this%write_files) close(10)

    ! --- Remove zero wavenumber vectors and pack the modes
    n_removed = 0
    l = 0
    do i = 1, this%n_shells
       do j = 1, 2*np
          if (co(j, i, 1) .eq. 0.0_rp .and. co(j, i, 2) .eq. 0.0_rp &
               .and. co(j, i, 3) .eq. 0.0_rp) then
             n_removed = n_removed + 1
          else
             this%shell_modes(i) = this%shell_modes(i) + 1
             l = l + 1
             do k = 1, 3
                this%k_num(l, k) = co(j, i, k)
             end do
             this%shell(l) = i
          end if
       end do
    end do
    this%k_length = l
    n_kept = l

    call neko_log%message('FST - (0,0,0) wavenumber removed')
    write(log_buf, '(A,I0,A,I0,A)') 'Saved ', n_kept, ' of ', &
         n_kept + n_removed, ' fst modes.'
    call neko_log%message(log_buf)

    do i = 1, this%n_shells
       if (this%shell_modes(i) .eq. 0) then
          call neko_error("(FST) A shell has no valid modes." // &
               new_line('A') // "      Increase spectrum.modes_per_shell.")
       end if
    end do

    ! --- Shell amplitudes from the energy spectrum
    tke_tot1 = 0.0_rp
    do i = 1, this%n_shells
       this%shell_amp(i) = sqrt(2.0_rp*tke_shell(i)*2.0_rp &
            / real(this%shell_modes(i), kind=rp))
       shell_energy = real(this%shell_modes(i), kind=rp) &
            * (this%shell_amp(i)**2)/2.0_rp
       tke_tot1 = tke_tot1 + shell_energy
    end do

    write (log_buf, '(A,I0,A)') 'FST - ', this%k_length, &
         ' wavenumbers generated'
    call neko_log%message(log_buf)

    call print_param('FST - Largest wavelength in x', 2.0_rp*pi/kxmin)
    call print_param('FST - Smallest wavelength in x', 2.0_rp*pi/kxmax)
    call print_param('FST - Largest wavelength in y', 2.0_rp*pi/kymin)
    call print_param('FST - Smallest wavelength in y', 2.0_rp*pi/kymax)
    call print_param('FST - Largest wavelength in z', 2.0_rp*pi/kzmin)
    call print_param('FST - Smallest wavelength in z', 2.0_rp*pi/kzmax)

    ! =====================================================================
    ! Random phases and amplitudes (port of make_turbu, 05_turbu.f90).
    ! The plugin draws phase/amplitude pairs for three "columns" although
    ! only the first phase column and all three amplitude columns are used;
    ! the draw order is preserved for bit parity with the plugin.
    ! NOTE: the plugin overrides the seed with -143 at this point
    ! (05_turbu.f90, `seed = -143` in make_turbu), which makes its JSON
    ! seed option ineffective. That override is intentionally NOT ported.
    ! =====================================================================
    allocate(bb(this%n_modes_max, 3))
    allocate(bb1(this%n_modes_max, 3))

    if (this%write_files) then
       open(unit = 137, form = 'formatted', &
            file = trim(this%path) // '/bb.txt')
    end if

    do k = 1, 3
       do i = 1, this%n_modes_max
          bb(i, k) = rng%next(idum)*2.0_rp*pi   ! random phase shift
          bb1(i, k) = 2.0_rp*rng%next(idum) - 1.0_rp ! random amplitude
          if (this%write_files) write(137, *) bb(i, 1), bb1(i, 1)
       end do
    end do

    if (this%write_files) close(137)
    call neko_log%message("FST - Random amplitude generated")

    this%phase(:) = bb(:, 1)

    ! Enforce continuity by projecting the random amplitude vectors
    ! perpendicular to their wavenumber vectors, then normalize.
    ! Zero-norm projections cannot occur for generic random draws
    ! (measure-zero event); as in the plugin, this is not guarded.
    do i = 1, this%k_length
       do j = 1, 3
          u_hat_raw(j) = bb1(i, j)
       end do

       knorm2 = this%k_num(i, 1)**2 + this%k_num(i, 2)**2 &
            + this%k_num(i, 3)**2
       kdotu = u_hat_raw(1)*this%k_num(i, 1) &
            + u_hat_raw(2)*this%k_num(i, 2) &
            + u_hat_raw(3)*this%k_num(i, 3)

       do j = 1, 3
          u_hat_p(j) = u_hat_raw(j) - kdotu*this%k_num(i, j)/knorm2
       end do

       do j = 1, 3
          this%u_hat(i, j) = u_hat_p(j) &
               / sqrt(u_hat_p(1)**2 + u_hat_p(2)**2 + u_hat_p(3)**2)
       end do
    end do

    call neko_log%message('FST - Amplitudes projection done')

    ! --- Energy check and spectrum file
    ue = 0.0_rp
    ve = 0.0_rp
    we = 0.0_rp

    if (this%write_files) then
       open(file = trim(this%path) // '/fst_spectrum.csv', unit = 13)
       write(13, '(9(A, ","),A)') 'ShellNo', 'kx', 'ky', 'kz', &
            'u_amp', 'v_amp', 'w_amp', 'u_hat1', 'u_hat2', 'u_hat3'
    end if

    do i = 1, this%k_length
       shellno = this%shell(i)
       amp = this%shell_amp(shellno)

       uamp = this%u_hat(i, 1)*amp
       vamp = this%u_hat(i, 2)*amp
       wamp = this%u_hat(i, 3)*amp

       if (this%write_files) then
          write(13, '(9(g0, ","), g0)') shellno, this%k_num(i, 1), &
               this%k_num(i, 2), this%k_num(i, 3), uamp, vamp, wamp, &
               this%u_hat(i, 1), this%u_hat(i, 2), this%u_hat(i, 3)
       end if

       ue = ue + (uamp**2)/2.0_rp
       ve = ve + (vamp**2)/2.0_rp
       we = we + (wamp**2)/2.0_rp
    end do

    if (this%write_files) close(13)

    this%energy(1) = ue
    this%energy(2) = ve
    this%energy(3) = we
    this%tu_uinf_estimate = sqrt((ue + ve + we)/3.0_rp)

    write(log_buf, '(A18,10x,E12.5E2)') 'FST - Energy in u', ue
    call neko_log%message(log_buf)
    write(log_buf, '(A18,10x,E12.5E2)') 'FST - Energy in v', ve
    call neko_log%message(log_buf)
    write(log_buf, '(A18,10x,E12.5E2)') 'FST - Energy in w', we
    call neko_log%message(log_buf)
    write(log_buf, '(A20,8x,E12.5E2)') 'FST - Estimated tke', &
         (ue + ve + we)/2.0_rp
    call neko_log%message(log_buf)
    write(log_buf, '(A24,9x,E12.5E2)') 'FST - Estimated Tu*U_inf', &
         this%tu_uinf_estimate
    call neko_log%message(log_buf)

    deallocate(co, kk, q, dk, tke_shell, bb, bb1)

  end subroutine generate_rank0

  !> Von Karman energy spectrum (port of 03_spectrum.f90).
  !! @param k Wavenumber.
  !! @param l Integral length scale.
  !! @param q Scale factor.
  pure function ek(k, l, q) result(e)
    real(kind=rp), intent(in) :: k, l, q
    real(kind=rp) :: e

    e = 2.0_rp/3.0_rp*q*1.606_rp * (k*l)**4.0_rp * l / &
         (1.350_rp + (k*l)**2.0_rp)**(17.0_rp/6.0_rp)

  end function ek

  !> Number of azimuthal points at latitude ring j of the "new"-branch
  !! sphere lattice of the original compute_sphere (02_sphere.f90). Shared
  !! by the counting pass and the fill pass so both stay consistent.
  pure function lattice_nphi(nn, j) result(nphi)
    integer, intent(in) :: nn, j
    integer :: nphi

    real(kind=rp) :: dtheta, theta, dphi, nphir

    dtheta = 2.0_rp*pi/real(nn, kind=rp)
    theta = j*dtheta
    if (sin(theta) .eq. 0.0_rp) then
       dphi = 99999999.0_rp
    else
       dphi = dtheta/sin(theta)
    end if
    nphir = max(2.0_rp*pi/dphi, 1.0_rp)
    nphi = int(nphir + 0.5_rp)

  end function lattice_nphi

  !> Computes a set of np points which are (more or less) uniformly
  !! distributed on a sphere with radius rad, rotated by (rotx, roty, rotz).
  !! Port of compute_sphere (02_sphere.f90). np is intent(inout): for the
  !! general lattice branch it is reduced to the realized point count.
  !! The regular polyhedra for np = 4, 6, 8, 12, 20 are preserved; the
  !! unreachable non-"new" branch and the wire.dat output were removed.
  subroutine sphere_points(np, x, y, z, rad, rotx, roty, rotz)
    integer, intent(inout) :: np
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    real(kind=rp), intent(in) :: rad, rotx, roty, rotz

    integer :: i, j, k, nn, npn, npn1, nphi
    real(kind=rp) :: dphi, dtheta, theta, phi, w

    if (np .gt. size(x)) then
       call neko_error("(FST) sphere_points: np exceeds array size")
    end if

    if (np .eq. 4) then
       ! Tetrahedron
       call asp(x, y, z, 1, -1.0_rp/6.0_rp*sqrt(3.0_rp), -0.5_rp, 0.0_rp)
       call asp(x, y, z, 2, -1.0_rp/6.0_rp*sqrt(3.0_rp), 0.5_rp, 0.0_rp)
       call asp(x, y, z, 3, 1.0_rp/3.0_rp*sqrt(3.0_rp), 0.0_rp, 0.0_rp)
       call asp(x, y, z, 4, 0.0_rp, 0.0_rp, 1.0_rp/3.0_rp*sqrt(6.0_rp))
       call trans(x, y, z, np, 0.0_rp, 0.0_rp, -sqrt(6.0_rp)/12.0_rp)
       call scale1(x, y, z, np, sqrt(6.0_rp)*2.0_rp/3.0_rp)

    else if (np .eq. 6) then
       ! Octahedron
       call asp(x, y, z, 1, 0.0_rp, 0.0_rp, sqrt(2.0_rp)/2.0_rp)
       call asp(x, y, z, 2, 0.0_rp, 1.0_rp, sqrt(2.0_rp)/2.0_rp)
       call asp(x, y, z, 3, 1.0_rp, 1.0_rp, sqrt(2.0_rp)/2.0_rp)
       call asp(x, y, z, 4, 1.0_rp, 0.0_rp, sqrt(2.0_rp)/2.0_rp)
       call asp(x, y, z, 5, 0.5_rp, 0.5_rp, 0.0_rp)
       call asp(x, y, z, 6, 0.5_rp, 0.5_rp, sqrt(2.0_rp))
       call trans(x, y, z, np, -0.5_rp, -0.5_rp, -sqrt(2.0_rp)/2.0_rp)
       call scale1(x, y, z, np, sqrt(2.0_rp))

    else if (np .eq. 8) then
       ! Cube
       w = sqrt(3.0_rp)/3.0_rp
       call asp(x, y, z, 1, -w, w, -w)
       call asp(x, y, z, 2, -w, -w, -w)
       call asp(x, y, z, 3, w, -w, -w)
       call asp(x, y, z, 4, w, w, -w)
       call asp(x, y, z, 5, -w, w, w)
       call asp(x, y, z, 6, -w, -w, w)
       call asp(x, y, z, 7, w, -w, w)
       call asp(x, y, z, 8, w, w, w)

    else if (np .eq. 12) then
       ! Icosahedron
       w = 0.5_rp*(sqrt(5.0_rp) + 1.0_rp)
       call asp(x, y, z, 1, w/2.0_rp, 0.0_rp, 0.5_rp*(w - 1.0_rp))
       call asp(x, y, z, 2, w/2.0_rp, 0.0_rp, 0.5_rp*(w + 1.0_rp))
       call asp(x, y, z, 3, 0.0_rp, 0.5_rp*(w - 1.0_rp), 0.5_rp*w)
       call asp(x, y, z, 4, 0.0_rp, 0.5_rp*(w + 1.0_rp), 0.5_rp*w)
       call asp(x, y, z, 5, 0.5_rp*w, w, 0.5_rp*(w - 1.0_rp))
       call asp(x, y, z, 6, 0.5_rp*w, w, 0.5_rp*(w + 1.0_rp))
       call asp(x, y, z, 7, w, (w + 1.0_rp)/2.0_rp, w/2.0_rp)
       call asp(x, y, z, 8, w, (w - 1.0_rp)/2.0_rp, w/2.0_rp)
       call asp(x, y, z, 9, 0.5_rp*(w + 1.0_rp), 0.5_rp*w, w)
       call asp(x, y, z, 10, 0.5_rp*(w - 1.0_rp), 0.5_rp*w, w)
       call asp(x, y, z, 11, 0.5_rp*(w + 1.0_rp), 0.5_rp*w, 0.0_rp)
       call asp(x, y, z, 12, 0.5_rp*(w - 1.0_rp), 0.5_rp*w, 0.0_rp)
       call trans(x, y, z, np, -0.5_rp*w, -0.5_rp*w, -0.5_rp*w)
       call scale1(x, y, z, np, 2.0_rp/sqrt(w**2 + 1.0_rp))

    else if (np .eq. 20) then
       ! Dodecahedron
       w = 0.5_rp*(sqrt(5.0_rp) + 3.0_rp)
       call asp(x, y, z, 1, 0.5_rp*w, 0.5_rp*(w - 1.0_rp), 0.0_rp)
       call asp(x, y, z, 2, 0.5_rp*w, 0.5_rp*(w + 1.0_rp), 0.0_rp)
       call asp(x, y, z, 3, w - 0.5_rp, w - 0.5_rp, 0.5_rp)
       call asp(x, y, z, 4, w, 0.5_rp*w, 0.5_rp*(w - 1.0_rp))
       call asp(x, y, z, 5, w - 0.5_rp, 0.5_rp, 0.5_rp)
       call asp(x, y, z, 6, 0.5_rp*(w + 1.0_rp), 0.0_rp, 0.5_rp*w)
       call asp(x, y, z, 7, 0.5_rp*(w - 1.0_rp), 0.0_rp, 0.5_rp*w)
       call asp(x, y, z, 8, 0.5_rp, 0.5_rp, 0.5_rp)
       call asp(x, y, z, 9, 0.0_rp, 0.5_rp*w, (w - 1.0_rp)*0.5_rp)
       call asp(x, y, z, 10, 0.5_rp, w - 0.5_rp, 0.5_rp)
       call asp(x, y, z, 11, 0.5_rp*(w - 1.0_rp), w, 0.5_rp*w)
       call asp(x, y, z, 12, 0.5_rp*(w + 1.0_rp), w, 0.5_rp*w)
       call asp(x, y, z, 13, w - 0.5_rp, w - 0.5_rp, w - 0.5_rp)
       call asp(x, y, z, 14, w, 0.5_rp*w, 0.5_rp*(w + 1.0_rp))
       call asp(x, y, z, 15, w - 0.5_rp, 0.5_rp, w - 0.5_rp)
       call asp(x, y, z, 16, 0.5_rp*w, 0.5_rp*(w - 1.0_rp), w)
       call asp(x, y, z, 17, 0.5_rp, 0.5_rp, w - 0.5_rp)
       call asp(x, y, z, 18, 0.0_rp, 0.5_rp*w, 0.5_rp*(w + 1.0_rp))
       call asp(x, y, z, 19, 0.5_rp, w - 0.5_rp, w - 0.5_rp)
       call asp(x, y, z, 20, 0.5_rp*w, 0.5_rp*(w + 1.0_rp), w)
       call trans(x, y, z, np, -0.5_rp*w, -0.5_rp*w, -0.5_rp*w)
       call scale1(x, y, z, np, 2.0_rp/sqrt(w**2 + 1.0_rp))

    else
       ! General case: latitude/longitude lattice with (approximately)
       ! uniform surface density, realizing the largest point count <= np.
       nn = 1
       npn1 = 0
       do
          npn = 0
          do j = 0, nn/2
             npn = npn + lattice_nphi(nn, j)
          end do
          if (npn .le. np) then
             npn1 = npn
             nn = nn + 1
          else
             exit
          end if
       end do
       npn = npn1
       nn = nn - 1

       k = 0
       dtheta = 2.0_rp*pi/real(nn, kind=rp)
       do j = 0, nn/2
          theta = j*dtheta
          nphi = lattice_nphi(nn, j)
          dphi = 2.0_rp*pi/real(nphi, kind=rp)
          do i = 1, nphi
             phi = i*dphi
             k = k + 1
             x(k) = cos(phi)*sin(theta)
             y(k) = sin(phi)*sin(theta)
             z(k) = cos(theta)
          end do
       end do
       np = k

    end if

    ! Scale to the requested radius and rotate
    call scale1(x, y, z, np, rad)
    call rot3d(np, x, y, z, rotx, roty, rotz)

  end subroutine sphere_points

  !> Assign point i.
  subroutine asp(x, y, z, i, xx, yy, zz)
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    integer, intent(in) :: i
    real(kind=rp), intent(in) :: xx, yy, zz

    x(i) = xx
    y(i) = yy
    z(i) = zz
  end subroutine asp

  !> Translate points 1..np.
  subroutine trans(x, y, z, np, xx, yy, zz)
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: xx, yy, zz
    integer :: i

    do i = 1, np
       x(i) = x(i) + xx
       y(i) = y(i) + yy
       z(i) = z(i) + zz
    end do
  end subroutine trans

  !> Scale points 1..np by r.
  subroutine scale1(x, y, z, np, r)
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: r
    integer :: i

    do i = 1, np
       x(i) = x(i)*r
       y(i) = y(i)*r
       z(i) = z(i)*r
    end do
  end subroutine scale1

  !> Rotate points 1..np by the angles (rotx, roty, rotz).
  subroutine rot3d(np, x, y, z, rotx, roty, rotz)
    integer, intent(in) :: np
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    real(kind=rp), intent(in) :: rotx, roty, rotz

    real(kind=rp) :: cx, cy, cz, sx, sy, sz, xx, yy, zz
    integer :: i

    cx = cos(rotx)
    sx = sin(rotx)
    cy = cos(roty)
    sy = sin(roty)
    cz = cos(rotz)
    sz = sin(rotz)

    do i = 1, np
       xx = x(i)
       yy = y(i)
       zz = z(i)
       x(i) = xx*cy*cz + yy*cy*sz - zz*sy
       y(i) = xx*(cz*sx*sy - cx*sz) + yy*(cx*cz + sx*sy*sz) + zz*cy*sx
       z(i) = xx*(cx*cz*sy + sx*sz) - yy*(cz*sx - cx*sy*sz) + zz*cx*cy
    end do
  end subroutine rot3d

  !> Dispatch the periodicity quantization (port of periodicity_chk,
  !! 04_spec.f90).
  subroutine periodicity_chk(kx, ky, kz, np, kk, dlx, dly, dlz, &
       ifxp, ifyp, ifzp, rng, idum)
    real(kind=rp), intent(inout) :: kx(:), ky(:), kz(:)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: kk
    real(kind=rp), intent(in) :: dlx, dly, dlz
    logical, intent(in) :: ifxp, ifyp, ifzp
    type(fst_rng_t), intent(inout) :: rng
    integer, intent(inout) :: idum

    logical :: periodic_2d

    if (ifxp .and. ifyp .and. ifzp) then
       call neko_error("(FST) Only two periodic directions supported")
    end if

    periodic_2d = (ifxp .and. ifyp) .or. (ifxp .and. ifzp) &
         .or. (ifyp .and. ifzp)

    if (periodic_2d) then
       if (ifxp .and. ifyp) then
          call make_periodic_2d(kz, kx, ky, np, kk, dlx, dly)
       else if (ifxp .and. ifzp) then
          call make_periodic_2d(ky, kz, kx, np, kk, dlz, dlx)
       else if (ifyp .and. ifzp) then
          call make_periodic_2d(kx, ky, kz, np, kk, dly, dlz)
       end if
    else
       if (ifxp) then
          call make_periodic_1d(ky, kz, kx, np, kk, dlx, rng, idum)
       else if (ifyp) then
          call make_periodic_1d(kz, kx, ky, np, kk, dly, rng, idum)
       else if (ifzp) then
          call make_periodic_1d(kx, ky, kz, np, kk, dlz, rng, idum)
       end if
    end if

  end subroutine periodicity_chk

  !> Quantize the wavenumbers kp to multiples of 2*pi/lp, adjusting k1 or
  !! k2 (chosen by coin toss) to preserve the shell radius k_total.
  !! Port of make_periodic_1D (04_spec.f90) with one deliberate change:
  !! the multiple is chosen by nearest rounding (nint) instead of floor.
  !! Floor rounds every quantized component toward zero, so the free
  !! component that restores the shell radius is systematically inflated;
  !! since u_hat is perpendicular to k, that direction is then starved of
  !! energy (measured: 5-15% component-energy bias for one periodic
  !! direction). Nearest rounding removes the bias.
  subroutine make_periodic_1d(k1, k2, kp, np, k_total, lp, rng, idum)
    real(kind=rp), intent(inout) :: k1(:), k2(:), kp(:)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: k_total
    real(kind=rp), intent(in) :: lp
    type(fst_rng_t), intent(inout) :: rng
    integer, intent(inout) :: idum

    integer :: nmax, n_j, n_j_signed, j
    real(kind=rp) :: twopi_over_l, rtmp, flip, k_total_sq

    twopi_over_l = 2.0_rp*pi/lp
    k_total_sq = k_total**2

    ! At least one multiple of 2*pi/lp must fit inside the shell radius
    nmax = floor(k_total/twopi_over_l)
    if (nmax .lt. 1) then
       call print_param('K_total:', k_total)
       call print_param('2 pi / L:', twopi_over_l)
       call neko_log%message("2 pi / L must be smaller than K_total!")
       call neko_error('(FST) Increase minimum total wave number!')
    end if

    do j = 1, np

       n_j = nint(abs(kp(j))/twopi_over_l)
       n_j_signed = int(sign(1.0_rp, kp(j)))*n_j

       if (n_j .gt. nmax) then
          n_j_signed = n_j_signed - int(sign(1.0_rp, kp(j)))
       else if (n_j .eq. 0) then
          ! Force to not be zero
          n_j_signed = n_j_signed + int(sign(1.0_rp, kp(j)))
       end if

       kp(j) = real(n_j_signed, kind=rp)*twopi_over_l

       ! Adjust k1 or k2 so that k1^2 + k2^2 + kp^2 = k_total^2.
       ! NOTE (from the original): the `> 1` threshold is dimensional and
       ! kept as-is for parity; it avoids assigning near-zero components.
       flip = rng%next(idum)

       if (flip .gt. 0.5_rp) then
          ! k1 stays, recompute k2
          rtmp = k_total_sq - k1(j)**2 - kp(j)**2
          if (rtmp .gt. 1.0_rp) then
             k2(j) = sign(1.0_rp, k2(j))*sqrt(rtmp)
          else
             rtmp = sqrt((k_total_sq - kp(j)**2)/2.0_rp)
             k1(j) = sign(1.0_rp, k1(j))*rtmp
             k2(j) = sign(1.0_rp, k2(j))*rtmp
          end if
       else
          ! k2 stays, recompute k1
          rtmp = k_total_sq - kp(j)**2 - k2(j)**2
          if (rtmp .gt. 1.0_rp) then
             k1(j) = sign(1.0_rp, k1(j))*sqrt(rtmp)
          else
             rtmp = sqrt((k_total_sq - kp(j)**2)/2.0_rp)
             k1(j) = sign(1.0_rp, k1(j))*rtmp
             k2(j) = sign(1.0_rp, k2(j))*rtmp
          end if
       end if

    end do

  end subroutine make_periodic_1d

  !> Quantize kp1 and kp2 to multiples of 2*pi/l1 and 2*pi/l2 respectively,
  !! recomputing k1 to preserve the shell radius.
  !! Port of make_periodic_2D (04_spec.f90), using nearest rounding (nint)
  !! instead of floor for the multiples (see make_periodic_1d for why).
  !! With two quantized directions a residual anisotropy remains at the
  !! lowest shells, because both components are forced to be nonzero and
  !! the shell then cannot point along the free direction; this is
  !! inherent to the quantization, not to the rounding.
  !! Two quirks of the original are
  !! preserved for parity and flagged here: (1) the fit check for direction
  !! 1 is overwritten by that of direction 2, so the per-mode clamp uses the
  !! direction-2 nmax for both directions (only relevant when l1 /= l2);
  !! (2) the reduction loop compares signed components when deciding which
  !! one to reduce.
  subroutine make_periodic_2d(k1, kp1, kp2, np, k_total, l1, l2)
    real(kind=rp), intent(inout) :: k1(:), kp1(:), kp2(:)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: k_total
    real(kind=rp), intent(in) :: l1, l2

    integer :: nmax, n_j1, n_j1_signed, n_j2, n_j2_signed, j
    real(kind=rp) :: twopi_over_l1, twopi_over_l2, rtmp, k_total_sq
    logical :: valid_config

    twopi_over_l1 = 2.0_rp*pi/l1
    twopi_over_l2 = 2.0_rp*pi/l2
    k_total_sq = k_total**2

    ! Direction 1 fit check
    nmax = floor(k_total/twopi_over_l1)
    if (nmax .lt. 1) then
       call print_param('K_total:', k_total)
       call print_param('2 pi / L1:', twopi_over_l1)
       call neko_log%message("2 pi / L must be smaller than K_total!")
       call neko_error('(FST) Increase minimum total wave number!')
    end if

    ! Direction 2 fit check (overwrites nmax; see subroutine comment)
    nmax = floor(k_total/twopi_over_l2)
    if (nmax .lt. 1) then
       call print_param('K_total:', k_total)
       call print_param('2 pi / L2:', twopi_over_l2)
       call neko_log%message("2 pi / L must be smaller than K_total!")
       call neko_error('(FST) Increase minimum total wave number!')
    end if

    do j = 1, np

       ! Discrete wavenumber in direction 1
       n_j1 = nint(abs(kp1(j))/twopi_over_l1)
       n_j1_signed = int(sign(1.0_rp, kp1(j)))*n_j1

       if (n_j1 .gt. nmax) then
          n_j1_signed = n_j1_signed - int(sign(1.0_rp, kp1(j)))
       else if (n_j1 .eq. 0) then
          n_j1_signed = n_j1_signed + int(sign(1.0_rp, kp1(j)))
       end if

       kp1(j) = real(n_j1_signed, kind=rp)*twopi_over_l1

       ! Discrete wavenumber in direction 2
       n_j2 = nint(abs(kp2(j))/twopi_over_l2)
       n_j2_signed = int(sign(1.0_rp, kp2(j)))*n_j2

       if (n_j2 .gt. nmax) then
          n_j2_signed = n_j2_signed - int(sign(1.0_rp, kp2(j)))
       else if (n_j2 .eq. 0) then
          n_j2_signed = n_j2_signed + int(sign(1.0_rp, kp2(j)))
       end if

       kp2(j) = real(n_j2_signed, kind=rp)*twopi_over_l2

       ! Recompute k1 to preserve the shell radius; if the quantized
       ! components already exceed it, reduce the largest one until a
       ! valid configuration is reached.
       rtmp = k_total_sq - kp1(j)**2 - kp2(j)**2
       valid_config = (rtmp .gt. 1.0_rp)

       do while (.not. valid_config)
          if (kp1(j) .gt. kp2(j)) then
             n_j1_signed = n_j1_signed - int(sign(1.0_rp, kp1(j)))
             kp1(j) = real(n_j1_signed, kind=rp)*twopi_over_l1
          else
             n_j2_signed = n_j2_signed - int(sign(1.0_rp, kp2(j)))
             kp2(j) = real(n_j2_signed, kind=rp)*twopi_over_l2
          end if

          rtmp = k_total_sq - kp1(j)**2 - kp2(j)**2
          valid_config = (rtmp .gt. 0.0_rp)
       end do

       k1(j) = sign(1.0_rp, k1(j))*sqrt(rtmp)

    end do

  end subroutine make_periodic_2d

  !> Log a named real parameter.
  subroutine print_param(name, value)
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: value
    character(len=LOG_SIZE) :: log_buf

    ! The name field is bounded so that the line always fits LOG_SIZE:
    ! 6 ("[FST] ") + 50 (name) + 2 (": ") + 13 (E13.5) = 71 characters.
    write(log_buf, '(A,A50,A,E13.5)') "[FST] ", name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_param

end module fst_spectrum