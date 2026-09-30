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
!> Mode set for synthetic free-stream turbulence: wavenumbers on spherical
!! shells, divergence-free directions, random phases and amplitudes from a
!! von Karman spectrum. Generated on rank 0 and broadcast.
!!
!! Ported from the FST plugin of V. Baconnet, E. Kluesberg, P. Negi and
!! P. Schlatter (neko-plugins). For the same inputs it gives the plugin's
!! mode set, except for these deliberate changes:
!!  - all sizes and parameters are runtime inputs,
!!  - the seed is honoured (the plugin always used -143),
!!  - periodic wavenumbers are rounded to the nearest multiple instead of
!!    down, which removes a 5-15% energy bias between components, and
!!    are never reduced to zero,
!!  - unused plugin code is dropped.
module fst_spectrum
  use num_types, only : rp
  use math, only : pi
  use utils, only : neko_error
  use logger, only : neko_log, LOG_SIZE
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use mpi_f08, only : MPI_Bcast, MPI_INTEGER
  implicit none
  private

  !> Numerical Recipes ran2, kept identical to the plugin for parity.
  type :: fst_rng_t
     integer :: ir(97) = 0
     integer :: iy = 0
     integer :: iff = 0
   contains
     procedure, pass(this) :: next => fst_rng_next
  end type fst_rng_t

  type, public :: fst_spectrum_t
     !> Shells, requested points per shell, and array size 2*n_shells*npmax.
     integer :: n_shells = 0
     integer :: npmax = 0
     integer :: n_modes_max = 0
     !> Wavenumber band.
     real(kind=rp) :: k_start = 0.0_rp
     real(kind=rp) :: k_end = 0.0_rp
     !> Turbulence intensity (fraction), integral length scale and the
     !! speed it refers to.
     real(kind=rp) :: ti = 0.0_rp
     real(kind=rp) :: il = 0.0_rp
     real(kind=rp) :: u_ref = 0.0_rp
     !> Directions whose wavenumbers are quantized to 2*pi*n/L.
     logical :: periodic(3) = .false.
     integer :: seed = -143
     !> Write sphere.dat, bb.txt and fst_spectrum.csv on rank 0.
     logical :: write_files = .false.
     character(len=:), allocatable :: path

     !> Number of valid modes; the arrays below are valid up to k_length.
     integer :: k_length = 0
     real(kind=rp), allocatable :: k_num(:,:)
     real(kind=rp), allocatable :: u_hat(:,:)
     real(kind=rp), allocatable :: phase(:)
     integer, allocatable :: shell(:)
     real(kind=rp), allocatable :: shell_amp(:)
     integer, allocatable :: shell_modes(:)

     !> Realized Tu*U and component energies, valid on rank 0 only.
     real(kind=rp) :: tu_uinf_estimate = 0.0_rp
     real(kind=rp) :: energy(3) = 0.0_rp
   contains
     procedure, pass(this) :: init => fst_spectrum_init
     procedure, pass(this) :: generate => fst_spectrum_generate
     procedure, pass(this) :: free => fst_spectrum_free
  end type fst_spectrum_t

contains

  !> Next uniform number in [0, 1). A negative idum reseeds.
  function fst_rng_next(this, idum) result(r)
    class(fst_rng_t), intent(inout) :: this
    integer, intent(inout) :: idum
    real(kind=rp) :: r

    integer, parameter :: m = 714025, ia = 1366, ic = 150889
    real, parameter :: rm = 1./m
    integer :: j

    if (idum .lt. 0 .or. this%iff .eq. 0) then
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

  !> Check the parameters and allocate the mode arrays.
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

    this%k_num = 0.0_rp
    this%u_hat = 0.0_rp
    this%phase = 0.0_rp
    this%shell = 0
    this%shell_amp = 0.0_rp
    this%shell_modes = 0

  end subroutine fst_spectrum_init

  !> Destructor.
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

  !> Generate on rank 0 and broadcast.
  !! @param lx, ly, lz Domain lengths, used in periodic directions.
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

  !> Shells, sphere points, amplitudes and continuity projection.
  !! The order of random draws matches the plugin.
  subroutine generate_rank0(this, dlx, dly, dlz)
    class(fst_spectrum_t), intent(inout) :: this
    real(kind=rp), intent(in) :: dlx, dly, dlz

    type(fst_rng_t) :: rng
    integer :: idum, np
    real(kind=rp), allocatable :: co(:,:,:), tke_shell(:)

    idum = this%seed
    allocate(co(2*this%npmax, this%n_shells, 3))
    allocate(tke_shell(this%n_shells))
    co = 0.0_rp

    call print_param('integral length scale', this%il)
    call shell_energies(this, tke_shell)
    call shell_points(this, dlx, dly, dlz, rng, idum, co, np)
    call pack_points(this, co, np)
    call shell_amplitudes(this, tke_shell)
    call log_wavelengths(co, np)
    call random_directions(this, rng, idum)
    call energy_check(this)

    deallocate(co, tke_shell)

  end subroutine generate_rank0

  !> Radius of shell i.
  pure function shell_radius(this, i) result(kk)
    class(fst_spectrum_t), intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp) :: kk, k2

    k2 = (this%k_start + (i-1)*(this%k_end - this%k_start) &
         / real(this%n_shells - 1, kind=rp))**2
    kk = sqrt(k2)

  end function shell_radius

  !> Energy of each shell, scaled so the band holds 3/2 (Ti U)^2.
  subroutine shell_energies(this, tke_shell)
    class(fst_spectrum_t), intent(in) :: this
    real(kind=rp), intent(out) :: tke_shell(:)

    real(kind=rp) :: dkint, dk, tke_tot, tke_band, tke_scaled, q
    integer :: i, ndk
    character(len=LOG_SIZE) :: log_buf

    tke_scaled = 1.5_rp * (this%ti * this%u_ref)**2

    ! Fine integral over the band, for the log only
    ndk = 5000
    dkint = (this%k_end - this%k_start)/real(ndk, kind=rp)
    tke_band = ek(this%k_start, this%il, 1.0_rp) &
         + ek(this%k_end, this%il, 1.0_rp)
    do i = 1, ndk - 1
       tke_band = tke_band + ek(this%k_start + i*dkint, this%il, 1.0_rp)
    end do
    tke_band = tke_band*dkint
    call print_param('FST - integrated energy in spectrum', tke_band)

    dkint = (this%k_end - this%k_start)/real(this%n_shells - 1, kind=rp)
    tke_tot = 0.0_rp
    do i = 1, this%n_shells
       tke_tot = tke_tot + ek(this%k_start + (i-1)*dkint, this%il, 1.0_rp)
    end do
    tke_tot = tke_tot*dkint
    write (log_buf, '(A,I0,A,E13.5)') 'FST - discretized on ', &
         this%n_shells, ' shells : ', tke_tot
    call neko_log%message(log_buf)
    call print_param("Truncated TKE", tke_scaled/tke_tot)

    do i = 1, this%n_shells
       dk = (this%k_end - this%k_start)/real(this%n_shells - 1, kind=rp)
       q = ek(shell_radius(this, i), this%il, tke_scaled/tke_tot)
       tke_shell(i) = q*dk
    end do

  end subroutine shell_energies

  !> Randomly rotated point set on every shell, quantized in periodic
  !! directions and mirrored about the x axis. np returns the points per
  !! shell actually used, before mirroring.
  subroutine shell_points(this, dlx, dly, dlz, rng, idum, co, np)
    class(fst_spectrum_t), intent(in) :: this
    real(kind=rp), intent(in) :: dlx, dly, dlz
    type(fst_rng_t), intent(inout) :: rng
    integer, intent(inout) :: idum
    real(kind=rp), intent(inout) :: co(:,:,:)
    integer, intent(out) :: np

    real(kind=rp) :: kk, rotx, roty, rotz
    integer :: i, j

    np = this%npmax

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

    do i = 1, this%n_shells
       kk = shell_radius(this, i)

       rotx = rng%next(idum)*2.0_rp*pi
       roty = rng%next(idum)*2.0_rp*pi
       rotz = rng%next(idum)*2.0_rp*pi
       call sphere_points(np, co(:, i, 1), co(:, i, 2), co(:, i, 3), &
            kk, rotx, roty, rotz)

       call periodicity_chk(co(:, i, 1), co(:, i, 2), co(:, i, 3), np, &
            kk, dlx, dly, dlz, this%periodic(1), this%periodic(2), &
            this%periodic(3), rng, idum)

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
    end do

    if (this%write_files) close(10)

  end subroutine shell_points

  !> Store the nonzero points as modes, with their shell index.
  subroutine pack_points(this, co, np)
    class(fst_spectrum_t), intent(inout) :: this
    real(kind=rp), intent(in) :: co(:,:,:)
    integer, intent(in) :: np

    integer :: i, j, k, l, n_removed
    character(len=LOG_SIZE) :: log_buf

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

    call neko_log%message('FST - (0,0,0) wavenumber removed')
    write(log_buf, '(A,I0,A,I0,A)') 'Saved ', l, ' of ', &
         l + n_removed, ' fst modes.'
    call neko_log%message(log_buf)

    do i = 1, this%n_shells
       if (this%shell_modes(i) .eq. 0) then
          call neko_error("(FST) A shell has no valid modes." // &
               new_line('A') // "      Increase spectrum.modes_per_shell.")
       end if
    end do

  end subroutine pack_points

  !> Amplitude of each shell, shared equally by its modes.
  subroutine shell_amplitudes(this, tke_shell)
    class(fst_spectrum_t), intent(inout) :: this
    real(kind=rp), intent(in) :: tke_shell(:)

    integer :: i
    character(len=LOG_SIZE) :: log_buf

    do i = 1, this%n_shells
       this%shell_amp(i) = sqrt(2.0_rp*tke_shell(i)*2.0_rp &
            / real(this%shell_modes(i), kind=rp))
    end do

    write (log_buf, '(A,I0,A)') 'FST - ', this%k_length, &
         ' wavenumbers generated'
    call neko_log%message(log_buf)

  end subroutine shell_amplitudes

  !> Log the longest and shortest wavelength in each direction.
  subroutine log_wavelengths(co, np)
    real(kind=rp), intent(in) :: co(:,:,:)
    integer, intent(in) :: np

    real(kind=rp) :: kmin(3), kmax(3)
    character(len=1), parameter :: dir_char(3) = ['x', 'y', 'z']
    integer :: d

    do d = 1, 3
       kmax(d) = max(1.0e-20_rp, maxval(abs(co(1:2*np, :, d))))
       kmin(d) = min(1.0e+20_rp, minval(abs(co(1:2*np, :, d))))
       call print_param('FST - Largest wavelength in ' // dir_char(d), &
            2.0_rp*pi/kmin(d))
       call print_param('FST - Smallest wavelength in ' // dir_char(d), &
            2.0_rp*pi/kmax(d))
    end do

  end subroutine log_wavelengths

  !> Random phases, and random directions projected normal to k so every
  !! mode is divergence free.
  subroutine random_directions(this, rng, idum)
    class(fst_spectrum_t), intent(inout) :: this
    type(fst_rng_t), intent(inout) :: rng
    integer, intent(inout) :: idum

    real(kind=rp), allocatable :: bb(:,:), bb1(:,:)
    real(kind=rp) :: u_hat_raw(3), u_hat_p(3), kdotu, knorm2
    integer :: i, j, k

    allocate(bb(this%n_modes_max, 3))
    allocate(bb1(this%n_modes_max, 3))

    if (this%write_files) then
       open(unit = 137, form = 'formatted', &
            file = trim(this%path) // '/bb.txt')
    end if

    ! Three columns are drawn though only one phase column is used, to keep
    ! the plugin's random sequence
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

    deallocate(bb, bb1)

  end subroutine random_directions

  !> Energy per component of the generated modes, and the realized Tu.
  subroutine energy_check(this)
    class(fst_spectrum_t), intent(inout) :: this

    real(kind=rp) :: ue, ve, we, uamp, vamp, wamp, amp
    integer :: i, shellno
    character(len=LOG_SIZE) :: log_buf

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

  end subroutine energy_check

  !> Von Karman spectrum; integrates to q over [0, inf).
  pure function ek(k, l, q) result(e)
    real(kind=rp), intent(in) :: k, l, q
    real(kind=rp) :: e

    e = 2.0_rp/3.0_rp*q*1.606_rp * (k*l)**4.0_rp * l / &
         (1.350_rp + (k*l)**2.0_rp)**(17.0_rp/6.0_rp)

  end function ek

  !> Points on latitude ring j of the sphere lattice.
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

  !> np points spread evenly on a sphere of radius rad, then rotated.
  !! Regular polyhedra for np = 4, 6, 8, 12, 20, otherwise a lattice that
  !! may use fewer points; np returns the number used.
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
       call asp(x, y, z, 1, -1.0_rp/6.0_rp*sqrt(3.0_rp), -0.5_rp, 0.0_rp)
       call asp(x, y, z, 2, -1.0_rp/6.0_rp*sqrt(3.0_rp), 0.5_rp, 0.0_rp)
       call asp(x, y, z, 3, 1.0_rp/3.0_rp*sqrt(3.0_rp), 0.0_rp, 0.0_rp)
       call asp(x, y, z, 4, 0.0_rp, 0.0_rp, 1.0_rp/3.0_rp*sqrt(6.0_rp))
       call trans(x, y, z, np, 0.0_rp, 0.0_rp, -sqrt(6.0_rp)/12.0_rp)
       call scale1(x, y, z, np, sqrt(6.0_rp)*2.0_rp/3.0_rp)

    else if (np .eq. 6) then
       call asp(x, y, z, 1, 0.0_rp, 0.0_rp, sqrt(2.0_rp)/2.0_rp)
       call asp(x, y, z, 2, 0.0_rp, 1.0_rp, sqrt(2.0_rp)/2.0_rp)
       call asp(x, y, z, 3, 1.0_rp, 1.0_rp, sqrt(2.0_rp)/2.0_rp)
       call asp(x, y, z, 4, 1.0_rp, 0.0_rp, sqrt(2.0_rp)/2.0_rp)
       call asp(x, y, z, 5, 0.5_rp, 0.5_rp, 0.0_rp)
       call asp(x, y, z, 6, 0.5_rp, 0.5_rp, sqrt(2.0_rp))
       call trans(x, y, z, np, -0.5_rp, -0.5_rp, -sqrt(2.0_rp)/2.0_rp)
       call scale1(x, y, z, np, sqrt(2.0_rp))

    else if (np .eq. 8) then
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

    call scale1(x, y, z, np, rad)
    call rot3d(np, x, y, z, rotx, roty, rotz)

  end subroutine sphere_points

  !> Set point i.
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

  !> Scale points 1..np.
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

  !> Rotate points 1..np.
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

  !> Quantize the periodic direction(s). At most two can be periodic.
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

  !> Snap kp to the nearest nonzero multiple of 2*pi/lp and adjust k1 or
  !! k2 (coin toss) so the shell radius is kept.
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
          n_j_signed = n_j_signed + int(sign(1.0_rp, kp(j)))
       end if

       kp(j) = real(n_j_signed, kind=rp)*twopi_over_l

       ! The > 1 test on rtmp below is dimensional; kept from the plugin
       flip = rng%next(idum)

       if (flip .gt. 0.5_rp) then
          rtmp = k_total_sq - k1(j)**2 - kp(j)**2
          if (rtmp .gt. 1.0_rp) then
             k2(j) = sign(1.0_rp, k2(j))*sqrt(rtmp)
          else
             rtmp = sqrt((k_total_sq - kp(j)**2)/2.0_rp)
             k1(j) = sign(1.0_rp, k1(j))*rtmp
             k2(j) = sign(1.0_rp, k2(j))*rtmp
          end if
       else
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

  !> Snap kp1 and kp2 to nonzero multiples of 2*pi/l1 and 2*pi/l2 and set
  !! k1 from the shell radius. The lowest shells stay slightly anisotropic
  !! since both quantized components must be nonzero.
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

    nmax = floor(k_total/twopi_over_l1)
    if (nmax .lt. 1) then
       call print_param('K_total:', k_total)
       call print_param('2 pi / L1:', twopi_over_l1)
       call neko_log%message("2 pi / L must be smaller than K_total!")
       call neko_error('(FST) Increase minimum total wave number!')
    end if

    ! As in the plugin, this overwrites nmax from direction 1, so the clamp
    ! below uses direction 2 for both (matters only when l1 /= l2)
    nmax = floor(k_total/twopi_over_l2)
    if (nmax .lt. 1) then
       call print_param('K_total:', k_total)
       call print_param('2 pi / L2:', twopi_over_l2)
       call neko_log%message("2 pi / L must be smaller than K_total!")
       call neko_error('(FST) Increase minimum total wave number!')
    end if

    do j = 1, np

       n_j1 = nint(abs(kp1(j))/twopi_over_l1)
       n_j1_signed = int(sign(1.0_rp, kp1(j)))*n_j1

       if (n_j1 .gt. nmax) then
          n_j1_signed = n_j1_signed - int(sign(1.0_rp, kp1(j)))
       else if (n_j1 .eq. 0) then
          n_j1_signed = n_j1_signed + int(sign(1.0_rp, kp1(j)))
       end if

       kp1(j) = real(n_j1_signed, kind=rp)*twopi_over_l1

       n_j2 = nint(abs(kp2(j))/twopi_over_l2)
       n_j2_signed = int(sign(1.0_rp, kp2(j)))*n_j2

       if (n_j2 .gt. nmax) then
          n_j2_signed = n_j2_signed - int(sign(1.0_rp, kp2(j)))
       else if (n_j2 .eq. 0) then
          n_j2_signed = n_j2_signed + int(sign(1.0_rp, kp2(j)))
       end if

       kp2(j) = real(n_j2_signed, kind=rp)*twopi_over_l2

       rtmp = k_total_sq - kp1(j)**2 - kp2(j)**2
       valid_config = (rtmp .gt. 1.0_rp)

       ! Step the larger component down until the pair fits the shell,
       ! but never to zero
       do while (.not. valid_config)
          if (abs(n_j1_signed) .gt. 1 .and. (abs(kp1(j)) .ge. abs(kp2(j)) &
               .or. abs(n_j2_signed) .eq. 1)) then
             n_j1_signed = n_j1_signed - int(sign(1.0_rp, kp1(j)))
             kp1(j) = real(n_j1_signed, kind=rp)*twopi_over_l1
          else if (abs(n_j2_signed) .gt. 1) then
             n_j2_signed = n_j2_signed - int(sign(1.0_rp, kp2(j)))
             kp2(j) = real(n_j2_signed, kind=rp)*twopi_over_l2
          else if (rtmp .gt. 0.0_rp) then
             exit
          else
             call neko_error("(FST) k_min is too small for two periodic " // &
                  "directions." // new_line('A') // &
                  "      It must satisfy k_min**2 > (2*pi/L1)**2 + " // &
                  "(2*pi/L2)**2.")
          end if

          rtmp = k_total_sq - kp1(j)**2 - kp2(j)**2
          valid_config = (rtmp .gt. 0.0_rp)
       end do

       k1(j) = sign(1.0_rp, k1(j))*sqrt(rtmp)

    end do

  end subroutine make_periodic_2d

  !> Log a named value.
  subroutine print_param(name, value)
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: value
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A,A50,A,E13.5)') "[FST] ", name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_param

end module fst_spectrum
