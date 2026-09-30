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
!> Implements `fst_source_term_t`: free-stream turbulence (FST) injected as
!! a volumetric fringe forcing
!!
!!    f_i = gain * ramp(t) * lambda(x,y,z) * (u_bf,i + u'_i(x,t) - u_i)
!!
!! where u' is a sum of Fourier modes sampled from a von Karman spectrum
!! (see `fst_spectrum`), convected as frozen turbulence at the constant
!! velocity vector U_c (phase k.(x - U_c t)), supporting oblique/cross
!! flow such as swept-wing configurations,
!! lambda is a per-direction SIMSON-type fringe with compact support and
!! ramp(t) is a linear ramp of length ramp_time starting at start_time.
!!
!! The perturbation and the fringe are evaluated at the CURRENT mesh
!! coordinates on every call, so the forcing is a purely Eulerian field and
!! remains consistent under ALE mesh deformation. The point-zone mask,
!! however, is built once from the initial mesh and follows material
!! points; the zone must therefore contain the fringe support [start, end]
!! at all times during deformation (make the zone generously larger than
!! the fringe; points with lambda = 0 are skipped at negligible cost).
!! The baseflow methods "initial_condition" and "field" also attach the
!! captured target velocity to material points; under mesh deformation use
!! "constant" unless this is the intended behavior.
!!
!! The compute kernels live in `fst_source_term_cpu` and
!! `fst_source_term_device`; this module holds the type, the JSON setup,
!! the diagnostics and the backend dispatch.
!!
!! JSON parameters (under case.fluid.source_terms, "type": "fst"):
!!   zone_name (string, mandatory)        Point zone of the forcing region.
!!   gain (real, mandatory)               Fringe gain [1/time].
!!   convection_velocity ([real x3], mandatory)
!!                                        Frozen-turbulence convection
!!                                        velocity U_c (nonzero magnitude).
!!                                        For flow aligned with x this is
!!                                        [U_inf, 0, 0]; for a swept wing,
!!                                        the oblique free-stream vector.
!!   turbulence_intensity (real, mandatory)
!!   integral_length_scale (real, mandatory)
!!   spectrum.n_shells (int, mandatory)
!!   spectrum.modes_per_shell (int, mandatory)
!!   spectrum.k_min (real, mandatory)
!!   spectrum.k_max (real, mandatory)
!!   baseflow.method (string, mandatory)  "initial_condition", "constant"
!!                                        (with baseflow.value = [u,v,w])
!!                                        or "field" (with
!!                                        baseflow.file_name, optional
!!                                        baseflow.interpolate,
!!                                        baseflow.mesh_file_name,
!!                                        baseflow.interpolation).
!!   fringe.x / fringe.y / fringe.z       Optional per direction. If absent
!!                                        the fringe is flat (lambda = 1) in
!!                                        that direction. If present, all of
!!                                        start, end, rise, fall are
!!                                        mandatory.
!!   periodic ([bool,bool,bool], optional, default [false,false,false])
!!                                        Wavenumber quantization to the
!!                                        full domain length per direction.
!!   start_time, end_time (optional)      Source term activity window.
!!   ramp_time (real, optional, 0)        Linear ramp length after
!!                                        start_time.
!!   seed (int, optional, -143)           RNG seed.
!!   write_files (bool, optional, false)  Write generation diagnostics.
!!   files_output_path (string, optional, "./fst_files")
!!   dump_fields (bool, optional, false)  Write fringe and u' snapshot to
!!                                        an fld file at the first compute.
!!   dump_file_name (string, optional, "fst_fields")
!!   validate_only (bool, optional, false) Stop after generation,
!!                                        diagnostics and field dump.
module fst_source_term
  use num_types, only : rp, sp
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, &
       json_get_subdict_or_empty
  use field, only : field_t
  use field_list, only : field_list_t
  use coefs, only : coef_t
  use source_term, only : source_term_t
  use fst_spectrum, only : fst_spectrum_t
  use point_zone, only : point_zone_t
  use point_zone_registry, only : neko_point_zone_registry
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  use fld_file_output, only : fld_file_output_t
  use import_field_utils, only : import_fields
  use time_state, only : time_state_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error, neko_warning, NEKO_FNAME_LEN
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, device_memcpy, device_free, &
       HOST_TO_DEVICE, DEVICE_TO_HOST
  use fst_source_term_cpu, only : fst_source_term_compute_cpu, &
       fst_source_term_preview_cpu
  use fst_source_term_device, only : fst_source_term_compute_device
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  use math, only : glmax, glmin, pi
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_MIN, MPI_MAX, &
       MPI_SUM, MPI_INTEGER
  implicit none
  private

  !> Free-stream turbulence fringe forcing source term.
  type, public, extends(source_term_t) :: fst_source_term_t
     !> The FST mode set.
     type(fst_spectrum_t) :: spectrum
     !> Point zone defining the forcing region.
     class(point_zone_t), pointer :: zone => null()
     !> Zone mask (local GLL indices, 1-based).
     integer, pointer :: mask(:) => null()
     !> Velocity fields.
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     !> Fringe gain [1/time].
     real(kind=rp) :: gain = 0.0_rp
     !> Frozen-turbulence convection velocity vector U_c.
     real(kind=rp) :: conv_vel(3) = 0.0_rp
     !> Reference speed |U_c| (energy scaling and diagnostics).
     real(kind=rp) :: u_ref = 0.0_rp
     !> Linear time ramp length after start_time.
     real(kind=rp) :: ramp_time = 0.0_rp
     !> Per-direction fringe: smooth (SIMSON step) or flat (lambda = 1).
     logical :: fringe_smooth(3) = .false.
     real(kind=rp) :: fringe_start(3) = 0.0_rp
     real(kind=rp) :: fringe_end(3) = 0.0_rp
     real(kind=rp) :: fringe_rise(3) = 0.0_rp
     real(kind=rp) :: fringe_fall(3) = 0.0_rp
     !> Baseflow method: "initial_condition", "constant" or "field".
     character(len=:), allocatable :: baseflow_method
     !> Number of modes (copy of spectrum%k_length, for the kernels).
     integer :: k_length = 0
     !> Mode data in a flat, kernel-friendly layout (length k_length):
     !! wavenumber components and amplitude-weighted direction vectors
     !! a_j(m) = u_hat(m,j) * shell_amp(shell(m)), so the kernels need no
     !! shell indirection. Shared by the CPU and device backends.
     real(kind=rp), allocatable :: kx(:), ky(:), kz(:)
     real(kind=rp), allocatable :: ax(:), ay(:), az(:)
     real(kind=rp), allocatable :: mode_phase(:)
     !> Device copies of the mode data (allocated only on device backends).
     type(c_ptr) :: kx_d = C_NULL_PTR
     type(c_ptr) :: ky_d = C_NULL_PTR
     type(c_ptr) :: kz_d = C_NULL_PTR
     type(c_ptr) :: ax_d = C_NULL_PTR
     type(c_ptr) :: ay_d = C_NULL_PTR
     type(c_ptr) :: az_d = C_NULL_PTR
     type(c_ptr) :: phase_d = C_NULL_PTR
     !> Device copies of the baseflow on the zone.
     type(c_ptr) :: u_bf_d = C_NULL_PTR
     type(c_ptr) :: v_bf_d = C_NULL_PTR
     type(c_ptr) :: w_bf_d = C_NULL_PTR
     !> Target velocity at the zone points (gathered, zone-local size).
     real(kind=rp), allocatable :: u_bf(:), v_bf(:), w_bf(:)
     !> Preview/validation flags.
     logical :: validate_only = .false.
     logical :: dump_flds = .false.
     character(len=NEKO_FNAME_LEN) :: dump_fname
     !> First-compute setup done (baseflow capture, dt diagnostics, dump).
     logical :: setup_done = .false.
     !> Latch for the gain*dt stability warning (checked every step to
     !! cover variable time steps; warn once).
     logical :: gain_dt_warned = .false.
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => fst_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => fst_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => fst_compute
  end type fst_source_term_t

contains

  !> Constructor from JSON.
  !! @param json The JSON object for the source term.
  !! @param fields The list of right-hand-side fields (f_x, f_y, f_z).
  !! @param coef The SEM coefficients.
  !! @param variable_name The name of the scheme owning this source term.
  subroutine fst_init_from_json(this, json, fields, coef, variable_name)
    class(fst_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

    character(len=:), allocatable :: zone_name, path, read_str, dump_name
    character(len=1), parameter :: dir_char(3) = ['x', 'y', 'z']
    logical, allocatable :: periodic_json(:)
    logical :: periodic(3), write_files
    real(kind=rp) :: start_time, end_time, ti, il, k_min, k_max
    real(kind=rp) :: lx_dom, ly_dom, lz_dom
    integer :: n_shells, modes_per_shell, seed
    integer :: d, n, n_zone_global, ierr
    character(len=LOG_SIZE) :: log_buf
    type(json_file) :: baseflow_subdict, spectrum_subdict

    call this%free()

    call neko_log%section("FST SOURCE TERM")

    if (fields%size() .ne. 3) then
       call neko_error("(FST) The fst source term is a momentum " // &
            "source and requires 3 right-hand-side fields." // &
            new_line('A') // "      Got the scheme '" // &
            trim(variable_name) // "' with a different field count.")
    end if
    if (coef%msh%gdim .ne. 3) then
       call neko_error("(FST) Only 3D meshes are supported")
    end if

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))
    call this%init_base(fields, coef, start_time, end_time)

    this%u => neko_registry%get_field_by_name("u")
    this%v => neko_registry%get_field_by_name("v")
    this%w => neko_registry%get_field_by_name("w")

    !
    ! Forcing region
    !
    call json_get(json, "zone_name", zone_name)
    this%zone => neko_point_zone_registry%get_point_zone(trim(zone_name))
    this%mask => this%zone%mask%get()

    n_zone_global = this%zone%size
    call MPI_Allreduce(MPI_IN_PLACE, n_zone_global, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    if (n_zone_global .eq. 0) then
       call neko_error("(FST) Point zone '" // trim(zone_name) // &
            "' contains no points")
    end if

    !
    ! Forcing parameters
    !
    call json_get(json, "gain", this%gain)
    if (this%gain .le. 0.0_rp) then
       call neko_error("(FST) gain must be > 0")
    end if
    if (json%valid_path("u_infty")) then
       call neko_error("(FST) 'u_infty' has been replaced by" // &
            new_line('A') // "      'convection_velocity': " // &
            "[Ucx, Ucy, Ucz]." // new_line('A') // &
            "      For the previous behavior use [u_infty, 0, 0].")
    end if
    block
      real(kind=rp), allocatable :: cv(:)
      call json_get(json, "convection_velocity", cv)
      if (size(cv) .ne. 3) then
         call neko_error("(FST) convection_velocity must have 3 elements")
      end if
      this%conv_vel = cv
      this%u_ref = norm2(this%conv_vel)
    end block
    call json_get_or_default(json, "ramp_time", this%ramp_time, 0.0_rp)
    if (this%ramp_time .lt. 0.0_rp) then
       call neko_error("(FST) ramp_time must be >= 0")
    end if

    !
    ! Fringe: per direction, either absent (flat, lambda = 1) or a smooth
    ! SIMSON fringe with all four parameters mandatory.
    !
    do d = 1, 3
       if (json%valid_path("fringe." // dir_char(d))) then
          this%fringe_smooth(d) = .true.
          call json_get(json, "fringe." // dir_char(d) // ".start", &
               this%fringe_start(d))
          call json_get(json, "fringe." // dir_char(d) // ".end", &
               this%fringe_end(d))
          call json_get(json, "fringe." // dir_char(d) // ".rise", &
               this%fringe_rise(d))
          call json_get(json, "fringe." // dir_char(d) // ".fall", &
               this%fringe_fall(d))

          if (this%fringe_end(d) .le. this%fringe_start(d)) then
             call neko_error("(FST) fringe." // dir_char(d) // &
                  ": end must be > start")
          end if
          if (this%fringe_rise(d) .le. 0.0_rp .or. &
               this%fringe_fall(d) .le. 0.0_rp) then
             call neko_error("(FST) fringe." // dir_char(d) // &
                  ": rise and fall must be > 0")
          end if
          if (this%fringe_end(d) - this%fringe_start(d) .lt. &
               this%fringe_rise(d) + this%fringe_fall(d)) then
             if (pe_rank .eq. 0) then
                call neko_warning("(FST) fringe." // dir_char(d) // &
                     ": rise + fall exceeds end - start." // &
                     new_line('A') // "      The ramps overlap and " // &
                     "lambda never reaches 1, so the effective" // &
                     new_line('A') // "      forcing is weaker than " // &
                     "gain suggests.")
             end if
          end if
       end if
    end do
    if (.not. any(this%fringe_smooth)) then
       if (pe_rank .eq. 0) then
          call neko_warning("(FST) No fringe direction given." // &
               new_line('A') // "      The forcing is applied at full " // &
               "strength on the entire zone.")
       end if
    end if

    !
    ! Periodicity for wavenumber quantization
    !
    periodic = .false.
    if (json%valid_path("periodic")) then
       call json_get(json, "periodic", periodic_json)
       if (size(periodic_json) .ne. 3) then
          call neko_error("(FST) periodic must have 3 elements")
       end if
       periodic = periodic_json
    end if

    !
    ! Spectrum
    !
    call json_get(json, "turbulence_intensity", ti)
    call json_get(json, "integral_length_scale", il)
    call json_get(json, "spectrum", spectrum_subdict)
    call json_get(spectrum_subdict, "n_shells", n_shells)
    call json_get(spectrum_subdict, "modes_per_shell", modes_per_shell)
    call json_get(spectrum_subdict, "k_min", k_min)
    call json_get(spectrum_subdict, "k_max", k_max)

    call json_get_or_default(json, "seed", seed, -143)
    call json_get_or_default(json, "write_files", write_files, .false.)
    call json_get_or_default(json, "files_output_path", path, "./fst_files")

    call json_get_or_default(json, "dump_fields", this%dump_flds, .false.)
    call json_get_or_default(json, "dump_file_name", dump_name, "fst_fields")
    this%dump_fname = trim(dump_name)
    call json_get_or_default(json, "validate_only", this%validate_only, &
         .false.)
    if (this%validate_only) this%dump_flds = .true.

    if (write_files .and. pe_rank .eq. 0) then
       call execute_command_line("mkdir -p " // trim(path))
    end if

    ! Full domain lengths for the periodic quantization
    n = coef%dof%size()
    lx_dom = glmax(coef%dof%x%x, n) - glmin(coef%dof%x%x, n)
    ly_dom = glmax(coef%dof%y%x, n) - glmin(coef%dof%y%x, n)
    lz_dom = glmax(coef%dof%z%x, n) - glmin(coef%dof%z%x, n)

    call this%spectrum%init(n_shells, modes_per_shell, k_min, k_max, &
         ti, il, this%u_ref, periodic, seed, write_files, path)
    call this%spectrum%generate(lx_dom, ly_dom, lz_dom)

    !
    ! Pack the mode data into the flat, kernel-friendly layout
    !
    call pack_modes(this)

    !
    ! Baseflow
    !
    allocate(this%u_bf(this%zone%size))
    allocate(this%v_bf(this%zone%size))
    allocate(this%w_bf(this%zone%size))
    this%u_bf = 0.0_rp
    this%v_bf = 0.0_rp
    this%w_bf = 0.0_rp

    call json_get(json, "baseflow", baseflow_subdict)
    call json_get(baseflow_subdict, "method", read_str)
    this%baseflow_method = trim(read_str)

    select case (this%baseflow_method)
    case ("initial_condition")
       ! Captured from the solution fields at the first compute call,
       ! since the initial condition is applied after source-term init.

    case ("constant")
       block
         real(kind=rp), allocatable :: values(:)
         call json_get(baseflow_subdict, "value", values)
         if (size(values) .ne. 3) then
            call neko_error("(FST) baseflow.value must have 3 elements")
         end if
         this%u_bf = values(1)
         this%v_bf = values(2)
         this%w_bf = values(3)
       end block

    case ("field")
       call baseflow_from_field(this, baseflow_subdict)

    case default
       call neko_error("(FST) '" // this%baseflow_method // &
            "' is not a valid baseflow method." // new_line('A') // &
            "      Use initial_condition, constant or field.")
    end select

    !
    ! Device: map the mode data and the baseflow. The baseflow for the
    ! "initial_condition" method is captured at the first compute call and
    ! copied there; the device buffers are allocated here in all cases.
    !
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%kx, this%kx_d, this%k_length)
       call device_map(this%ky, this%ky_d, this%k_length)
       call device_map(this%kz, this%kz_d, this%k_length)
       call device_map(this%ax, this%ax_d, this%k_length)
       call device_map(this%ay, this%ay_d, this%k_length)
       call device_map(this%az, this%az_d, this%k_length)
       call device_map(this%mode_phase, this%phase_d, this%k_length)

       call device_memcpy(this%kx, this%kx_d, this%k_length, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%ky, this%ky_d, this%k_length, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%kz, this%kz_d, this%k_length, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%ax, this%ax_d, this%k_length, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%ay, this%ay_d, this%k_length, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%az, this%az_d, this%k_length, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%mode_phase, this%phase_d, this%k_length, &
            HOST_TO_DEVICE, sync = .true.)

       if (this%zone%size .gt. 0) then
          call device_map(this%u_bf, this%u_bf_d, this%zone%size)
          call device_map(this%v_bf, this%v_bf_d, this%zone%size)
          call device_map(this%w_bf, this%w_bf_d, this%zone%size)
          call fst_baseflow_to_device(this)
       end if
    end if

    !
    ! Init-time diagnostics
    !
    call diagnostics_init(this, coef)

    write(log_buf, '(A,A)') "Baseflow method: ", this%baseflow_method
    call neko_log%message(log_buf)
    write(log_buf, '(A,E13.5)') "Gain: ", this%gain
    call neko_log%message(log_buf)
    write(log_buf, '(A,E13.5)') "Ramp time: ", this%ramp_time
    call neko_log%message(log_buf)

    call neko_log%end_section()

  end subroutine fst_init_from_json

  !> Import the baseflow from an fld file and gather it onto the zone.
  subroutine baseflow_from_field(this, baseflow_subdict)
    class(fst_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: baseflow_subdict

    character(len=:), allocatable :: read_str
    character(len=NEKO_FNAME_LEN) :: fname, mesh_fname
    logical :: interpolate
    type(json_file) :: interp_subdict
    type(field_t), target :: wku, wkv, wkw
    type(field_t), pointer :: pu, pv, pw

    call json_get(baseflow_subdict, "file_name", read_str)
    fname = trim(read_str)
    call json_get_or_default(baseflow_subdict, "interpolate", interpolate, &
         .false.)
    call json_get_or_default(baseflow_subdict, "mesh_file_name", read_str, &
         "none")
    mesh_fname = trim(read_str)
    call json_get_subdict_or_empty(baseflow_subdict, "interpolation", &
         interp_subdict)

    call wku%init(this%u%dof, "fst_bf_u")
    call wkv%init(this%u%dof, "fst_bf_v")
    call wkw%init(this%u%dof, "fst_bf_w")
    pu => wku
    pv => wkv
    pw => wkw

    call import_fields(trim(fname), interp_subdict, mesh_fname, &
         u = pu, v = pv, w = pw, interpolate = interpolate)

    call gather_masked(this%u_bf, wku%x, this%mask, wku%size(), &
         this%zone%size)
    call gather_masked(this%v_bf, wkv%x, this%mask, wkv%size(), &
         this%zone%size)
    call gather_masked(this%w_bf, wkw%x, this%mask, wkw%size(), &
         this%zone%size)

    call wku%free()
    call wkv%free()
    call wkw%free()

  end subroutine baseflow_from_field

  !> Packs the generated mode set into flat arrays used by both the CPU
  !! and the device kernels: the wavenumber components and the
  !! amplitude-weighted direction vectors a_j(m) = u_hat(m,j) *
  !! shell_amp(shell(m)). Folding the shell amplitude in here removes the
  !! shell indirection from the inner loop of every backend.
  subroutine pack_modes(this)
    class(fst_source_term_t), intent(inout) :: this

    integer :: m
    real(kind=rp) :: amp

    this%k_length = this%spectrum%k_length

    allocate(this%kx(this%k_length))
    allocate(this%ky(this%k_length))
    allocate(this%kz(this%k_length))
    allocate(this%ax(this%k_length))
    allocate(this%ay(this%k_length))
    allocate(this%az(this%k_length))
    allocate(this%mode_phase(this%k_length))

    do m = 1, this%k_length
       amp = this%spectrum%shell_amp(this%spectrum%shell(m))
       this%kx(m) = this%spectrum%k_num(m, 1)
       this%ky(m) = this%spectrum%k_num(m, 2)
       this%kz(m) = this%spectrum%k_num(m, 3)
       this%ax(m) = this%spectrum%u_hat(m, 1)*amp
       this%ay(m) = this%spectrum%u_hat(m, 2)*amp
       this%az(m) = this%spectrum%u_hat(m, 3)*amp
       this%mode_phase(m) = this%spectrum%phase(m)
    end do

  end subroutine pack_modes

  !> Copies the zone baseflow to the device.
  subroutine fst_baseflow_to_device(this)
    class(fst_source_term_t), intent(inout) :: this

    call device_memcpy(this%u_bf, this%u_bf_d, this%zone%size, &
         HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%v_bf, this%v_bf_d, this%zone%size, &
         HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%w_bf, this%w_bf_d, this%zone%size, &
         HOST_TO_DEVICE, sync = .true.)

  end subroutine fst_baseflow_to_device

  !> Destructor.
  subroutine fst_free(this)
    class(fst_source_term_t), intent(inout) :: this

    call this%spectrum%free()

    if (c_associated(this%kx_d)) call device_free(this%kx_d)
    if (c_associated(this%ky_d)) call device_free(this%ky_d)
    if (c_associated(this%kz_d)) call device_free(this%kz_d)
    if (c_associated(this%ax_d)) call device_free(this%ax_d)
    if (c_associated(this%ay_d)) call device_free(this%ay_d)
    if (c_associated(this%az_d)) call device_free(this%az_d)
    if (c_associated(this%phase_d)) call device_free(this%phase_d)
    if (c_associated(this%u_bf_d)) call device_free(this%u_bf_d)
    if (c_associated(this%v_bf_d)) call device_free(this%v_bf_d)
    if (c_associated(this%w_bf_d)) call device_free(this%w_bf_d)

    if (allocated(this%kx)) deallocate(this%kx)
    if (allocated(this%ky)) deallocate(this%ky)
    if (allocated(this%kz)) deallocate(this%kz)
    if (allocated(this%ax)) deallocate(this%ax)
    if (allocated(this%ay)) deallocate(this%ay)
    if (allocated(this%az)) deallocate(this%az)
    if (allocated(this%mode_phase)) deallocate(this%mode_phase)

    if (allocated(this%u_bf)) deallocate(this%u_bf)
    if (allocated(this%v_bf)) deallocate(this%v_bf)
    if (allocated(this%w_bf)) deallocate(this%w_bf)
    if (allocated(this%baseflow_method)) deallocate(this%baseflow_method)

    nullify(this%zone)
    nullify(this%mask)
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)

    this%fringe_smooth = .false.
    this%setup_done = .false.
    this%gain_dt_warned = .false.
    this%k_length = 0

    call this%free_base()

  end subroutine fst_free

  !> Computes the FST forcing and adds it to the right-hand-side fields.
  !! @param time The time state.
  subroutine fst_compute(this, time)
    class(fst_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: fu, fv, fw
    real(kind=rp) :: ramp, coeff
    character(len=LOG_SIZE) :: log_buf

    if (.not. this%setup_done) then

       if (this%baseflow_method .eq. "initial_condition") then
          ! On device backends the host copies of u, v, w may be stale
          ! after the initial condition was applied, so sync them before
          ! gathering, then push the captured baseflow back to the device.
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_memcpy(this%u%x, this%u%x_d, this%u%size(), &
                  DEVICE_TO_HOST, sync = .false.)
             call device_memcpy(this%v%x, this%v%x_d, this%v%size(), &
                  DEVICE_TO_HOST, sync = .false.)
             call device_memcpy(this%w%x, this%w%x_d, this%w%size(), &
                  DEVICE_TO_HOST, sync = .true.)
          end if

          call gather_masked(this%u_bf, this%u%x, this%mask, &
               this%u%size(), this%zone%size)
          call gather_masked(this%v_bf, this%v%x, this%mask, &
               this%v%size(), this%zone%size)
          call gather_masked(this%w_bf, this%w%x, this%mask, &
               this%w%size(), this%zone%size)

          if (NEKO_BCKND_DEVICE .eq. 1 .and. this%zone%size .gt. 0) then
             call fst_baseflow_to_device(this)
          end if
       end if

       ! Explicit-treatment strength diagnostic (heuristic): the fringe
       ! term is advanced explicitly, so gain*dt should stay O(1). The
       ! value is reported here; the threshold itself is checked every
       ! step below, since dt may vary during the run.
       write(log_buf, '(A,E13.5)') "[FST] gain * dt = ", this%gain*time%dt
       call neko_log%message(log_buf)

       if (this%dump_flds) then
          call dump_preview(this, time)
       end if

       if (this%validate_only) then
          call neko_error("(FST) validate_only requested: stopping" // &
               new_line('A') // "      after spectrum generation, " // &
               "diagnostics and field dump." // new_line('A') // &
               "      Review the log above and the dumped '" // &
               trim(this%dump_fname) // "' fld file," // &
               new_line('A') // "      then set validate_only = false " // &
               "to run.")
       end if

       this%setup_done = .true.
    end if

    ! Stability watch, every step (covers variable time steps): warn once
    ! the first time gain*dt exceeds the heuristic O(1) threshold.
    if (.not. this%gain_dt_warned .and. &
         this%gain*time%dt .gt. 1.0_rp) then
       this%gain_dt_warned = .true.
       if (pe_rank .eq. 0) then
          write(log_buf, '(A,E13.5,A,E13.5)') &
               "(FST) gain * dt exceeded 1 (gain * dt = ", &
               this%gain*time%dt, " at dt = ", time%dt
          call neko_warning(trim(log_buf) // ")." // new_line('A') // &
               "      The explicitly treated fringe forcing may be " // &
               "unstable at this" // new_line('A') // &
               "      time step (heuristic threshold 1). " // &
               "Reduce gain or dt.")
       end if
    end if

    ramp = time_ramp(time%t, this%start_time, this%ramp_time)
    if (ramp .le. 0.0_rp) return

    coeff = this%gain*ramp

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)

    if (this%zone%size .eq. 0) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call fst_source_term_compute_device(this%zone%size, &
            this%zone%mask%get_d(), &
            this%coef%dof%x%x_d, this%coef%dof%y%x_d, this%coef%dof%z%x_d, &
            this%u%x_d, this%v%x_d, this%w%x_d, fu%x_d, fv%x_d, fw%x_d, &
            this%u_bf_d, this%v_bf_d, this%w_bf_d, &
            this%k_length, this%kx_d, this%ky_d, this%kz_d, &
            this%ax_d, this%ay_d, this%az_d, this%phase_d, &
            this%conv_vel*time%t, coeff, &
            merge(1, 0, this%fringe_smooth), this%fringe_start, &
            this%fringe_end, this%fringe_rise, this%fringe_fall)
    else
       call fst_source_term_compute_cpu(this%u%dof%size(), this%zone%size, &
            this%mask, &
            this%coef%dof%x%x, this%coef%dof%y%x, this%coef%dof%z%x, &
            this%u%x, this%v%x, this%w%x, fu%x, fv%x, fw%x, &
            this%u_bf, this%v_bf, this%w_bf, &
            this%k_length, this%kx, this%ky, this%kz, &
            this%ax, this%ay, this%az, this%mode_phase, &
            this%conv_vel*time%t, coeff, &
            this%fringe_smooth, this%fringe_start, this%fringe_end, &
            this%fringe_rise, this%fringe_fall)
    end if

  end subroutine fst_compute






  !> Linear time ramp: 0 at start_time, 1 at start_time + ramp_time.
  pure function time_ramp(t, t_start, t_ramp) result(ramp)
    real(kind=rp), intent(in) :: t, t_start, t_ramp
    real(kind=rp) :: ramp

    if (t .le. t_start) then
       ramp = 0.0_rp
    else if (t_ramp .le. 0.0_rp) then
       ramp = 1.0_rp
    else
       ramp = min(1.0_rp, (t - t_start)/t_ramp)
    end if

  end function time_ramp

  !> Gather src at the masked points into dst (zone-local numbering).
  subroutine gather_masked(dst, src, mask, n, n_mask)
    integer, intent(in) :: n, n_mask
    real(kind=rp), intent(out) :: dst(n_mask)
    real(kind=rp), intent(in) :: src(n)
    integer, intent(in) :: mask(n_mask)

    integer :: idx

    do idx = 1, n_mask
       dst(idx) = src(mask(idx))
    end do

  end subroutine gather_masked

  !> Init-time diagnostics: grid resolution of the smallest FST wavelength
  !! in the zone, fringe support vs. zone bounding box, spectrum energy vs.
  !! target, and fringe strength. All thresholds are heuristic guidance and
  !! produce warnings, not errors.
  subroutine diagnostics_init(this, coef)
    class(fst_source_term_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    character(len=1), parameter :: dir_char(3) = ['x', 'y', 'z']
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: h_min, h_max, h_avg_max, h_avg, gap_sum
    real(kind=rp) :: lambda_min, ppw, bbox_min(3), bbox_max(3)
    integer :: gap_count
    real(kind=rp) :: tu_target, tu_rel_err, g_nondim, iso(3)
    integer :: idx, i, d, ierr

    call neko_log%message("--- FST validation diagnostics " // &
         "(thresholds are heuristic guidance) ---")

    !
    ! Grid resolution: the smallest FST wavelength 2*pi/k_max must be
    ! representable everywhere in the zone, so the binding measure is the
    ! COARSEST average GLL spacing of any element line in the zone
    ! (spectral criterion: ~pi points per wavelength of the average
    ! spacing; threshold 4 adds margin). Min/max adjacent gaps are
    ! reported for context only -- GLL gaps within one element vary by a
    ! factor ~3 at typical orders, so neither extreme is a valid basis.
    !
    call zone_gll_spacing(coef, this%mask, this%zone%size, h_min, h_max, &
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

    lambda_min = 2.0_rp*pi/this%spectrum%k_end
    ppw = lambda_min/h_avg_max

    write(log_buf, '(A,E13.5,A,E13.5)') "[FST] GLL gap in zone: min ", &
         h_min, " mean ", h_avg
    call neko_log%message(log_buf)
    write(log_buf, '(A,E13.5)') "[FST] GLL gap in zone: max ", h_max
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
            "zone spacing (spectral" // new_line('A') // &
            "      criterion ~pi with margin): the smallest scales may" // &
            new_line('A') // "      be under-resolved there. " // &
            "Reduce spectrum.k_max or refine.")
    end if

    !
    ! Zone bounding box (at the initial mesh) vs. fringe support. The
    ! fringe support [start, end] must be inside the zone, also at every
    ! later time if the mesh deforms.
    !
    bbox_min = huge(0.0_rp)
    bbox_max = -huge(0.0_rp)
    do idx = 1, this%zone%size
       i = this%mask(idx)
       bbox_min(1) = min(bbox_min(1), coef%dof%x%x(i, 1, 1, 1))
       bbox_max(1) = max(bbox_max(1), coef%dof%x%x(i, 1, 1, 1))
       bbox_min(2) = min(bbox_min(2), coef%dof%y%x(i, 1, 1, 1))
       bbox_max(2) = max(bbox_max(2), coef%dof%y%x(i, 1, 1, 1))
       bbox_min(3) = min(bbox_min(3), coef%dof%z%x(i, 1, 1, 1))
       bbox_max(3) = max(bbox_max(3), coef%dof%z%x(i, 1, 1, 1))
    end do
    call MPI_Allreduce(MPI_IN_PLACE, bbox_min, 3, MPI_REAL_PRECISION, &
         MPI_MIN, NEKO_COMM, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, bbox_max, 3, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)

    do d = 1, 3
       if (this%fringe_smooth(d)) then
          if (this%fringe_start(d) .lt. bbox_min(d) .or. &
               this%fringe_end(d) .gt. bbox_max(d)) then
             if (pe_rank .eq. 0) then
                call neko_warning("(FST) fringe." // dir_char(d) // &
                     " support [start, end] extends beyond" // &
                     new_line('A') // "      the zone bounding box at " // &
                     "the initial mesh: the forcing is" // &
                     new_line('A') // "      truncated where the zone " // &
                     "does not cover the fringe.")
             end if
          end if

          ! Fringe strength: gain * residence time of the flow crossing
          ! this fringe direction. The crossing speed is the convection
          ! velocity component along d. Heuristic: >= 5 for the forcing
          ! to imprint the target turbulence within one fringe passage.
          if (abs(this%conv_vel(d)) .gt. 0.0_rp) then
             g_nondim = this%gain &
                  * (this%fringe_end(d) - this%fringe_start(d)) &
                  / abs(this%conv_vel(d))
             write(log_buf, '(A,A,A,F10.3)') &
                  "[FST] gain * L_fringe/|U_c." , dir_char(d), "|: ", &
                  g_nondim
             call neko_log%message(log_buf)
             if (g_nondim .lt. 5.0_rp .and. pe_rank .eq. 0) then
                call neko_warning("(FST) gain * L_fringe/|U_c." // &
                     dir_char(d) // "| < 5." // new_line('A') // &
                     "      The fringe may be too weak to imprint the " // &
                     "target turbulence" // new_line('A') // &
                     "      within one fringe passage (heuristic " // &
                     "threshold 5)." // new_line('A') // &
                     "      Increase gain or widen the fringe.")
             end if
          else
             call neko_log%message("[FST] no U_c." // dir_char(d) // &
                  " component: residence check skipped for fringe." // &
                  dir_char(d))
          end if
       end if
    end do

    !
    ! Realized turbulence intensity of the generated mode set vs. target
    ! (diagnostics are valid on rank 0, where generation ran).
    !
    if (pe_rank .eq. 0) then
       tu_target = this%spectrum%ti*this%u_ref
       tu_rel_err = abs(this%spectrum%tu_uinf_estimate - tu_target) &
            / tu_target
       write(log_buf, '(A,E13.5,A,E13.5)') &
            "[FST] Tu*U target: ", tu_target, &
            " realized: ", this%spectrum%tu_uinf_estimate
       call neko_log%message(log_buf)
       ! Isotropy of the generated mode set: with M modes the component
       ! energies scatter by roughly 1/sqrt(M); a larger imbalance points
       ! to too few modes or to periodic quantization at low shells.
       iso = 3.0_rp*this%spectrum%energy/sum(this%spectrum%energy)
       write(log_buf, '(A,3F7.3)') &
            "[FST] component energy ratio E_u:E_v:E_w (x3/sum): ", iso
       call neko_log%message(log_buf)
       if (maxval(abs(iso - 1.0_rp)) .gt. 0.2_rp) then
          call neko_warning("(FST) Component energies deviate more " // &
               "than 20% from isotropy." // new_line('A') // &
               "      Increase spectrum.modes_per_shell / n_shells, or " // &
               "raise k_min relative" // new_line('A') // &
               "      to 2*pi/L in periodic directions.")
       end if

       if (tu_rel_err .gt. 0.1_rp) then
          call neko_warning("(FST) Realized Tu deviates more than 10%" // &
               new_line('A') // "      from the target: increase " // &
               "spectrum.n_shells and/or" // new_line('A') // &
               "      spectrum.modes_per_shell, or check that " // &
               "[k_min, k_max]" // new_line('A') // &
               "      captures the spectrum peak for the given " // &
               "integral_length_scale.")
       end if
    end if

    call neko_log%message("--- End FST validation diagnostics ---")

  end subroutine diagnostics_init

  !> GLL spacing statistics over every element that contains at least one
  !! zone point: minimum and maximum adjacent-node gap (context), and the
  !! largest per-line AVERAGE spacing h_avg_max (the basis for the
  !! spectral resolution criterion). Adjacent-node distances are Euclidean
  !! in physical space, so deformed/stretched/curved elements are handled
  !! correctly. The per-line average is (chord length of the GLL line) /
  !! (points - 1); the maximum over all lines and directions is the
  !! coarsest average spacing anywhere in the zone. Local (per-rank)
  !! results; reduce with MPI_MIN / MPI_MAX outside.
  subroutine zone_gll_spacing(coef, mask, n_mask, h_min, h_max, &
       h_avg_max, gap_sum, gap_count)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n_mask
    integer, intent(in) :: mask(n_mask)
    real(kind=rp), intent(out) :: h_min, h_max, h_avg_max
    !> Sum and count of all adjacent gaps (for the zone-wide mean
    !! spacing; reduce both with MPI_SUM and divide outside).
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

  end subroutine zone_gll_spacing

  !> Write the fringe field and a snapshot of the FST perturbation u'
  !! (without gain and ramp) at the current time and coordinates to an fld
  !! file for inspection, as fields 1-4: fringe, u', v', w'.
  subroutine dump_preview(this, time)
    class(fst_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(fld_file_output_t) :: fout
    type(field_t), pointer :: f_lam, f_up, f_vp, f_wp
    integer :: i1, i2, i3, i4

    call neko_log%message("[FST] Writing preview fields 1-4 " // &
         "(fringe, u', v', w') to '" // trim(this%dump_fname) // "'")

    call neko_scratch_registry%request_field(f_lam, i1, .true.)
    call neko_scratch_registry%request_field(f_up, i2, .true.)
    call neko_scratch_registry%request_field(f_vp, i3, .true.)
    call neko_scratch_registry%request_field(f_wp, i4, .true.)

    call fst_source_term_preview_cpu(this%u%dof%size(), this%zone%size, &
         this%mask, &
         this%coef%dof%x%x, this%coef%dof%y%x, this%coef%dof%z%x, &
         f_lam%x, f_up%x, f_vp%x, f_wp%x, &
         this%k_length, this%kx, this%ky, this%kz, &
         this%ax, this%ay, this%az, this%mode_phase, &
         this%conv_vel*time%t, &
         this%fringe_smooth, this%fringe_start, this%fringe_end, &
         this%fringe_rise, this%fringe_fall)

    ! The preview is built on the host; sync it so the sampler writes the
    ! same data on device backends.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(f_lam%x, f_lam%x_d, f_lam%size(), &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(f_up%x, f_up%x_d, f_up%size(), &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(f_vp%x, f_vp%x_d, f_vp%size(), &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(f_wp%x, f_wp%x_d, f_wp%size(), &
            HOST_TO_DEVICE, sync = .true.)
    end if

    call fout%init(sp, trim(this%dump_fname), 4)
    call fout%fields%assign_to_ptr(1, f_lam)
    call fout%fields%assign_to_ptr(2, f_up)
    call fout%fields%assign_to_ptr(3, f_vp)
    call fout%fields%assign_to_ptr(4, f_wp)
    call fout%sample(time%t)

    call neko_scratch_registry%relinquish_field(i1)
    call neko_scratch_registry%relinquish_field(i2)
    call neko_scratch_registry%relinquish_field(i3)
    call neko_scratch_registry%relinquish_field(i4)

  end subroutine dump_preview

end module fst_source_term