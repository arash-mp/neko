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
!> Free-stream turbulence (FST) injected as a volume forcing.
!! Inside a point zone the velocity is relaxed towards a base flow plus a
!! synthetic turbulent field,
!!   f = gain * ramp(t) * lambda(x) * (u_bf + u' - u),
!! where u' is a sum of divergence-free Fourier modes from a von Karman
!! spectrum, convected with a constant velocity (frozen turbulence), and
!! lambda is a smooth fringe. u' and lambda are evaluated at the current
!! mesh coordinates every step, so the forcing is valid with ALE.
module fst_source_term
  use num_types, only : rp, sp
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, &
       json_get_or_lookup, json_get_or_lookup_or_default, &
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
  use device, only : device_map, device_unmap, device_memcpy, &
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

  type, public, extends(source_term_t) :: fst_source_term_t
     !> The term does nothing unless enabled.
     logical :: enabled = .false.
     !> Generated mode set.
     type(fst_spectrum_t) :: spectrum
     !> Forcing region. Its mask follows the mesh (material points).
     class(point_zone_t), pointer :: zone => null()
     integer, pointer :: mask(:) => null()
     !> Velocity components.
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     !> Relaxation rate [1/time].
     real(kind=rp) :: gain = 0.0_rp
     !> Convection velocity of the frozen turbulence and its magnitude.
     real(kind=rp) :: conv_vel(3) = 0.0_rp
     real(kind=rp) :: u_ref = 0.0_rp
     !> Length of the linear ramp after start_time.
     real(kind=rp) :: ramp_time = 0.0_rp
     !> Fringe per direction; directions that are not smooth are flat.
     logical :: fringe_smooth(3) = .false.
     real(kind=rp) :: fringe_start(3) = 0.0_rp
     real(kind=rp) :: fringe_end(3) = 0.0_rp
     real(kind=rp) :: fringe_rise(3) = 0.0_rp
     real(kind=rp) :: fringe_fall(3) = 0.0_rp
     !> One of initial_condition, constant or field.
     character(len=:), allocatable :: baseflow_method
     !> Modes in the layout used by the kernels, a_j = u_hat_j * amplitude.
     integer :: k_length = 0
     real(kind=rp), allocatable :: kx(:), ky(:), kz(:)
     real(kind=rp), allocatable :: ax(:), ay(:), az(:)
     real(kind=rp), allocatable :: mode_phase(:)
     !> Device copies of the modes and of the base flow.
     type(c_ptr) :: kx_d = C_NULL_PTR
     type(c_ptr) :: ky_d = C_NULL_PTR
     type(c_ptr) :: kz_d = C_NULL_PTR
     type(c_ptr) :: ax_d = C_NULL_PTR
     type(c_ptr) :: ay_d = C_NULL_PTR
     type(c_ptr) :: az_d = C_NULL_PTR
     type(c_ptr) :: phase_d = C_NULL_PTR
     type(c_ptr) :: u_bf_d = C_NULL_PTR
     type(c_ptr) :: v_bf_d = C_NULL_PTR
     type(c_ptr) :: w_bf_d = C_NULL_PTR
     !> Base flow at the zone points.
     real(kind=rp), allocatable :: u_bf(:), v_bf(:), w_bf(:)
     !> Preview: dump the forcing fields, then stop.
     logical :: validate_only = .false.
     logical :: dump_flds = .false.
     character(len=NEKO_FNAME_LEN) :: dump_fname
     !> First-step setup done, and gain*dt warning issued.
     logical :: setup_done = .false.
     logical :: gain_dt_warned = .false.
   contains
     procedure, pass(this) :: init => fst_init_from_json
     procedure, pass(this) :: free => fst_free
     procedure, pass(this) :: compute_ => fst_compute
  end type fst_source_term_t

contains

  !> Constructor from json.
  !! @param json The json object for the source term.
  !! @param fields The right-hand-side fields f_x, f_y, f_z.
  !! @param coef SEM coefficients.
  !! @param variable_name Name of the scheme the term belongs to.
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

    call json_get_or_lookup_or_default(json, "start_time", start_time, &
         0.0_rp)
    call json_get_or_lookup_or_default(json, "end_time", end_time, &
         huge(0.0_rp))
    call this%init_base(fields, coef, start_time, end_time)

    call json_get_or_default(json, "enabled", this%enabled, .false.)
    if (.not. this%enabled) then
       call neko_log%message("Disabled (set enabled = true to use it)")
       call neko_log%end_section()
       return
    end if

    if (fields%size() .ne. 3) then
       call neko_error("(FST) Needs the 3 momentum right-hand sides, " // &
            "got scheme '" // trim(variable_name) // "'")
    end if
    if (coef%msh%gdim .ne. 3) then
       call neko_error("(FST) Only 3D meshes are supported")
    end if

    this%u => neko_registry%get_field_by_name("u")
    this%v => neko_registry%get_field_by_name("v")
    this%w => neko_registry%get_field_by_name("w")

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

    call json_get_or_lookup(json, "gain", this%gain)
    if (this%gain .le. 0.0_rp) then
       call neko_error("(FST) gain must be > 0")
    end if
    block
      real(kind=rp), allocatable :: cv(:)
      call json_get_or_lookup(json, "convection_velocity", cv)
      if (size(cv) .ne. 3) then
         call neko_error("(FST) convection_velocity must have 3 elements")
      end if
      this%conv_vel = cv
      this%u_ref = norm2(this%conv_vel)
    end block
    call json_get_or_lookup_or_default(json, "ramp_time", this%ramp_time, &
         0.0_rp)
    if (this%ramp_time .lt. 0.0_rp) then
       call neko_error("(FST) ramp_time must be >= 0")
    end if

    do d = 1, 3
       if (json%valid_path("fringe." // dir_char(d))) then
          this%fringe_smooth(d) = .true.
          call json_get_or_lookup(json, "fringe." // dir_char(d) // &
               ".start", &
               this%fringe_start(d))
          call json_get_or_lookup(json, "fringe." // dir_char(d) // &
               ".end", &
               this%fringe_end(d))
          call json_get_or_lookup(json, "fringe." // dir_char(d) // &
               ".rise", &
               this%fringe_rise(d))
          call json_get_or_lookup(json, "fringe." // dir_char(d) // &
               ".fall", &
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

    periodic = .false.
    if (json%valid_path("periodic")) then
       call json_get(json, "periodic", periodic_json)
       if (size(periodic_json) .ne. 3) then
          call neko_error("(FST) periodic must have 3 elements")
       end if
       periodic = periodic_json
    end if

    call json_get_or_lookup(json, "turbulence_intensity", ti)
    call json_get_or_lookup(json, "integral_length_scale", il)
    call json_get(json, "spectrum", spectrum_subdict)
    call json_get_or_lookup(spectrum_subdict, "n_shells", n_shells)
    call json_get_or_lookup(spectrum_subdict, "modes_per_shell", &
         modes_per_shell)
    call json_get_or_lookup(spectrum_subdict, "k_min", k_min)
    call json_get_or_lookup(spectrum_subdict, "k_max", k_max)

    call json_get_or_lookup_or_default(json, "seed", seed, -143)
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

    n = coef%dof%size()
    lx_dom = glmax(coef%dof%x%x, n) - glmin(coef%dof%x%x, n)
    ly_dom = glmax(coef%dof%y%x, n) - glmin(coef%dof%y%x, n)
    lz_dom = glmax(coef%dof%z%x, n) - glmin(coef%dof%z%x, n)

    call this%spectrum%init(n_shells, modes_per_shell, k_min, k_max, &
         ti, il, this%u_ref, periodic, seed, write_files, path)
    call this%spectrum%generate(lx_dom, ly_dom, lz_dom)

    call pack_modes(this)

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
       ! Captured at the first step, the initial condition is set after init
    case ("constant")
       block
         real(kind=rp), allocatable :: values(:)
         call json_get_or_lookup(baseflow_subdict, "value", values)
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

    call diagnostics_init(this, coef)

    write(log_buf, '(A,A)') "Baseflow method: ", this%baseflow_method
    call neko_log%message(log_buf)
    write(log_buf, '(A,E13.5)') "Gain: ", this%gain
    call neko_log%message(log_buf)
    write(log_buf, '(A,E13.5)') "Ramp time: ", this%ramp_time
    call neko_log%message(log_buf)

    call neko_log%end_section()

  end subroutine fst_init_from_json

  !> Read the base flow from an fld file and keep its zone values.
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

    ! On device the imported values are only on the device
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call wku%copy_from(DEVICE_TO_HOST, .false.)
       call wkv%copy_from(DEVICE_TO_HOST, .false.)
       call wkw%copy_from(DEVICE_TO_HOST, .true.)
    end if

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

  !> Copy the modes into flat arrays and fold the shell amplitude into a_j.
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

  !> Copy the base flow to the device.
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

    call free_mapped(this%kx, this%kx_d)
    call free_mapped(this%ky, this%ky_d)
    call free_mapped(this%kz, this%kz_d)
    call free_mapped(this%ax, this%ax_d)
    call free_mapped(this%ay, this%ay_d)
    call free_mapped(this%az, this%az_d)
    call free_mapped(this%mode_phase, this%phase_d)
    call free_mapped(this%u_bf, this%u_bf_d)
    call free_mapped(this%v_bf, this%v_bf_d)
    call free_mapped(this%w_bf, this%w_bf_d)
    if (allocated(this%baseflow_method)) deallocate(this%baseflow_method)

    nullify(this%zone)
    nullify(this%mask)
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)

    this%enabled = .false.
    this%fringe_smooth = .false.
    this%setup_done = .false.
    this%gain_dt_warned = .false.
    this%k_length = 0

    call this%free_base()

  end subroutine fst_free

  !> Unmap (on device) and deallocate a host array created with device_map.
  subroutine free_mapped(x, x_d)
    real(kind=rp), allocatable, intent(inout) :: x(:)
    type(c_ptr), intent(inout) :: x_d

    if (allocated(x)) then
       if (c_associated(x_d)) call device_unmap(x, x_d)
       deallocate(x)
    end if

  end subroutine free_mapped

  !> Add the forcing to the right-hand side.
  !! @param time Current time state.
  subroutine fst_compute(this, time)
    class(fst_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: fu, fv, fw
    real(kind=rp) :: ramp, coeff
    character(len=LOG_SIZE) :: log_buf

    if (.not. this%enabled) return

    if (.not. this%setup_done) then

       if (this%baseflow_method .eq. "initial_condition") then
          ! The host copy of the velocity is stale on device backends
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

       write(log_buf, '(A,E13.5)') "[FST] gain * dt = ", this%gain*time%dt
       call neko_log%message(log_buf)

       if (this%dump_flds) then
          call dump_preview(this, time)
       end if

       if (this%validate_only) then
          call neko_error("(FST) validate_only is set, stopping here." // &
               new_line('A') // "      Check the log and the fields in '" // &
               trim(this%dump_fname) // "'." // new_line('A') // &
               "      Set validate_only to false to run.")
       end if

       this%setup_done = .true.
    end if

    ! Checked every step since dt can change
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

    ramp = time_ramp(real(time%t, kind=rp), this%start_time, &
         this%ramp_time)
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
            this%conv_vel*real(time%t, kind=rp), coeff, &
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
            this%conv_vel*real(time%t, kind=rp), coeff, &
            this%fringe_smooth, this%fringe_start, this%fringe_end, &
            this%fringe_rise, this%fringe_fall)
    end if

  end subroutine fst_compute

  !> Linear ramp from 0 at t_start to 1 at t_start + t_ramp.
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

  !> dst(i) = src(mask(i)).
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

  !> Log resolution, fringe strength and spectrum checks, and warn when a
  !! value looks wrong. The thresholds are rules of thumb.
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

    ! The smallest wavelength has to be resolved where the grid is coarsest.
    ! Spectral elements need about pi points per wavelength of the average
    ! spacing, so the check uses the coarsest line average (4 for margin).
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

    ! The fringe support should lie inside the zone, at t = 0 and, with a
    ! deforming mesh, at all later times.
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

          ! gain times the time the flow spends crossing the fringe
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

    ! Spectrum energies are only known on rank 0, where they were generated
    if (pe_rank .eq. 0) then
       tu_target = this%spectrum%ti*this%u_ref
       tu_rel_err = abs(this%spectrum%tu_uinf_estimate - tu_target) &
            / tu_target
       write(log_buf, '(A,E13.5,A,E13.5)') &
            "[FST] Tu*U target: ", tu_target, &
            " realized: ", this%spectrum%tu_uinf_estimate
       call neko_log%message(log_buf)
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

  !> GLL spacing over the elements touched by the zone: smallest and largest
  !! gap, the sum and count of all gaps, and the largest line average.
  !! Local values; reduce over ranks outside.
  subroutine zone_gll_spacing(coef, mask, n_mask, h_min, h_max, &
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

  end subroutine zone_gll_spacing

  !> Write fringe, u', v', w' to an fld file as fields 1-4.
  subroutine dump_preview(this, time)
    class(fst_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(fld_file_output_t) :: fout
    type(field_t), pointer :: f_lam, f_up, f_vp, f_wp
    integer :: i1, i2, i3, i4

    call neko_log%message("[FST] Writing preview fields 1-4 " // &
         "(fringe, u', v', w') to '" // trim(this%dump_fname) // "'")

    call neko_scratch_registry%request_field(f_lam, i1, .false.)
    call neko_scratch_registry%request_field(f_up, i2, .false.)
    call neko_scratch_registry%request_field(f_vp, i3, .false.)
    call neko_scratch_registry%request_field(f_wp, i4, .false.)

    f_lam%x = 0.0_rp
    f_up%x = 0.0_rp
    f_vp%x = 0.0_rp
    f_wp%x = 0.0_rp

    call fst_source_term_preview_cpu(this%u%dof%size(), this%zone%size, &
         this%mask, &
         this%coef%dof%x%x, this%coef%dof%y%x, this%coef%dof%z%x, &
         f_lam%x, f_up%x, f_vp%x, f_wp%x, &
         this%k_length, this%kx, this%ky, this%kz, &
         this%ax, this%ay, this%az, this%mode_phase, &
         this%conv_vel*real(time%t, kind=rp), &
         this%fringe_smooth, this%fringe_start, this%fringe_end, &
         this%fringe_rise, this%fringe_fall)

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
    call fout%free()

    call neko_scratch_registry%relinquish_field(i1)
    call neko_scratch_registry%relinquish_field(i2)
    call neko_scratch_registry%relinquish_field(i3)
    call neko_scratch_registry%relinquish_field(i4)

  end subroutine dump_preview

end module fst_source_term
