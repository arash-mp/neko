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
       json_get_or_lookup, json_get_subdict_or_empty
  use field, only : field_t
  use field_list, only : field_list_t
  use coefs, only : coef_t
  use source_term, only : source_term_t
  use fst_modes, only : fst_modes_t, fst_free_mapped
  use fst_fringe, only : fst_fringe_t, fst_ramp_t
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
  use device, only : device_map, device_memcpy, &
       HOST_TO_DEVICE, DEVICE_TO_HOST
  use fst_source_term_cpu, only : fst_source_term_compute_cpu
  use fst_source_term_device, only : fst_source_term_compute_device
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  use math, only : masked_gather_copy
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_MIN, MPI_MAX, &
       MPI_SUM, MPI_INTEGER
  implicit none
  private

  type, public, extends(source_term_t) :: fst_source_term_t
     !> The term does nothing unless enabled.
     logical :: enabled = .false.
     !> Mode set and fringe.
     type(fst_modes_t) :: modes
     type(fst_fringe_t) :: fringe
     !> Forcing region. Its mask follows the mesh (material points).
     class(point_zone_t), pointer :: zone => null()
     integer, pointer :: mask(:) => null()
     !> Velocity components.
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     !> Relaxation rate [1/time].
     real(kind=rp) :: gain = 0.0_rp
     !> Active window and ramp in time.
     type(fst_ramp_t) :: ramp
     !> One of initial_condition, constant or field.
     character(len=:), allocatable :: baseflow_method
     !> Device copies of the base flow.
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

    character(len=:), allocatable :: zone_name, read_str, dump_name
    integer :: n_zone_global, ierr
    character(len=LOG_SIZE) :: log_buf
    type(json_file) :: baseflow_subdict

    call this%free()

    call neko_log%section("FST SOURCE TERM")

    call this%ramp%init(json)
    call this%init_base(fields, coef, this%ramp%start, this%ramp%end)

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

    if (.not. json%valid_path("fringe") .and. pe_rank .eq. 0) then
       call neko_warning("(FST) No fringe given." // new_line('A') // &
            "      The forcing is applied at full strength on the " // &
            "entire zone.")
    end if
    call this%fringe%init(json, "fringe", .false.)

    call this%modes%init(json, coef)

    call json_get_or_default(json, "dump_fields", this%dump_flds, .false.)
    call json_get_or_default(json, "dump_file_name", dump_name, "fst_fields")
    this%dump_fname = trim(dump_name)
    call json_get_or_default(json, "validate_only", this%validate_only, &
         .false.)
    if (this%validate_only) this%dump_flds = .true.

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
    write(log_buf, '(A,E13.5)') "Ramp time: ", this%ramp%ramp
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
    call wku%copy_from(DEVICE_TO_HOST, .false.)
    call wkv%copy_from(DEVICE_TO_HOST, .false.)
    call wkw%copy_from(DEVICE_TO_HOST, .true.)

    call masked_gather_copy(this%u_bf, wku%x, this%mask, wku%size(), &
         this%zone%size)
    call masked_gather_copy(this%v_bf, wkv%x, this%mask, wkv%size(), &
         this%zone%size)
    call masked_gather_copy(this%w_bf, wkw%x, this%mask, wkw%size(), &
         this%zone%size)

    call wku%free()
    call wkv%free()
    call wkw%free()

  end subroutine baseflow_from_field

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

    call this%modes%free()
    call fst_free_mapped(this%u_bf, this%u_bf_d)
    call fst_free_mapped(this%v_bf, this%v_bf_d)
    call fst_free_mapped(this%w_bf, this%w_bf_d)
    if (allocated(this%baseflow_method)) deallocate(this%baseflow_method)

    nullify(this%zone)
    nullify(this%mask)
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)

    this%enabled = .false.
    this%setup_done = .false.
    this%gain_dt_warned = .false.

    call this%free_base()

  end subroutine fst_free

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
          call this%u%copy_from(DEVICE_TO_HOST, .false.)
          call this%v%copy_from(DEVICE_TO_HOST, .false.)
          call this%w%copy_from(DEVICE_TO_HOST, .true.)

          call masked_gather_copy(this%u_bf, this%u%x, this%mask, &
               this%u%size(), this%zone%size)
          call masked_gather_copy(this%v_bf, this%v%x, this%mask, &
               this%v%size(), this%zone%size)
          call masked_gather_copy(this%w_bf, this%w%x, this%mask, &
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

    ramp = this%ramp%value(real(time%t, kind=rp))
    if (ramp .le. 0.0_rp) return

    coeff = this%gain*ramp

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)

    call fst_apply(this, this%u, this%v, this%w, fu, fv, fw, &
         this%u_bf, this%v_bf, this%w_bf, &
         this%u_bf_d, this%v_bf_d, this%w_bf_d, &
         coeff, this%fringe%smooth, real(time%t, kind=rp))

  end subroutine fst_compute

  !> Run the kernel of the active backend at the zone points,
  !! f += coeff * lambda * (bf + u' - u).
  !! @param u, v, w Velocity.
  !! @param fu, fv, fw Fields the forcing is added to.
  !! @param u_bf, v_bf, w_bf Base flow at the zone points (host).
  !! @param u_bf_d, v_bf_d, w_bf_d The same on the device.
  !! @param coeff gain * ramp(t).
  !! @param smooth Directions with a smooth fringe, the others are flat.
  !! @param t Time.
  subroutine fst_apply(this, u, v, w, fu, fv, fw, u_bf, v_bf, w_bf, &
       u_bf_d, v_bf_d, w_bf_d, coeff, smooth, t)
    class(fst_source_term_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: fu, fv, fw
    real(kind=rp), intent(in) :: u_bf(*), v_bf(*), w_bf(*)
    type(c_ptr) :: u_bf_d, v_bf_d, w_bf_d
    real(kind=rp), intent(in) :: coeff, t
    logical, intent(in) :: smooth(3)

    if (this%zone%size .eq. 0) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call fst_source_term_compute_device(this%zone%size, &
            this%zone%mask%get_d(), &
            this%coef%dof%x%x_d, this%coef%dof%y%x_d, this%coef%dof%z%x_d, &
            u%x_d, v%x_d, w%x_d, fu%x_d, fv%x_d, fw%x_d, &
            u_bf_d, v_bf_d, w_bf_d, &
            this%modes%k_length, this%modes%kx_d, this%modes%ky_d, &
            this%modes%kz_d, this%modes%ax_d, this%modes%ay_d, &
            this%modes%az_d, this%modes%phase_d, &
            this%modes%conv_vel*t, coeff, &
            merge(1, 0, smooth), this%fringe%start, &
            this%fringe%end, this%fringe%rise, this%fringe%fall)
    else
       call fst_source_term_compute_cpu(u%dof%size(), this%zone%size, &
            this%mask, &
            this%coef%dof%x%x, this%coef%dof%y%x, this%coef%dof%z%x, &
            u%x, v%x, w%x, fu%x, fv%x, fw%x, &
            u_bf, v_bf, w_bf, &
            this%modes%k_length, this%modes%kx, this%modes%ky, this%modes%kz, &
            this%modes%ax, this%modes%ay, this%modes%az, this%modes%phase, &
            this%modes%conv_vel*t, coeff, &
            smooth, this%fringe%start, this%fringe%end, &
            this%fringe%rise, this%fringe%fall)
    end if

  end subroutine fst_apply

  !> Log resolution, fringe strength and spectrum checks, and warn when a
  !! value looks wrong. The thresholds are rules of thumb.
  subroutine diagnostics_init(this, coef)
    class(fst_source_term_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    character(len=1), parameter :: dir_char(3) = ['x', 'y', 'z']
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: bbox_min(3), bbox_max(3)
    real(kind=rp) :: g_nondim
    integer :: idx, i, d, ierr

    call neko_log%message("--- FST validation diagnostics " // &
         "(thresholds are heuristic guidance) ---")

    call this%modes%check_resolution(coef, this%mask, this%zone%size)

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
       if (this%fringe%smooth(d)) then
          if (this%fringe%start(d) .lt. bbox_min(d) .or. &
               this%fringe%end(d) .gt. bbox_max(d)) then
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
          if (abs(this%modes%conv_vel(d)) .gt. 0.0_rp) then
             g_nondim = this%gain &
                  * (this%fringe%end(d) - this%fringe%start(d)) &
                  / abs(this%modes%conv_vel(d))
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

    call neko_log%message("--- End FST validation diagnostics ---")

  end subroutine diagnostics_init

  !> Write fringe, u', v', w' to an fld file as fields 1-4.
  subroutine dump_preview(this, time)
    class(fst_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(fld_file_output_t) :: fout
    type(field_t), pointer :: f_lam, f_up, f_vp, f_wp, f_zero
    integer :: i1, i2, i3, i4, i5
    logical, parameter :: flat(3) = .false.

    call neko_log%message("[FST] Writing preview fields 1-4 " // &
         "(fringe, u', v', w') to '" // trim(this%dump_fname) // "'")

    call neko_scratch_registry%request_field(f_up, i2, .true.)
    call neko_scratch_registry%request_field(f_vp, i3, .true.)
    call neko_scratch_registry%request_field(f_wp, i4, .true.)
    call neko_scratch_registry%request_field(f_zero, i5, .true.)

    ! The fringe is cheap and built on the host.
    call neko_scratch_registry%request_field(f_lam, i1, .false.)
    f_lam%x = 0.0_rp
    call this%fringe%fill(f_lam%size(), this%zone%size, this%mask, &
         this%coef%dof%x%x, this%coef%dof%y%x, this%coef%dof%z%x, f_lam%x)
    call f_lam%copy_from(HOST_TO_DEVICE, .true.)

    ! u' with a flat fringe, gain 1 and
    ! zero base flow and velocity it adds 1 * 1 * (0 + u' - 0) = u'
    call fst_apply(this, f_zero, f_zero, f_zero, f_up, f_vp, f_wp, &
         f_zero%x, f_zero%x, f_zero%x, &
         f_zero%x_d, f_zero%x_d, f_zero%x_d, &
         1.0_rp, flat, real(time%t, kind=rp))

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
    call neko_scratch_registry%relinquish_field(i5)

  end subroutine dump_preview

end module fst_source_term
