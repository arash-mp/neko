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
!> Writes the motion and the loads of the FSI bodies to CSV files, one file
!! per body.
module fsi_results
  use num_types, only : rp, dp
  use fsi_dynamics, only : fsi_body_t, fsi_body_acceleration
  use ale_manager, only : ale_manager_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get_or_lookup, &
       json_get_or_lookup_or_default
  use time_state, only : time_state_t
  use time_based_controller, only : time_based_controller_t
  use file, only : file_t
  use csv_file, only : csv_file_t
  use vector, only : vector_t
  use utils, only : neko_error
  implicit none
  private

  !> Names of the six DOFs in the column names.
  character(len=2), parameter, public :: FSI_DOF_NAMES(6) = &
       [character(len=2) :: 'x', 'y', 'z', 'rx', 'ry', 'rz']

  !> Longest header the CSV writer can hold.
  integer, parameter :: MAX_HEADER_LEN = 1024

  !> CSV output of the FSI bodies.
  type, public :: fsi_results_t
     !> Whether the results are also printed to the log.
     logical :: log_results = .true.
     !> Decides at which steps the files are written.
     type(time_based_controller_t) :: controller
     !> One file and one row per body.
     type(file_t), allocatable :: files(:)
     type(vector_t), allocatable :: rows(:)
   contains
     procedure, pass(this) :: init => fsi_results_init
     procedure, pass(this) :: restart => fsi_results_restart
     procedure, pass(this) :: write => fsi_results_write
     procedure, pass(this) :: free => fsi_results_free
  end type fsi_results_t

contains

  !> Constructor.
  !! @param params The case file.
  !! @param bodies The FSI bodies.
  !! @param fsi_dof_map Row of each DOF in the structural system, 0 if the
  !! DOF is not active; (nbodies, 6).
  !! @param extra_names Names of the scheme-specific columns of each body;
  !! (max(n_extra), nbodies).
  !! @param n_extra Number of scheme-specific columns of each body.
  subroutine fsi_results_init(this, params, bodies, fsi_dof_map, extra_names, &
       n_extra)
    class(fsi_results_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(fsi_body_t), intent(in) :: bodies(:)
    integer, intent(in) :: fsi_dof_map(:,:)
    character(len=*), intent(in) :: extra_names(:,:)
    integer, intent(in) :: n_extra(:)
    character(len=:), allocatable :: control, dir, fname, header
    character(len=1), parameter :: xyz(3) = ['x', 'y', 'z']
    real(kind=dp) :: value, start_time, end_time
    integer :: i, k, ncol, n

    call this%free()

    call json_get_or_default(params, 'case.fluid.fsi.log_results', &
         this%log_results, .true.)
    call json_get_or_default(params, 'case.fluid.fsi.output_control', &
         control, 'tsteps')
    call json_get_or_default(params, 'case.fluid.fsi.output_value', value, &
         1.0_dp)
    call json_get_or_lookup_or_default(params, 'case.time.start_time', &
         start_time, 0.0_dp)
    call json_get_or_lookup(params, 'case.time.end_time', end_time)
    call this%controller%init(start_time, end_time, control, value, &
         write_at_start = .false.)

    call json_get_or_default(params, 'case.output_directory', dir, '')
    if (len_trim(dir) .gt. 0) then
       if (dir(len_trim(dir):len_trim(dir)) .ne. '/') dir = trim(dir) // '/'
    end if

    allocate(this%files(size(bodies)))
    allocate(this%rows(size(bodies)))

    do i = 1, size(bodies)
       ! The time is written by the CSV writer; the row holds the rest
       header = 'time'
       ncol = 0
       call add_active('disp_')
       do k = 1, 3
          call add('total_disp_' // xyz(k))
       end do
       do k = 1, 3
          call add('total_rot_' // xyz(k))
       end do
       call add_active('vel_')
       call add_active('acc_')
       do k = 1, 3
          call add('force_' // xyz(k))
       end do
       do k = 1, 3
          call add('torque_' // xyz(k))
       end do
       do k = 1, 3
          call add('visc_force_' // xyz(k))
       end do
       do k = 1, 3
          call add('visc_torque_' // xyz(k))
       end do
       do k = 1, n_extra(i)
          call add(trim(extra_names(k, i)))
       end do
       if (len(header) .gt. MAX_HEADER_LEN) then
          call neko_error("FSI results: the header of body " // &
               trim(bodies(i)%name) // " is too long for the CSV writer")
       end if

       ! The name as given, with the csv suffix added if it is missing
       fname = trim(bodies(i)%output_filename)
       n = len(fname)
       if (n .lt. 4) then
          fname = fname // '.csv'
       else if (fname(n-3:n) .ne. '.csv') then
          fname = fname // '.csv'
       end if
       do k = 1, i - 1
          if (trim(this%files(k)%get_fname()) .eq. dir // fname) then
             call neko_error("FSI results: bodies " // &
                  trim(bodies(k)%name) // " and " // trim(bodies(i)%name) // &
                  " write to the same file " // fname)
          end if
       end do

       call this%files(i)%init(dir // fname)
       call this%files(i)%set_header(header)
       call this%files(i)%set_overwrite(.true.)
       call this%rows(i)%init(ncol)
    end do

  contains

    !> Add one column to the header.
    subroutine add(name)
      character(len=*), intent(in) :: name

      header = header // ',' // name
      ncol = ncol + 1
    end subroutine add

    !> Add one column per active DOF of body `i`.
    subroutine add_active(prefix)
      character(len=*), intent(in) :: prefix
      integer :: j

      do j = 1, 6
         if (fsi_dof_map(i, j) .gt. 0) then
            call add(prefix // trim(FSI_DOF_NAMES(j)))
         end if
      end do
    end subroutine add_active

  end subroutine fsi_results_init

  !> Continue the files of the run that wrote the checkpoint, instead of
  !! starting new ones.
  !! @param time The time of the restart.
  subroutine fsi_results_restart(this, time)
    class(fsi_results_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    logical :: exists
    integer :: i

    if (.not. allocated(this%files)) return

    call this%controller%set_counter(time)

    do i = 1, size(this%files)
       call this%files(i)%set_overwrite(.false.)
       inquire(file = trim(this%files(i)%get_fname()), exist = exists)
       if (exists) then
          select type (f => this%files(i)%file_type)
          type is (csv_file_t)
             f%header_is_written = .true.
          end select
       end if
    end do
  end subroutine fsi_results_restart

  !> Write one row per body, if a row is due at this step.
  !! @param time The current time.
  !! @param bodies The FSI bodies.
  !! @param fsi_dof_map Row of each DOF in the structural system, 0 if the
  !! DOF is not active; (nbodies, 6).
  !! @param ale The ALE manager, for the pivot and the orientation.
  !! @param force Fluid force and torque on each body; (6, nbodies).
  !! @param viscous Their viscous part; (6, nbodies).
  !! @param extra Values of the scheme-specific columns of each body.
  subroutine fsi_results_write(this, time, bodies, fsi_dof_map, ale, force, &
       viscous, extra)
    class(fsi_results_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(fsi_body_t), intent(in) :: bodies(:)
    integer, intent(in) :: fsi_dof_map(:,:)
    type(ale_manager_t), intent(in) :: ale
    real(kind=rp), intent(in) :: force(:,:), viscous(:,:)
    real(kind=rp), intent(in) :: extra(:,:)
    real(kind=rp) :: acc(6), R(3,3)
    integer :: i, k, n, a, n_extra

    if (.not. allocated(this%files)) return
    if (.not. this%controller%check(time)) return

    do i = 1, size(bodies)
       associate (x => this%rows(i)%x)
         n = 0
         do k = 1, 6
            if (fsi_dof_map(i, k) .gt. 0) then
               n = n + 1
               x(n) = bodies(i)%disp_rel(k)
            end if
         end do

         ! Total motion: the pivot displacement, and roll, pitch and yaw of
         ! the body
         a = bodies(i)%ale_id
         x(n+1:n+3) = ale%ale_pivot(a)%pos - ale%config%bodies(a)%rot_center
         R = ale%body_rot_matrices(:, :, a)
         x(n+4) = atan2(R(3,2), R(3,3))
         x(n+5) = atan2(-R(3,1), sqrt(R(3,2)**2 + R(3,3)**2))
         x(n+6) = atan2(R(2,1), R(1,1))
         n = n + 6

         do k = 1, 6
            if (fsi_dof_map(i, k) .gt. 0) then
               n = n + 1
               x(n) = bodies(i)%body_vel(k)
            end if
         end do

         acc = fsi_body_acceleration(bodies(i))
         do k = 1, 6
            if (fsi_dof_map(i, k) .gt. 0) then
               n = n + 1
               x(n) = acc(k)
            end if
         end do

         x(n+1:n+6) = force(:, i)
         x(n+7:n+12) = viscous(:, i)
         n = n + 12

         n_extra = size(x) - n
         x(n+1:n+n_extra) = extra(1:n_extra, i)
       end associate

       call this%files(i)%write(this%rows(i), time%t)
    end do

    call this%controller%register_execution(time)
  end subroutine fsi_results_write

  !> Destructor.
  subroutine fsi_results_free(this)
    class(fsi_results_t), intent(inout) :: this
    integer :: i

    if (allocated(this%rows)) then
       do i = 1, size(this%rows)
          call this%rows(i)%free()
       end do
       deallocate(this%rows)
    end if
    if (allocated(this%files)) then
       do i = 1, size(this%files)
          call this%files(i)%free()
       end do
       deallocate(this%files)
    end if
    call this%controller%free()
  end subroutine fsi_results_free

end module fsi_results
