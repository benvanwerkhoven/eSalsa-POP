!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module tavg

!BOP
! !MODULE: tavg
! !DESCRIPTION:
!  This module contains data types and routines for computing running
!  time-averages of selected fields and writing this data to files.
!
! !REVISION HISTORY:
!  CVS:$Id: tavg.F90,v 1.31 2003/12/23 22:32:16 pwjones Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use blocks
   use distribution
   use domain
   use constants
   use prognostic
   use grid
   use time_management
   use global_reductions
   use gather_scatter
   use broadcast
   use io
   use exit_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: read_tavg_namelist,    &
             init_tavg,             &
             define_tavg_field,     &
             accumulate_tavg_field, &
             tavg_requested,        &
             write_tavg,            &
             read_tavg,             &
             tavg_set_flag

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public :: &
      ltavg_on      = .false., & ! tavg file output wanted
      ltavg_restart = .false., & ! run started from restart
      tavg_do_MOC_Strength = .false., &
      tavg_no_write = .false., &
      skip_tavg_dump = .false.

   integer (int_kind), parameter, public :: &
      tavg_method_unknown = 0,              &
      tavg_method_avg     = 1,              &
      tavg_method_min     = 2,              &
      tavg_method_max     = 3

   integer (i4), public ::     &
      tavg_freq,       &! frequency of tavg output
      tavg_start        ! start tavg after tavg_start

   character (char_len), public ::    &
      tavg_infile,            & ! filename for restart input
      tavg_outfile,           & ! root filename for tavg output
      tavg_fmt_in,            & ! format (nc or bin) for reading
      tavg_fmt_out              ! format (nc or bin) for writing

   character (char_len), public :: &
      tavg_freq_opt,       &! choice for frequency of tavg output
      tavg_start_opt,      &! choice for starting averaging
      tavg_contents,       &! filename for choosing fields for output
      char_temp             ! temporary for manipulating fields

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  tavg field descriptor data type and array of such types
!
!-----------------------------------------------------------------------

   type :: tavg_field_desc
      character(char_len)     :: short_name     ! short name for field
      character(char_len)     :: long_name      ! long descriptive name
      character(char_len)     :: units          ! units
      character(4)            :: grid_loc       ! location in grid
      real (r4)               :: missing_value  ! value on land pts
      real (r4), dimension(2) :: valid_range    ! min/max
      integer (i4)            :: ndims          ! num dims (2 or 3)
      integer (i4)            :: buf_loc        ! location in buffer
      integer (i4)            :: method         ! method for averaging
      integer (i4)            :: field_loc      ! grid location and field
      integer (i4)            :: field_type     !  type for io, ghost cells
   end type

   integer (int_kind), parameter :: &
      max_avail_tavg_fields = 200    ! limit on available fields - can
                                     !   be pushed as high as necessary

   integer (int_kind) ::           &
      num_avail_tavg_fields = 0,   &! current number of defined fields
      num_requested_tavg_fields,   &! number of fields requested
      tavg_flag                     ! time flag for writing tavg files

   type (tavg_field_desc), dimension(max_avail_tavg_fields) :: &
      avail_tavg_fields

!-----------------------------------------------------------------------
!
!  AMOC diagnostic
!
!-----------------------------------------------------------------------

   real (r8), dimension(:),public,allocatable ::  &
      lat_aux_center,     &! cell center latitude values (degrees north)
      lat_aux_edge         ! cell edge   latitude values (degrees north)

   integer (int_kind), dimension(:,:), allocatable ::  &
      ATLANTIC_MASK_LAT_AUX   ! Global atlantic mask (master_task only)

   integer (int_kind), dimension(nx_block,ny_block,max_blocks_clinic), public :: &
      ATLANTIC_MASK_LAT       ! Local (blocks) atlantic mask

   integer (int_kind) :: &
      lat_aux_region_start    ! starting latitude indices for AMOC

   real (r4), dimension(:), allocatable :: &  ! Array for read-in atlantic mask
      ATLANTIC_MASK_LONG, &
      ATLANTIC_MASK_LONG_MIN, &
      ATLANTIC_MASK_LONG_MAX

   character (char_len), public ::    &       ! Read in atlantic mask file
      atlantic_mask_file

   real (r8), dimension(:,:), allocatable ::  &
      TLATD_G              ! latitude of T points in degrees (global array)

   real (r8), dimension(:,:,:), public, allocatable ::  &
      TAVG_MOC_G             ! meridional overturning circulation (master_task only)

   integer (int_kind), public :: &
      n_lat_aux_grid          ! Number of lat grid points in aux grid

   real (r8), public ::          &
      lat_aux_begin,        & ! beginning latitude for the auxilary
                              !    grid (degrees north)
      lat_aux_end             ! ending latitude for the auxilary
                              !    grid (degrees north)

   integer (int_kind) :: &
      amoc_strength_target_lat_idx, &   ! Latitude of the target AMOC scalar probe
      amoc_strength_target_depth_idx    ! Depth of the target AMOC scalar probe

   real (r8), public :: &               ! Scalar estimate of the AMOC strength
      amoc_strength

   integer (int_kind) ::  & ! indices needed for tavg AMOC diagnostics
      tavg_id_WVEL,       &
      tavg_id_VVEL,       &
      tavg_loc_WVEL,      &
      tavg_loc_VVEL

!-----------------------------------------------------------------------
!
!  buffers for holding running tavg variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_bufsize_2d,   &    ! size of buffer for 2d fields
      tavg_bufsize_3d         ! size of buffer for 3d fields

   real (r4), dimension(:,:,:,:), allocatable :: &
      TAVG_BUF_2D         ! buffer for holding accumulated sums

   real (r4), dimension(:,:,:,:,:), allocatable :: &
      TAVG_BUF_3D         ! buffer for holding accumulated sums

   integer (int_kind), dimension(:), allocatable :: &
      TAVG_BUF_2D_METHOD,  &! method for each requested 2d field
      TAVG_BUF_3D_METHOD    ! method for each requested 3d field

!-----------------------------------------------------------------------
!
!  variables for writing data
!
!-----------------------------------------------------------------------

   integer (i4) ::     &
      tavg_freq_iopt,  &! frequency option for writing tavg
      tavg_start_iopt   ! start after option

   type (datafile) :: tavg_file_desc    ! IO file descriptor

   type (io_field_desc), target :: &
      TAVG_iodesc                  ! io descriptor for tavg fields

!-----------------------------------------------------------------------
!
!  scalars
!
!-----------------------------------------------------------------------

   real (r8) ::  &
      tavg_sum,  &     ! accumulated time (in seconds)
      dtavg            ! current time step

   character (10) :: &
      beg_date       ! date on which the current accumulated sum
                     ! was started (not the tavg_start date)
!-----------------------------------------------------------------------
!
!  namelist options
!
!-----------------------------------------------------------------------

   character (33), parameter :: &
      freq_fmt = "('tavg diagnostics every ',i6,a8)"

   character (44), parameter :: &
      start_fmt = "('tavg sums accumulated starting at ',a5,i8)"

   character (char_len) :: exit_string

!EOC
!***********************************************************************

 contains

!***********************************************************************
!EOP
! !IROUTINE: init_tavg, split into read_tavg_namelist
! !INTERFACE:

 subroutine read_tavg_namelist

! !DESCRIPTION:
!  This routine initializes tavg options and reads in contents file to
!  determine which fields for which the user wants time-averaged data.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::      &
      nml_error          ! namelist i/o error flag

   namelist /tavg_nml/ tavg_freq_opt, tavg_freq, tavg_infile,       &
                       tavg_outfile, tavg_contents, tavg_start_opt, &
                       tavg_start, tavg_fmt_in, tavg_fmt_out,       &
                       skip_tavg_dump,                              &
                       tavg_do_MOC_Strength, n_lat_aux_grid, lat_aux_begin, &
                       lat_aux_end, atlantic_mask_file

!-----------------------------------------------------------------------
!
!  read tavg file output frequency and filenames from namelist
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a12)') 'Tavg options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

   tavg_freq_iopt = freq_opt_never
   tavg_freq      = 100000
   tavg_start_iopt = start_opt_nstep
   tavg_start      = 0
   tavg_infile    = 'unknown_tavg_infile'
   tavg_outfile   = 't'
   tavg_contents  = 'unknown_tavg_contents'
   lat_aux_begin      = -90.0_r8
   lat_aux_end        =  90.0_r8
   n_lat_aux_grid     =  180
   atlantic_mask_file = 'unknown_atlantic_mask_file'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=tavg_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading tavg_nml')
   endif

end subroutine read_tavg_namelist
!split here

subroutine init_tavg
! !DESCRIPTION:
!  This routine initializes tavg options and reads in contents file to
!  determine which fields for which the user wants time-averaged data.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::         &
      n,                   &! dummy index
      iblock,              &! local block index
      loc,                 &! location of field in buffer
      nu,                  &! unit for contents input file
      cindex,              &! character index for manipulating strings
      contents_error        ! error flag for contents file read

   if (my_task == master_task) then
      select case (tavg_freq_opt)
      case ('never')
         tavg_freq_iopt = freq_opt_never
         write(stdout,'(a20)') 'tavg diagnostics off'
      case ('nyear')
         tavg_freq_iopt = freq_opt_nyear
         write(stdout,freq_fmt) tavg_freq,' years  '
      case ('nmonth')
         tavg_freq_iopt = freq_opt_nmonth
         write(stdout,freq_fmt) tavg_freq,' months '
      case ('nday')
         tavg_freq_iopt = freq_opt_nday
         write(stdout,freq_fmt) tavg_freq,' days   '
      case ('nhour')
         tavg_freq_iopt = freq_opt_nhour
         write(stdout,freq_fmt) tavg_freq,' hours  '
      case ('nsecond')
         tavg_freq_iopt = freq_opt_nsecond
         write(stdout,freq_fmt) tavg_freq,' seconds'
      case ('nstep')
         tavg_freq_iopt = freq_opt_nstep
         write(stdout,freq_fmt) tavg_freq,' steps  '
      case default
         tavg_freq_iopt = -1000
      end select

      if (tavg_freq_iopt /= freq_opt_never) then
         select case (tavg_start_opt)
         case ('nstep')
            tavg_start_iopt = start_opt_nstep
            write(stdout,start_fmt) 'step ', tavg_start
         case ('nday')
            tavg_start_iopt = start_opt_nday
            write(stdout,start_fmt) 'day  ', tavg_start
         case ('nyear')
            tavg_start_iopt = start_opt_nyear
            write(stdout,start_fmt) 'year ', tavg_start
         case ('date')
            tavg_start_iopt = start_opt_date
            write(stdout,start_fmt) '     ', tavg_start
         case default
            tavg_start_iopt = -1000
         end select
      endif

   endif

   call broadcast_scalar(tavg_freq_iopt, master_task)

   if (tavg_freq_iopt == -1000) then
      call exit_POP(sigAbort,'unknown option for tavg file frequency')
   else if (tavg_freq_iopt /= freq_opt_never) then
      call broadcast_scalar(tavg_freq,            master_task)
      call broadcast_scalar(tavg_start_iopt,      master_task)
      call broadcast_scalar(tavg_start,           master_task)
      call broadcast_scalar(tavg_infile,          master_task)
      call broadcast_scalar(tavg_outfile,         master_task)
      call broadcast_scalar(tavg_contents,        master_task)
      call broadcast_scalar(tavg_fmt_in,          master_task)
      call broadcast_scalar(tavg_fmt_out,         master_task)
      call broadcast_scalar(skip_tavg_dump,       master_task)
      call broadcast_scalar(tavg_do_MOC_Strength, master_task)
      call broadcast_scalar(n_lat_aux_grid,       master_task)
      call broadcast_scalar(lat_aux_begin,        master_task)
      call broadcast_scalar(lat_aux_end,          master_task)

      if (tavg_start_iopt == -1000) then
         call exit_POP(sigAbort,'unknown option for tavg start option')
      endif

   endif

!-----------------------------------------------------------------------
!
!  initialize time flag for writing tavg files
!
!-----------------------------------------------------------------------

   call init_time_flag('tavg', tavg_flag,              &
                        owner        = 'init_tavg',    &
                        default      =.false.,         &
                        freq_opt     = tavg_freq_iopt, &
                        freq         = tavg_freq       )

!-----------------------------------------------------------------------
!
!  read contents file to determine which fields to dump
!
!-----------------------------------------------------------------------

   if (tavg_freq_iopt /= freq_opt_never) then

      tavg_bufsize_2d = 0
      tavg_bufsize_3d = 0

      call get_unit(nu)

      if (my_task == master_task) then
         open(nu, file=tavg_contents, status='old')
         read(nu,*) num_requested_tavg_fields
         write(stdout,'(a38)') 'tavg diagnostics requested for fields:'
      endif

      call broadcast_scalar(num_requested_tavg_fields, master_task)

      contents_error = 0

      do n=1,num_requested_tavg_fields
         if (my_task == master_task) then
            read(nu,'(a80)',iostat=contents_error) char_temp
            char_temp = adjustl(char_temp)
            cindex = index(char_temp,' ')
            char_temp(cindex:) = ' '
            write(stdout,*) '  ',trim(char_temp)
         endif

         call broadcast_scalar(contents_error, master_task)
         if (contents_error /= 0) then
            call exit_POP(sigAbort,'error reading tavg contents')
         endif

         call broadcast_scalar(char_temp, master_task)
         call request_tavg_field(trim(char_temp))
      end do

      call release_unit(nu)

      !*** allocate and initialize running tavg buffers

      allocate(                                                            &
         TAVG_BUF_2D(nx_block,ny_block,   nblocks_clinic,tavg_bufsize_2d), &
         TAVG_BUF_3D(nx_block,ny_block,km,nblocks_clinic,tavg_bufsize_3d), &
         TAVG_BUF_2D_METHOD(tavg_bufsize_2d),                              &
         TAVG_BUF_3D_METHOD(tavg_bufsize_3d))

      tavg_sum = c0
      call time_stamp('now','ymd',date_string=beg_date)
      if (tavg_freq_iopt == freq_opt_nstep) &
         write(beg_date,'(i10)') nsteps_total

      !*** initialize buffers based on requested method

      !$OMP PARALLEL DO PRIVATE(n,loc)

      do iblock=1,nblocks_clinic
      do n = 1,num_avail_tavg_fields  ! check all available fields

         loc = abs(avail_tavg_fields(n)%buf_loc)
         if (loc /= 0) then  ! field is actually requested and in buffer

            if (avail_tavg_fields(n)%ndims == 2) then

               TAVG_BUF_2D_METHOD(loc) = avail_tavg_fields(n)%method

               select case (TAVG_BUF_2D_METHOD(loc))
               case (tavg_method_avg)
                  TAVG_BUF_2D(:,:,  iblock,loc) = c0
               case (tavg_method_min)
                  TAVG_BUF_2D(:,:,  iblock,loc) = bignum
               case (tavg_method_max)
                  TAVG_BUF_2D(:,:,  iblock,loc) = -bignum
               case default
                  TAVG_BUF_2D(:,:,  iblock,loc) = c0
               end select

            else if (avail_tavg_fields(n)%ndims == 3) then

               TAVG_BUF_3D_METHOD(loc) = avail_tavg_fields(n)%method
               select case (TAVG_BUF_3D_METHOD(loc))
               case (tavg_method_avg)
                  TAVG_BUF_3D(:,:,:,iblock,loc) = c0
               case (tavg_method_min)
                  TAVG_BUF_3D(:,:,:,iblock,loc) = bignum
               case (tavg_method_max)
                  TAVG_BUF_3D(:,:,:,iblock,loc) = -bignum
               case default
                  TAVG_BUF_3D(:,:,:,iblock,loc) = c0
               end select

            endif

         endif
      end do
      end do
      !$OMP END PARALLEL DO

      if (tavg_do_MOC_Strength) then
        call init_tavg_amoc
      endif

   endif

!-----------------------------------------------------------------------
!
!  read restart file if necessary
!
!-----------------------------------------------------------------------

   !*** make sure tavg flag is set correctly
   call tavg_set_flag(flagonly=.true.)
   call eval_time_flag(tavg_flag) ! evaluates time_flag(tavg_flag)%value via time_to_do

   if (ltavg_on .and. (ltavg_restart .and. (.not. skip_tavg_dump))) then
      !*** do not read restart if last restart was at a tavg dump
      !*** interval (should start new tavg sums in this case)
      if (.not. check_time_flag(tavg_flag)) call read_tavg
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine init_tavg

!***********************************************************************
!BOP
! !IROUTINE: tavg_set_flag
! !INTERFACE:

 subroutine tavg_set_flag(flagonly)

! !DESCRIPTION:
!  This routine checks the time avg option and tavg start condition
!  to see whether tavg sums should be accumulated.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in), optional :: &
      flagonly        ! if true, only sets ltavg_on without advancing
                      !  time interval

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: n ! loop index

   logical (log_kind) :: update_time

!-----------------------------------------------------------------------
!
!  if tavg requested and tavg not already turned on, check to see
!  if it is time to start time averaging
!
!-----------------------------------------------------------------------

   if (.not. ltavg_on .and. tavg_freq_iopt /= freq_opt_never) then

      ltavg_on = time_to_start(tavg_start_iopt, tavg_start)
      call time_stamp('now','ymd',date_string=beg_date)
      if (tavg_freq_iopt == freq_opt_nstep) &
         write(beg_date,'(i10)') nsteps_total

      !*** if it is time to start, make sure requested fields
      !*** get triggered by the requested function

      if (ltavg_on) then
         do n=1,num_avail_tavg_fields
            if (avail_tavg_fields(n)%buf_loc < 0) &
                avail_tavg_fields(n)%buf_loc =    &
                abs(avail_tavg_fields(n)%buf_loc)
         end do
      endif
   endif

!-----------------------------------------------------------------------
!
!  setup time step and total integrated time for time average
!  adjust for averaging timesteps: if this is an averaging timestep,
!  the past values only contribute for 1/4 of a step and the
!  values for the step just before an averaging timestep contribute
!  for 1 1/4 steps.
!
!-----------------------------------------------------------------------

   update_time = .true.
   if (present(flagonly)) then
      if (flagonly) update_time = .false.
   endif

   if (ltavg_on .and. update_time) then
      if (avg_ts .or. back_to_back) then
         dtavg = p5*dtt
      else
         dtavg = dtt
      endif

      tavg_sum = tavg_sum + dtavg
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_set_flag

!***********************************************************************
!BOP
! !IROUTINE: write_tavg
! !INTERFACE:

 subroutine write_tavg(restart_type)

! !DESCRIPTION:
!  This routine writes requested tavg fields to a file.  The fields are
!  normalized by the time interval before writing.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      restart_type           ! tells tavg whether to write restart

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      nu,          &! i/o unit for output file
      iblock,      &! dummy block index
      nfield,      &! dummy field index
      loc,         &! buffer location for field
      io_phase      !'define' or 'write'

   character (char_len) ::  &
      file_suffix,          &! suffix to append to tavg file name
      hist_string,          &! string containing file history
      tavg_filename,        &! filename for tavg data
      tavg_pointer_file      ! filename for pointer file containing
                             !   location/name of last restart file

   character (8) :: &
      date_created   ! string with (real) date this file created

   character (10) :: &
      time_created   ! string with (real) date this file created

   type (io_field_desc), dimension(:), allocatable :: &
      tavg_fields

   type (io_dim) :: &
      i_dim, j_dim, &! dimension descriptors for horiz dims
      k_dim          ! dimension descriptor  for vertical levels

   logical (log_kind) :: &
      ltavg_write,       &! time to write a file
      lreset_tavg         ! time to reset time averages (reg tavg dump)

!-----------------------------------------------------------------------
!
!  is it time to write a file - if yes, create a file suffix
!
!-----------------------------------------------------------------------

   ltavg_write = .false.
   lreset_tavg = .false.

   if (ltavg_on) then
      ltavg_write = check_time_flag(tavg_flag)

      !*** regular tavg dump
      if (ltavg_write) then
         call create_suffix_tavg(file_suffix)
         lreset_tavg = .true.
      endif

      !*** tavg restart
      if (trim(restart_type) /= 'none') then
         if (.not. ltavg_write) then
            ltavg_write = .true.

            select case (trim(restart_type))
            case('even')
               file_suffix = trim(runid)/&
                                         &/'.even'
            case('odd')
               file_suffix = trim(runid)/&
                                         &/'.odd'
            case('end')
               file_suffix = trim(runid)/&
                                         &/'.end'
            case default
               call create_suffix_tavg(file_suffix)
               file_suffix = trim(file_suffix)/&
                                               &/'.restart'
            end select
         endif
      endif
   endif

!-----------------------------------------------------------------------
!
!  do the rest only if it is time to do a tavg dump
!
!-----------------------------------------------------------------------

   if (ltavg_write) then

!-----------------------------------------------------------------------
!
!     compute global averages of tavg fields
!     do this before normalization
!
!-----------------------------------------------------------------------

      call tavg_global

!-----------------------------------------------------------------------
!
!     normalize time averages
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(nfield)
      do iblock=1,nblocks_clinic
         do nfield=1,tavg_bufsize_2d
            if (TAVG_BUF_2D_METHOD(nfield) == tavg_method_avg) then
               TAVG_BUF_2D(:,:,  iblock,nfield) = &
               TAVG_BUF_2D(:,:,  iblock,nfield)/tavg_sum
            endif
         end do
         do nfield=1,tavg_bufsize_3d
            if (TAVG_BUF_3D_METHOD(nfield) == tavg_method_avg) then
               TAVG_BUF_3D(:,:,:,iblock,nfield) = &
               TAVG_BUF_3D(:,:,:,iblock,nfield)/tavg_sum
            endif
         end do
      end do
      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     AMOC strength calculation
!
!-----------------------------------------------------------------------
      if (tavg_do_MOC_Strength) then
        if (my_task == master_task) then
          write(*,*) "Averaging time: ", tavg_sum
        endif
        call tavg_amoc_diags
      endif

!-----------------------------------------------------------------------
!
!     create data file descriptor
!
!-----------------------------------------------------------------------

      if (.not. skip_tavg_dump) then

        call date_and_time(date=date_created, time=time_created)
        hist_string = char_blank
        write(hist_string,'(a23,a8,1x,a10)') &
           'POP TAVG file created: ',date_created,time_created

        tavg_file_desc = construct_file(tavg_fmt_out,                    &
                                     root_name  = trim(tavg_outfile),    &
                                     file_suffix= trim(file_suffix),     &
                                     title      ='POP TAVG file',        &
                                     conventions='POP TAVG conventions', &
                                     history    = trim(hist_string),     &
                                     record_length = rec_type_real,      &
                                     recl_words=nx_global*ny_global)

!-----------------------------------------------------------------------
!
!     add scalar fields to file as file attributes
!
!-----------------------------------------------------------------------

        call add_attrib_file(tavg_file_desc, 'tavg_sum'    , tavg_sum)
        call add_attrib_file(tavg_file_desc, 'nsteps_total', nsteps_total)
        call add_attrib_file(tavg_file_desc, 'tday'        , tday)
        call add_attrib_file(tavg_file_desc, 'iyear'       , iyear)
        call add_attrib_file(tavg_file_desc, 'imonth'      , imonth)
        call add_attrib_file(tavg_file_desc, 'iday'        , iday)
        call add_attrib_file(tavg_file_desc, 'beg_date'    , beg_date)

!-----------------------------------------------------------------------
!
!     open output file
!
!-----------------------------------------------------------------------

        call data_set (tavg_file_desc, 'open')

!-----------------------------------------------------------------------
!
!     write fields to file - this requires two phases
!     in this first phase, we define all the fields to be written
!
!-----------------------------------------------------------------------

        !*** define dimensions

        i_dim = construct_io_dim('i',nx_global)
        j_dim = construct_io_dim('j',ny_global)
        k_dim = construct_io_dim('k',km)

        allocate(tavg_fields(num_avail_tavg_fields))

        do nfield = 1,num_avail_tavg_fields  ! check all available fields

           loc = avail_tavg_fields(nfield)%buf_loc ! locate field in buffer

           if (loc > 0) then  ! field is actually requested and in buffer

              !*** construct io_field descriptors for each field

              if (avail_tavg_fields(nfield)%ndims == 2) then

                 tavg_fields(nfield) = construct_io_field(               &
                                avail_tavg_fields(nfield)%short_name,    &
                                i_dim, j_dim,                            &
                      long_name=avail_tavg_fields(nfield)%long_name,     &
                      units    =avail_tavg_fields(nfield)%units    ,     &
                      grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                     field_loc =avail_tavg_fields(nfield)%field_loc,     &
                    field_type =avail_tavg_fields(nfield)%field_type,    &
                  missing_value=avail_tavg_fields(nfield)%missing_value, &
                    valid_range=avail_tavg_fields(nfield)%valid_range,   &
                     r2d_array =TAVG_BUF_2D(:,:,:,loc) )

              else if (avail_tavg_fields(nfield)%ndims == 3) then

                 tavg_fields(nfield) = construct_io_field(               &
                                avail_tavg_fields(nfield)%short_name,    &
                                i_dim, j_dim, dim3=k_dim,                &
                      long_name=avail_tavg_fields(nfield)%long_name,     &
                      units    =avail_tavg_fields(nfield)%units    ,     &
                      grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                     field_loc =avail_tavg_fields(nfield)%field_loc,     &
                    field_type =avail_tavg_fields(nfield)%field_type,    &
                  missing_value=avail_tavg_fields(nfield)%missing_value, &
                    valid_range=avail_tavg_fields(nfield)%valid_range,   &
                     r3d_array =TAVG_BUF_3D(:,:,:,:,loc) )

              endif

              call data_set (tavg_file_desc, 'define', tavg_fields(nfield))
           endif
        end do

!-----------------------------------------------------------------------
!
!     write fields to file
!     in this second phase, we actually write the data for all the fields
!     after writing a field, the field descriptor is destroyed and the
!     file can be closed
!
!-----------------------------------------------------------------------

        do nfield = 1,num_avail_tavg_fields  ! check all available fields

           loc = avail_tavg_fields(nfield)%buf_loc ! locate field in buffer

           if (loc > 0) then  ! field is actually requested and in buffer
              call data_set (tavg_file_desc, 'write', tavg_fields(nfield))
              call destroy_io_field(tavg_fields(nfield))
           endif
        end do

        deallocate(tavg_fields)
        call data_set (tavg_file_desc, 'close')

        if (my_task == master_task) then
           write(stdout,blank_fmt)
           write(stdout,*) 'Wrote file: ', trim(tavg_file_desc%full_name)
        endif

!-----------------------------------------------------------------------
!
!     if pointer files are used, write tavg filenames to pointer file
!     do this only for tavg restarts - not tavg dumps
!
!-----------------------------------------------------------------------

        if (luse_pointer_files .and. .not. lreset_tavg) then
           call get_unit(nu)
           if (my_task == master_task) then
              tavg_pointer_file = trim(pointer_filename)/&
                                                         &/'.tavg'

              open(nu,file=tavg_pointer_file,form='formatted', &
                      status='unknown')
              write(nu,'(a)') trim(tavg_file_desc%full_name)
              close(nu)
           endif
           call release_unit(nu)
        endif

      endif ! test on skip_tavg_dump

!-----------------------------------------------------------------------
!
!     reset time averages if this is a regular tavg dump (as opposed
!     to a restart dump)  if this is a restart dump, unnormalize
!     in case a normal restart dump is written on the same timestep
!
!-----------------------------------------------------------------------

      if (lreset_tavg) then
         tavg_sum = c0
         call time_stamp('now','ymd',date_string=beg_date)
         if (tavg_freq_iopt == freq_opt_nstep) &
            write(beg_date,'(i10)') nsteps_total

         !$OMP PARALLEL DO PRIVATE(nfield)
         do iblock=1,nblocks_clinic
            do nfield=1,tavg_bufsize_2d
               select case (TAVG_BUF_2D_METHOD(nfield))
               case (tavg_method_avg)
                  TAVG_BUF_2D(:,:,  iblock,nfield) = c0
               case (tavg_method_min)
                  TAVG_BUF_2D(:,:,  iblock,nfield) = bignum
               case (tavg_method_max)
                  TAVG_BUF_2D(:,:,  iblock,nfield) = -bignum
               case default
                  TAVG_BUF_2D(:,:,  iblock,nfield) = c0
               end select
            end do
            do nfield=1,tavg_bufsize_3d
               select case (TAVG_BUF_3D_METHOD(nfield))
               case (tavg_method_avg)
                  TAVG_BUF_3D(:,:,:,iblock,nfield) = c0
               case (tavg_method_min)
                  TAVG_BUF_3D(:,:,:,iblock,nfield) = bignum
               case (tavg_method_max)
                  TAVG_BUF_3D(:,:,:,iblock,nfield) = -bignum
               case default
                  TAVG_BUF_3D(:,:,:,iblock,nfield) = c0
               end select
            end do
         end do
         !$OMP END PARALLEL DO

      else

         !$OMP PARALLEL DO PRIVATE(nfield)
         do iblock=1,nblocks_clinic
            do nfield=1,tavg_bufsize_2d
               if (TAVG_BUF_2D_METHOD(nfield) == tavg_method_avg) then
                  TAVG_BUF_2D(:,:,  iblock,nfield) = &
                  TAVG_BUF_2D(:,:,  iblock,nfield)*tavg_sum
               endif
            end do
            do nfield=1,tavg_bufsize_3d
               if (TAVG_BUF_3D_METHOD(nfield) == tavg_method_avg) then
                  TAVG_BUF_3D(:,:,:,iblock,nfield) = &
                  TAVG_BUF_3D(:,:,:,iblock,nfield)*tavg_sum
               endif
            end do
         end do
         !$OMP END PARALLEL DO

      endif

!-----------------------------------------------------------------------
!
!     get rid of file descriptor
!
!-----------------------------------------------------------------------

      call destroy_file(tavg_file_desc)
   endif ! lwrite_tavg

!-----------------------------------------------------------------------
!EOC

 end subroutine write_tavg

!***********************************************************************
!BOP
! !IROUTINE: read_tavg
! !INTERFACE:

 subroutine read_tavg

! !DESCRIPTION:
!  This routine reads a time average restart dump to continue
!  running time averages of requested fields.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::     &
     nu,               &   ! i/o unit
     iblock,           &   ! dummy block index
     n,                &   ! dummy for indexing character string
     in_fields,        &   ! num of fields in restart file
     nfield,           &   ! dummy field counter
     hdr_error,        &   ! error file for reading restart hdr
     in_nsteps_total,  &   ! nsteps_total according to tavg file
     in_iyear,         &   ! iyear according to tavg file
     in_imonth,        &   ! imonth according to tavg file
     in_iday,          &   ! iday according to tavg file
     loc,              &   ! buffer location
     errVal                ! internal error flag

   real (r8) ::        &
     in_tday               ! tday according to tavg file

   logical (log_kind  ) ::  &
     file_exists

   character (char_len) ::  &
     string,            &
     header_filename,   &! filename for restart contents
     char_temp,         &! for string manipulation
     tavg_pointer_file   ! filename for pointer file containing
                         !   location/name of last restart file

   type (io_field_desc), dimension(:), allocatable :: &
      tavg_fields          ! io field description for each field in file

   type (io_dim) :: &
      i_dim, j_dim, &! dimension descriptors for horiz dims
      k_dim          ! dimension descriptor  for vertical levels

!-----------------------------------------------------------------------
!
!  if pointer files are used, pointer file and must be read to get
!  actual filenames
!
!-----------------------------------------------------------------------

   errVal = 0

   call get_unit(nu)

   if (luse_pointer_files) then

      if (my_task == master_task) then
         tavg_pointer_file = char_blank
         tavg_pointer_file = trim(pointer_filename)/&
                                                   &/'.tavg'
         write(stdout,*) 'Reading pointer file: ', &
                         trim(tavg_pointer_file)
         inquire(file=trim(tavg_pointer_file),exist=file_exists)
         if (file_exists) then
            open(nu, file=trim(tavg_pointer_file), form='formatted', &
                     status='old')
            read(nu,'(a)') tavg_infile
            close(nu)
         else
            errVal = -1
         endif

      endif
      call broadcast_scalar(errVal, master_task)

      if (errVal /= 0) then
        write(stdout,*) char_blank
        write(string,*) 'FATAL ERROR: tavg_rpointer_file does not exist. pop2 model will exit'
        write(stdout,*) string
#ifdef CCSMCOUPLED
         call shr_sys_flush(stdout)
#endif
        call exit_POP(sigAbort,'(read_tavg) ERROR: tavg_rpointer_file does not exist ')
      endif

      call broadcast_scalar(tavg_infile, master_task)

   endif

   call release_unit(nu)

!-----------------------------------------------------------------------
!
!  define input file
!
!-----------------------------------------------------------------------

   tavg_file_desc = construct_file (tavg_fmt_in,                   &
                                    full_name=trim(tavg_infile),   &
                                    record_length = rec_type_real, &
                                    recl_words=nx_global*ny_global)

!-----------------------------------------------------------------------
!
!  define scalar fields in file as file attributes to be read during
!  open
!
!-----------------------------------------------------------------------

   call add_attrib_file(tavg_file_desc, 'tavg_sum'    , tavg_sum)
   call add_attrib_file(tavg_file_desc, 'nsteps_total', nsteps_total)
   call add_attrib_file(tavg_file_desc, 'tday'        , tday)
   call add_attrib_file(tavg_file_desc, 'iyear'       , iyear)
   call add_attrib_file(tavg_file_desc, 'imonth'      , imonth)
   call add_attrib_file(tavg_file_desc, 'iday'        , iday)
   call add_attrib_file(tavg_file_desc, 'beg_date'    , beg_date)

!-----------------------------------------------------------------------
!
!  open input file
!  this will also extract scalar variables which are defined as
!  file attributes
!
!-----------------------------------------------------------------------

   call data_set (tavg_file_desc, 'open_read')

   call extract_attrib_file(tavg_file_desc, 'nsteps_total', &
                                          in_nsteps_total)
   call extract_attrib_file(tavg_file_desc, 'tavg_sum', tavg_sum)
   call extract_attrib_file(tavg_file_desc, 'beg_date', beg_date)
   !call extract_attrib_file(tavg_file_desc, 'tday', in_tday)
   !call extract_attrib_file(tavg_file_desc, 'iyear', in_iyear)
   !call extract_attrib_file(tavg_file_desc, 'imonth', in_imonth)
   !call extract_attrib_file(tavg_file_desc, 'iday', in_iday)

   !*** check nsteps total for validity
   if (in_nsteps_total /= nsteps_total) then
      write(stdout,'(i6,a29,i6,a35)') &
         in_nsteps_total,' nsteps_total in tavg restart', &
         nsteps_total,   ' nsteps_total in current simulation'
#ifdef CCSMCOUPLED
         call shr_sys_flush(stdout)
#endif
      call exit_POP(sigAbort,'TAVG:restart file has wrong time step?')
   endif

!-----------------------------------------------------------------------
!
!  define requested fields to read in from file
!  NOTE: This requires that the tavg_contents file is consistent
!  with the tavg restart file.  There are currently no checks on this.
!
!-----------------------------------------------------------------------

   !*** define dimensions

   i_dim = construct_io_dim('i',nx_global)
   j_dim = construct_io_dim('j',ny_global)
   k_dim = construct_io_dim('k',km)

   allocate(tavg_fields(num_avail_tavg_fields))

   do nfield = 1,num_avail_tavg_fields
      loc = avail_tavg_fields(nfield)%buf_loc
      if (loc > 0) then
         if (avail_tavg_fields(nfield)%ndims == 2) then

            tavg_fields(nfield) = construct_io_field(                &
                            avail_tavg_fields(nfield)%short_name,    &
                            i_dim, j_dim,                            &
                  long_name=avail_tavg_fields(nfield)%long_name,     &
                  units    =avail_tavg_fields(nfield)%units    ,     &
                  grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                 field_loc =avail_tavg_fields(nfield)%field_loc,     &
                field_type =avail_tavg_fields(nfield)%field_type,    &
              missing_value=avail_tavg_fields(nfield)%missing_value, &
                valid_range=avail_tavg_fields(nfield)%valid_range,   &
                 r2d_array =TAVG_BUF_2D(:,:,:,loc) )

         else if (avail_tavg_fields(nfield)%ndims == 3) then

            tavg_fields(nfield) = construct_io_field(                &
                            avail_tavg_fields(nfield)%short_name,    &
                            i_dim, j_dim, dim3=k_dim,                &
                  long_name=avail_tavg_fields(nfield)%long_name,     &
                  units    =avail_tavg_fields(nfield)%units    ,     &
                  grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                 field_loc =avail_tavg_fields(nfield)%field_loc,     &
                field_type =avail_tavg_fields(nfield)%field_type,    &
              missing_value=avail_tavg_fields(nfield)%missing_value, &
                valid_range=avail_tavg_fields(nfield)%valid_range,   &
                 r3d_array =TAVG_BUF_3D(:,:,:,:,loc) )

         endif

         call data_set (tavg_file_desc, 'define', tavg_fields(nfield))
      endif
   end do

!-----------------------------------------------------------------------
!
!  now we actually read each field
!  after reading, get rid of io field descriptors and close file
!
!-----------------------------------------------------------------------

   do nfield = 1,num_avail_tavg_fields
      loc = avail_tavg_fields(nfield)%buf_loc
      if (loc > 0) then
         call data_set (tavg_file_desc, 'read', tavg_fields(nfield))
         call destroy_io_field(tavg_fields(nfield))
      endif
   end do

   deallocate(tavg_fields)
   call data_set (tavg_file_desc, 'close')

   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) ' file read: ', tavg_infile
   endif

   call destroy_file(tavg_file_desc)
   call release_unit(nu)

!-----------------------------------------------------------------------
!
!  de-normalize sums
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(nfield)
   do iblock=1,nblocks_clinic
      do nfield=1,tavg_bufsize_2d
         if (TAVG_BUF_2D_METHOD(nfield) == tavg_method_avg) then
            TAVG_BUF_2D(:,:,  iblock,nfield) = &
            TAVG_BUF_2D(:,:,  iblock,nfield)*tavg_sum
         endif
      end do
      do nfield=1,tavg_bufsize_3d
         if (TAVG_BUF_3D_METHOD(nfield) == tavg_method_avg) then
            TAVG_BUF_3D(:,:,:,iblock,nfield) = &
            TAVG_BUF_3D(:,:,:,iblock,nfield)*tavg_sum
         endif
      end do
   end do
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------

   call tavg_global   ! print global sums of time averages

!-----------------------------------------------------------------------
!EOC

 end subroutine read_tavg

!***********************************************************************
!BOP
! !IROUTINE: tavg_global
! !INTERFACE:

 subroutine tavg_global

! !DESCRIPTION:
!  Calculates and print global integrals of time average fields
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::     &
      k,               &   ! vertical level index
      ifield,          &   ! field identifier
      iblock,          &   ! block index
      nfield,          &   ! dummy field index
      field_loc,       &   ! field location (center,Nface,Eface,NEcorner)
      field_type           ! field type (scalar, vector, angle)

   real (r8) ::        &
      tavg_field_sum,  &   ! sum of tavg field
      tavg_norm            ! normalization for average

   real (r8), dimension (:,:,:), allocatable ::  &
      WORK               ! temp for holding area_weighted field

   real (r8), dimension (:,:), allocatable ::  &
      RMASK              ! topography mask for global sum

!-----------------------------------------------------------------------
!
!  calculate globally-integrated time average of each chosen 2d field
!
!-----------------------------------------------------------------------

   allocate (RMASK(nx_block,ny_block), &
             WORK (nx_block,ny_block,nblocks_clinic))

   if (my_task == master_task) then
     write (stdout,blank_fmt)
     write (stdout,'(a22)') 'Global Time Averages: '
   endif

   do nfield=1,num_avail_tavg_fields
      ifield = avail_tavg_fields(nfield)%buf_loc
      if (ifield > 0) then

         field_loc  = avail_tavg_fields(nfield)%field_loc
         field_type = avail_tavg_fields(nfield)%field_type

         if (avail_tavg_fields(nfield)%method == tavg_method_avg) then
            tavg_norm = tavg_sum
         else
            tavg_norm = c1
         endif

         !*** 2-d fields

         if (avail_tavg_fields(nfield)%ndims == 2) then

            !$OMP PARALLEL DO
            do iblock = 1,nblocks_clinic
               select case(field_loc)
               case(field_loc_center)
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    TAREA(:,:,iblock)*RCALCT(:,:,iblock)
               case(field_loc_NEcorner)
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    UAREA(:,:,iblock)*RCALCU(:,:,iblock)
               case default ! make U cell the default for all other cases
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    UAREA(:,:,iblock)*RCALCU(:,:,iblock)
               end select
            end do
            !$OMP END PARALLEL DO

            tavg_field_sum = global_sum(WORK, distrb_clinic, field_loc)

            select case(field_loc)
            case(field_loc_center)
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_t)
            case(field_loc_NEcorner)
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_u)
            case default ! make U cell the default for all other cases
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_u)
            end select

         !*** 3-d fields

         else

            !$OMP PARALLEL DO PRIVATE(k)
            do iblock = 1,nblocks_clinic
               WORK(:,:,iblock) = c0

               select case(field_loc)

               case(field_loc_center)
                  do k=1,km
                     RMASK = merge(c1, c0, k <= KMT(:,:,iblock))
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        TAREA(:,:,iblock)*RMASK
                  end do

               case(field_loc_NEcorner)
                  do k=1,km
                     RMASK = merge(c1, c0, k <= KMU(:,:,iblock))
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        UAREA(:,:,iblock)*RMASK
                  end do

               case default ! make U cell the default for all other cases
                  do k=1,km
                     RMASK = merge(c1, c0, k <= KMU(:,:,iblock))
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        UAREA(:,:,iblock)*RMASK
                  end do

               end select
            end do
            !$OMP END PARALLEL DO

            tavg_field_sum = global_sum(WORK, distrb_clinic, field_loc)

            select case(field_loc)
            case(field_loc_center)
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_t)
            case(field_loc_NEcorner)
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_u)
            case default ! make U cell the default for all other cases
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_u)
            end select

         endif

         if (my_task == master_task) then
            write (stdout,*) trim(avail_tavg_fields(nfield)%short_name), &
                             ': ', tavg_field_sum
         endif
      endif
   end do

   deallocate (RMASK, WORK)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_global

!***********************************************************************
!BOP
! !IROUTINE: accumulate_tavg_field
! !INTERFACE:

 subroutine accumulate_tavg_field(ARRAY,field_id,block,k)

! !DESCRIPTION:
!  This routine updates a tavg field.  If the time average of the
!  field is requested, it accumulates a time sum of a field by
!  multiplying by the time step and accumulating the sum into the
!  tavg buffer array.  If the min or max of a field is requested, it
!  checks the current value and replaces the min, max if the current
!  value is less than or greater than the stored value.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block,           &! local block address (in baroclinic distribution)
      k,               &! vertical level
      field_id          ! index into available fields for tavg field info

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      ARRAY             ! array of data for this block to add to
                        !  accumulated sum in tavg buffer
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bufloc,            &! location of field in tavg buffer
      ndims               ! rank of field (2=2d,3=3d)

!-----------------------------------------------------------------------
!
!  get buffer location and field info from avail_tavg_field array
!
!-----------------------------------------------------------------------

   bufloc = avail_tavg_fields(field_id)%buf_loc
   if (bufloc <= 0) &
     call exit_POP(sigAbort, &
                    'tavg: attempt to accumulate bad tavg field')

   ndims = avail_tavg_fields(field_id)%ndims

!-----------------------------------------------------------------------
!
!  update the field into the tavg buffer
!
!-----------------------------------------------------------------------

   select case (avail_tavg_fields(field_id)%method)

   case (tavg_method_avg)  ! accumulate running time sum for time avg
      if (ndims == 2) then
         TAVG_BUF_2D(:,:,block,bufloc) = &
         TAVG_BUF_2D(:,:,block,bufloc) + dtavg*ARRAY
      else
         TAVG_BUF_3D(:,:,k,block,bufloc) = &
         TAVG_BUF_3D(:,:,k,block,bufloc) + dtavg*ARRAY
      endif

   case (tavg_method_min)  ! replace with current minimum value
      if (ndims == 2) then
         where (ARRAY < TAVG_BUF_2D(:,:,block,bufloc))
            TAVG_BUF_2D(:,:,block,bufloc) = ARRAY
         end where
      else
         where (ARRAY < TAVG_BUF_3D(:,:,k,block,bufloc))
            TAVG_BUF_3D(:,:,k,block,bufloc) = ARRAY
         end where
      endif

   case (tavg_method_max)  ! replace with current minimum value
      if (ndims == 2) then
         where (ARRAY > TAVG_BUF_2D(:,:,block,bufloc))
            TAVG_BUF_2D(:,:,block,bufloc) = ARRAY
         end where
      else
         where (ARRAY > TAVG_BUF_3D(:,:,k,block,bufloc))
            TAVG_BUF_3D(:,:,k,block,bufloc) = ARRAY
         end where
      endif

   case default
   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine accumulate_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: define_tavg_field
! !INTERFACE:

 subroutine define_tavg_field(id, short_name, ndims, tavg_method, &
                                  long_name, units, &
                                  grid_loc, missing_value, valid_range, &
                                  field_loc, field_type)

! !DESCRIPTION:
!  Initializes description of an available field and returns location
!  in the available fields array for use in later tavg calls.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      id                ! location in avail_fields array for use in
                        ! later tavg routines

! !INPUT PARAMETERS:

   character(*), intent(in) :: &
      short_name               ! short name for field

   integer (i4), intent(in) :: &
      ndims                     ! number of dims (2 or 3) of tavg field

   integer (i4), intent(in), optional :: &
      field_loc,              &! location in grid
      field_type,             &! type of field (scalar, vector, angle)
      tavg_method              ! id for method of averaging
                               ! default is tavg_method_avg

   character(*), intent(in), optional :: &
      long_name,              &! long descriptive name for field
      units                    ! physical units for field

   character(4), intent(in), optional :: &
      grid_loc                 ! location in grid (in 4-digit code)

   real (r4), intent(in), optional :: &
      missing_value            ! value on land pts

   real (r4), dimension(2), intent(in), optional :: &
      valid_range              ! min/max

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  increment the number of defined fields and make sure it does not
!  exceed the maximum
!  return the id as the current number
!
!-----------------------------------------------------------------------

   num_avail_tavg_fields = num_avail_tavg_fields + 1
   if (num_avail_tavg_fields > max_avail_tavg_fields) then
      call exit_POP(sigAbort, &
                    'tavg: defined tavg fields > max allowed')
   endif

   id = num_avail_tavg_fields

!-----------------------------------------------------------------------
!
!  now fill the field descriptor
!
!-----------------------------------------------------------------------

   avail_tavg_fields(id)%ndims      = ndims
   avail_tavg_fields(id)%short_name = short_name
   avail_tavg_fields(id)%buf_loc    = 0  ! will be reset later

   if (present(long_name)) then
      avail_tavg_fields(id)%long_name = long_name
   else
      avail_tavg_fields(id)%long_name = char_blank
   endif

   if (present(units)) then
      avail_tavg_fields(id)%units = units
   else
      avail_tavg_fields(id)%units = char_blank
   endif

   if (present(grid_loc)) then
      avail_tavg_fields(id)%grid_loc = grid_loc
   else
      avail_tavg_fields(id)%grid_loc = '    '
   endif

   if (present(tavg_method)) then
      avail_tavg_fields(id)%method = tavg_method
   else
      avail_tavg_fields(id)%method = tavg_method_avg
   endif

   if (present(missing_value)) then
      avail_tavg_fields(id)%missing_value = missing_value
   else
      avail_tavg_fields(id)%missing_value = undefined
   endif

   if (present(valid_range)) then
      avail_tavg_fields(id)%valid_range = valid_range
   else
      avail_tavg_fields(id)%valid_range = undefined
   endif

   !*** set field location, field type used by i/o, ghost cell update
   !*** and other communication routines.  because ghost cells for tavg
   !*** fields are not typically used, the default is field_xxx_noupdate

   if (present(field_loc)) then
      avail_tavg_fields(id)%field_loc = field_loc
   else
      !*** try to decode field location from grid_loc
      if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '1') then
         avail_tavg_fields(id)%field_loc = field_loc_center
      else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '2') then
         avail_tavg_fields(id)%field_loc = field_loc_NEcorner
      else if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '2') then
         avail_tavg_fields(id)%field_loc = field_loc_Nface
      else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '1') then
         avail_tavg_fields(id)%field_loc = field_loc_Eface
      else
         avail_tavg_fields(id)%field_loc = field_loc_noupdate
      endif
   endif

   if (present(field_type)) then
      avail_tavg_fields(id)%field_type = field_type
   else
      avail_tavg_fields(id)%field_type = field_type_noupdate
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine define_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: request_tavg_field
! !INTERFACE:

 subroutine request_tavg_field(short_name)

! !DESCRIPTION:
!  This field marks an available field as requested and computes
!  the location in the tavg buffer array.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      short_name                ! the short name of the field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! loop index
      id                  ! location of field in avail_fields array

!-----------------------------------------------------------------------
!
!  search for field with same name
!
!-----------------------------------------------------------------------

   id = 0
   srch_loop: do n=1,num_avail_tavg_fields
      if (trim(avail_tavg_fields(n)%short_name) == short_name) then
         id = n
         exit srch_loop
      endif
   end do srch_loop

   if (id == 0) then
      if (my_task == master_task) write(stdout,*) 'Requested ', &
                                                  trim(short_name)
      call exit_POP(sigAbort,'tavg: requested field unknown')
   endif

!-----------------------------------------------------------------------
!
!  set the position in the buffer and advance the buffer position
!  for the next field
!
!-----------------------------------------------------------------------

   if (avail_tavg_fields(id)%ndims == 2) then
      tavg_bufsize_2d = tavg_bufsize_2d + 1
      avail_tavg_fields(id)%buf_loc = tavg_bufsize_2d
   else if (avail_tavg_fields(id)%ndims == 3) then
      tavg_bufsize_3d = tavg_bufsize_3d + 1
      avail_tavg_fields(id)%buf_loc = tavg_bufsize_3d
   endif

!-----------------------------------------------------------------------
!
!  if tavg is on, but not started yet, set the buf_loc to -buf_loc
!  to show that it is requested, but will not return true for
!  requested_tavg_field
!
!-----------------------------------------------------------------------

   if (.not. ltavg_on) then
      avail_tavg_fields(id)%buf_loc = -avail_tavg_fields(id)%buf_loc
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine request_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: tavg_requested
! !INTERFACE:

 function tavg_requested(id)

! !DESCRIPTION:
!  This function determines whether an available (defined) tavg field
!  has been requested by a user (through the input contents file) and
!  returns true if it has.  Note that if tavg has been turned off,
!  the function will always return false.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      id                   ! id returned by the define function which
                           !   gives the location of the field

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      tavg_requested     ! result of checking whether the field has
                         !   been requested

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check the buffer location - if zero, the field has not been
!  requested
!
!-----------------------------------------------------------------------

   if (id < 1 .or. id > num_avail_tavg_fields) then
      call exit_POP(sigAbort,'tavg_requested: invalid tavg id')
   endif

   if (avail_tavg_fields(id)%buf_loc > 0) then
      tavg_requested = .true.
   else
      tavg_requested = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function tavg_requested

!***********************************************************************
!BOP
! !IROUTINE: create_suffix_tavg
! !INTERFACE:

 subroutine create_suffix_tavg(file_suffix)

! !DESCRIPTION:
!  Creates a suffix to append to output filename based on frequency
!  option and averaging interval.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   character (char_len), intent(out) :: &
      file_suffix           ! suffix to append to root filename

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variable
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      cindx1, cindx2,    &! indices into character strings
      len_date            ! length of date string

   character (char_len) :: &
      char_temp            ! temp character space (for removing spaces)

   character (10) :: &
      cstep_beg,     &! beginning step  of this particular average
      cstep_end,     &! ending    step  of this particular average
      cdate           ! character string with yyyymmdd and optional
                      ! separator (eg yyyy-mm-dd)

   character (4) :: &
      cyear_beg,    &! beginning year  of this particular average
      cyear_end      ! end       year  of this particular average

   character (2) :: &
      cmonth_beg,   &! beginning month of this particular average
      cmonth_end,   &! end       month of this particular average
      cday_beg,     &! beginning day   of this particular average
      cday_end       ! end       day   of this particular average

!-----------------------------------------------------------------------
!
!  start suffix with runid
!
!-----------------------------------------------------------------------

   file_suffix = char_blank
   cindx2 = len_trim(runid) + 1
   file_suffix(1:cindx2) = trim(runid)/&
                                       &/'.'
   cindx1 = cindx2 + 1

!-----------------------------------------------------------------------
!
!  extract beginning year, month, day or time step from beg_date
!  and determine end date
!
!-----------------------------------------------------------------------

   cdate = adjustl(beg_date)

   !***
   !*** use step numbers if tavg freq option is nstep
   !***

   cstep_beg  = trim(cdate) ! in case beg_date actually step number
   write(cstep_end,'(i10)') nsteps_total - 1
   cdate  = adjustl(cstep_end)
   cstep_end = trim(cdate)

   call time_stamp('last','ymd',date_string=cdate)  ! last date

   if (date_separator == ' ') then  ! no date separator
      cyear_beg  = beg_date(1:4)
      cmonth_beg = beg_date(5:6)
      cday_beg   = beg_date(7:8)
      cyear_end  = cdate(1:4)
      cmonth_end = cdate(5:6)
      cday_end   = cdate(7:8)
   else
      cyear_beg  = beg_date(1:4)
      cmonth_beg = beg_date(6:7)
      cday_beg   = beg_date(9:10)
      cyear_end  = cdate(1:4)
      cmonth_end = cdate(6:7)
      cday_end   = cdate(9:10)
   endif

!-----------------------------------------------------------------------
!
!  create time portion of suffix based on frequency option
!  note that concatenation operator split across lines to avoid
!   problems with some cpp preprocessors
!
!-----------------------------------------------------------------------

   select case (tavg_freq_iopt)
   case (freq_opt_nyear)
      if (tavg_freq == 1) then
         cindx2 = cindx1 + 3
         file_suffix(cindx1:cindx2) = cyear_end
      else
         cindx2 = cindx1 + 8
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                                &/'-'/&
                                                      &/cyear_end
      endif

   case (freq_opt_nmonth)
      if (tavg_freq == 1) then
         cindx2 = cindx1 + 5
         file_suffix(cindx1:cindx2) = cyear_end/&
                                    &/cmonth_end
      else
         cindx2 = cindx1 + 12
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                    &/cmonth_beg/&
                                    &/'-'/&
                                    &/cyear_end/&
                                    &/cmonth_end
      endif

   case (freq_opt_nday)
      if (tavg_freq == 1) then
         cindx2 = cindx1 + 7
         file_suffix(cindx1:cindx2) = cyear_end/&
                                    &/cmonth_end/&
                                    &/cday_end
      else
         cindx2 = cindx1 + 16
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                    &/cmonth_beg/&
                                    &/cday_beg/&
                                    &/'-'/&
                                    &/cyear_end/&
                                    &/cmonth_end/&
                                    &/cday_end
      endif

   case (freq_opt_nstep)
      cindx2 = cindx1 + len_trim(cstep_beg) + len_trim(cstep_end)
      file_suffix(cindx1:cindx2) = trim(cstep_beg)/&
                                 &/'-'/&
                                 &/trim(cstep_end)

   case default  ! use nstep for other options
      cindx2 = len_trim(cstep_beg) + len_trim(cstep_end) + 1
      file_suffix(cindx1:cindx2) = trim(cstep_beg)/&
                                 &/'-'/&
                                 &/trim(cstep_end)

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine create_suffix_tavg

!***********************************************************************

!***********************************************************************
!BOP
! !IROUTINE: set_in_tavg_contents
! !INTERFACE:

 function set_in_tavg_contents(id)

! !DESCRIPTION:
!  This function determines whether a tavg field has been set in
!  the input contents file and returns true if it has.  This function is
!  different from tavg_requested in that ltavg_on status is irrelevent.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      id                   ! id returned by the define function which
                           !   gives the location of the field

! !OUTPUT PARAMETERS:

   logical (log_kind) ::      &
      set_in_tavg_contents    ! result of checking whether the field has
                               !   been requested

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check the buffer location - if not zero, then the field is in the
!  tavg contents file
!
!-----------------------------------------------------------------------

   if (id < 1 .or. id > num_avail_tavg_fields) then
     set_in_tavg_contents = .false.
   elseif (avail_tavg_fields(id)%buf_loc /= 0) then
     set_in_tavg_contents = .true.
   else
     set_in_tavg_contents = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function set_in_tavg_contents

!***********************************************************************
!BOP
! !IROUTINE: tavg_id
! !INTERFACE:

 function tavg_id(short_name,quiet)

! !DESCRIPTION:
!  This function determines whether a tavg field has been defined
!  by some module.  If so, then the id of that field is returned.
!  This function does not cause a model exit if the field has not
!  been defined; error-checking must be done by the calling routine.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)  :: &
      short_name                  ! the short name of the tavg field

   logical (log_kind), intent(in), optional :: &
      quiet                        ! do not print error message

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      tavg_id               ! id of the tavg field, if it exists

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::  &
     id,                  &
     n

   logical (log_kind) ::  &
     msg

   id = 0

   srch_loop: do n=1,num_avail_tavg_fields
      if (trim(avail_tavg_fields(n)%short_name) == trim(short_name)) then
         id = n
         exit srch_loop
      endif
   end do srch_loop

   msg = .true.
   if (present(quiet)) then
     if (quiet) msg = .false.
   endif

   if (id == 0 .and. msg) then
      if (my_task == master_task)  &
          write(stdout,*) 'Field ', trim(short_name), ' has not been defined.'
   endif

   tavg_id = id

!-----------------------------------------------------------------------
!EOC

 end function tavg_id

!***********************************************************************
!BOP
! !IROUTINE: init_tavg_amoc
! !INTERFACE:

  subroutine init_tavg_amoc

    ! !DESCRIPTION:
    !
    !  Compute atlantic meridional overturning circulation diagnostically from
    !  other tavg quantities
    !
    !
    ! !REVISION HISTORY:
    !  same as module

    integer (int_kind) ::   &
      i,j,k,ig,igp1,jg,n,right_idx,iblock,        & ! loop indices
      nu,            &
      ioerr

    logical :: is_found

    real (r8) :: &
      p_lat, p_lon, dlat, dlon, &
      long_deg, lat_deg      ! work variable for auxilary grid spacing (degrees)

    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      WORK

    type(block) :: &
       this_block  ! block info for current block
    !EOP
    !BOC

    !-------------------------------------------------------------------------
    ! Check for requireds field if AMOC strength is requested
    tavg_id_WVEL = tavg_id('WVEL')
    tavg_id_VVEL = tavg_id('VVEL')
    if (tavg_id_WVEL == 0 .or.  &
        tavg_id_VVEL == 0 )     then
      exit_string = 'FATAL ERROR: for amoc diagnostics, WVEL and VVEL must be requested in tavg_contents file'
      call exit_POP (sigAbort,exit_string)
    endif

    if (my_task == master_task) then
       write(stdout,'(A)') 'Computing AMOC strength from averaged quantities requested.'
    endif

    tavg_loc_WVEL = 0 ; tavg_loc_VVEL = 0
    tavg_loc_WVEL = abs(avail_tavg_fields(tavg_id_WVEL)%buf_loc)
    tavg_loc_VVEL = abs(avail_tavg_fields(tavg_id_VVEL)%buf_loc)

    ! Build auxiliary grid
    if (my_task == master_task) then
      allocate ( TLATD_G(nx_global,ny_global) )
    endif

    WORK = TLAT * radian
    call gather_global (TLATD_G, WORK, master_task, distrb_clinic)

    !-------------------------------------------------------------------------
    ! Allocate aux grid array
    allocate ( lat_aux_edge  (n_lat_aux_grid+1),  &
               lat_aux_center(n_lat_aux_grid  ) )

    dlat = (lat_aux_end - lat_aux_begin) / dble(n_lat_aux_grid)

    do j=1,n_lat_aux_grid+1
      lat_aux_edge(j) = lat_aux_begin + dble(j-1)*dlat
    enddo

    do j=1,n_lat_aux_grid
      lat_aux_center(j) = lat_aux_begin + p5*dlat + dble(j-1)*dlat
    enddo

    if ( my_task == master_task ) then
      allocate (TAVG_MOC_G(n_lat_aux_grid+1,km+1,1))
    endif

    !-------------------------------------------------------------------------
    ! Build atlantic mask
    allocate(ATLANTIC_MASK_LONG(n_lat_aux_grid))
    allocate(ATLANTIC_MASK_LONG_MIN(n_lat_aux_grid))
    allocate(ATLANTIC_MASK_LONG_MAX(n_lat_aux_grid))

    ! Read mask file data
    call get_unit(nu)
    if (my_task == master_task) then
       open(nu,file=atlantic_mask_file,status='old',form='formatted', iostat=ioerr)
    endif
    call broadcast_scalar(ioerr, master_task)
    if (ioerr /= 0) call exit_POP(sigAbort, 'Error opening atlantic_mask_file')

    if (my_task == master_task) then
       if (ioerr == 0) then ! successful open
          grid_read: do j = 1, n_lat_aux_grid
             read(nu,*,iostat=ioerr) ATLANTIC_MASK_LONG(j), ATLANTIC_MASK_LONG_MIN(j), ATLANTIC_MASK_LONG_MAX(j)
             if (ioerr /= 0) exit grid_read
          end do grid_read
          close(nu)
       endif
    endif
    call release_unit(nu)

    call broadcast_scalar(ioerr, master_task)
    if (ioerr /= 0) call exit_POP(sigAbort, 'Error reading atlantic_mask_file')

    call broadcast_array(ATLANTIC_MASK_LONG, master_task)
    call broadcast_array(ATLANTIC_MASK_LONG_MIN, master_task)
    call broadcast_array(ATLANTIC_MASK_LONG_MAX, master_task)

    if (my_task == master_task) then
      allocate (ATLANTIC_MASK_LAT_AUX(nx_global,ny_global) )
    endif

    !$OMP PARALLEL DO PRIVATE(iblock, this_block)
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)

       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             long_deg = TLON(i,j,iblock) * radian
             lat_deg = TLAT(i,j,iblock) * radian
             if (long_deg > 180.0) then
                long_deg = long_deg - 360.0
             endif

             if (lat_deg < ATLANTIC_MASK_LONG(1)) then
               right_idx = 1
             else if (lat_deg > ATLANTIC_MASK_LONG(n_lat_aux_grid)) then
               right_idx = n_lat_aux_grid-1
             else
               do n = 1, n_lat_aux_grid-1
                 if (lat_deg >= ATLANTIC_MASK_LONG(n) .and. &
                     lat_deg < ATLANTIC_MASK_LONG(n+1) ) then
                    right_idx = n
                    exit
                 endif
               enddo
             endif

             if (KMT(i,j,iblock) > 0 .and. &
                 lat_deg >= ATLANTIC_MASK_LONG(1) .and. &
                 lat_deg < ATLANTIC_MASK_LONG(n_lat_aux_grid) .and. &
                 long_deg >= ATLANTIC_MASK_LONG_MIN(right_idx) .and. &
                 long_deg <= ATLANTIC_MASK_LONG_MAX(right_idx+1))  then
                ATLANTIC_MASK_LAT(i,j,iblock) = 1
             else
                ATLANTIC_MASK_LAT(i,j,iblock) = 0
             endif
          enddo
       enddo
    enddo

    call gather_global (ATLANTIC_MASK_LAT_AUX(:,:), ATLANTIC_MASK_LAT(:,:,:),&
                        master_task,distrb_clinic)

    if (my_task == master_task) then
      !-------------------------------------------------------------------------
      ! Get data for final AMOC strength estimate location
      amoc_strength_target_lat_idx = -1
      do j = 1, n_lat_aux_grid+1
        if ( abs(lat_aux_edge(j)-26.0d0) < 0.01) then
           amoc_strength_target_lat_idx = j
           exit
        endif
      enddo
      if (amoc_strength_target_lat_idx < 0) then
        call exit_POP (sigAbort,"Can't find an aux grid latitude close to target 26 deg N")
      endif
      do k = 1, km
        if (zt(k) > 100000.0d0) then
          amoc_strength_target_depth_idx = k-1
          exit
        endif
      enddo

      !-------------------------------------------------------------------------
      ! Get the index of southern the southern boundary of the AMOC
      lat_aux_region_start = 0

      do j = 1, ny_global
        if (any(ATLANTIC_MASK_LAT_AUX(:,j) == 1)) then
          lat_aux_region_start = j
          exit
        endif
      enddo

      WRITE(*,*) " Southern j-index of the AMOC region ", lat_aux_region_start
    endif

    !-----------------------------------------------------------------------
    !EOC

  end subroutine init_tavg_amoc

!***********************************************************************
!BOP
! !IROUTINE: tavg_amoc_diags
! !INTERFACE:

  subroutine tavg_amoc_diags

    ! !DESCRIPTION:
    !
    !  Compute atlantic meridional overturning circulation diagnostically from
    !  other tavg quantities
    !
    !
    ! !REVISION HISTORY:
    !  same as module

    !-----------------------------------------------------------------------
    !
    !  begin computations
    !
    !-----------------------------------------------------------------------

    call compute_amoc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ))

    !-----------------------------------------------------------------------
    !EOC
  end subroutine tavg_amoc_diags

!***********************************************************************
!BOP
! !IROUTINE: compute_amoc
! !INTERFACE:
  subroutine compute_amoc ( W_E, V_E )

    ! !DESCRIPTION:
    ! This subroutine computes atlantic meridional overturning circulation
    !
    ! !REVISION HISTORY:
    !  same as module

    ! !INPUT PARAMETERS:
    !  input variables. the present design assumes that these are
    !  time-averaged inputs from tavg.F.
    !
    real (r4), dimension(:,:,:,:), intent(in) ::  &
      W_E,    &! Eulerian-mean vertical velocity component
      V_E      ! Eulerian-mean velocity component in the grid-y direction

!EOP
!BOC
    !-----------------------------------------------------------------------
    !
    !  local variables
    !
    !-----------------------------------------------------------------------
    integer (int_kind) ::  &
       i, j, k, n, iblock     ! loop indices

    real (r8), dimension(km,1) ::  &
       moc_s                  ! southern boundary values of moc

    real (r8), dimension(:,:,:), allocatable ::  &
       WORK1, WORK2           ! work arrays

    real (r8), dimension(:,:,:), allocatable ::  &
       WORK1_G                ! global work arrays

    allocate ( WORK1(nx_block,ny_block,nblocks_clinic), &
               WORK2(nx_block,ny_block,nblocks_clinic) )

    allocate ( WORK1_G(nx_global,ny_global,km) )

    ! Compute the w_e * area where there is ocean
    do k = 1, km
      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1,nblocks_clinic
         WORK1(:,:,iblock) = merge(W_E(:,:,k,iblock)*TAREA(:,:,iblock), c0, k <= KMT(:,:,iblock))
      enddo
      !$OMP END PARALLEL DO
      call gather_global (WORK1_G(:,:,k), WORK1, master_task,distrb_clinic)
    enddo

    if ( my_task == master_task )  then
      ! Reset AMOC data
      TAVG_MOC_G(:,:,1) = c0

      ! Integrate in longitude, accumulating on latitudes
      do n = 2, n_lat_aux_grid + 1
        TAVG_MOC_G(n,:,1) = TAVG_MOC_G(n-1,:,1)
        do j = 1, ny_global
          do i = 1, nx_global
            if ( TLATD_G(i,j) >= lat_aux_edge(n-1) .and.  &
                 TLATD_G(i,j) <  lat_aux_edge(n  ) .and.  &
                 ATLANTIC_MASK_LAT_AUX(i,j) == 1 ) then
              do k = 1, km
                TAVG_MOC_G(n,k,1) = TAVG_MOC_G(n,k,1) + WORK1_G(i,j,k)
              enddo
            endif
          enddo
        enddo
      enddo
    endif

    ! Compute Tgrid y-velocity*dx by averaging 2 Ugrid points
    do k = 1, km
      !$OMP PARALLEL DO PRIVATE(iblock,i,j)
      do iblock = 1, nblocks_clinic
        WORK1(:,:,iblock) = p5 * V_E(:,:,k,iblock) * DXU(:,:,iblock)
        do j=1,ny_block
          do i=2,nx_block
            WORK2(i,j,iblock) = WORK1(i-1,j,iblock)
          enddo
        enddo
        WORK2(1,:,iblock) = c0
        WORK1(:,:,iblock) = WORK1(:,:,iblock) + WORK2(:,:,iblock)
      enddo ! iblock
      !$OMP END PARALLEL DO
      call gather_global (WORK1_G(:,:,k), WORK1, master_task,distrb_clinic)
    enddo

    if ( my_task == master_task )  then
      moc_s = c0

      ! Compute sourthern boundary integral of V_E
      j = lat_aux_region_start
      do i = 1, nx_global
        if (ATLANTIC_MASK_LAT_AUX(i,j) == 1) then
          do k = 1, km
            moc_s(k,1) = moc_s(k,1) + WORK1_G(i,j,k)
          enddo ! k
        endif
      enddo
      moc_s(km,:) = - dz(km) * moc_s(km,:)
      do k = km-1, 1, -1
        moc_s(k,:) = moc_s(k+1,:) - dz(k) * moc_s(k,:)
      enddo

      ! Add sounthern boundary to MOC
      do k = 1, km
        do n = 1, n_lat_aux_grid+1
          TAVG_MOC_G(n,k,1) = TAVG_MOC_G(n,k,1) + moc_s(k,1)
        enddo
      enddo

     !-----------------------------------------------------------------------
     !
     !  convert MOC to Sverdrups, prior to masking
     !
     !-----------------------------------------------------------------------
      TAVG_MOC_G = TAVG_MOC_G*1.0e-12_r8

      amoc_strength = TAVG_MOC_G(amoc_strength_target_lat_idx,amoc_strength_target_depth_idx,1)

      Write(*,*) "  "
      Write(*,*) "==============================="
      Write(*,*) "Strength idx :", amoc_strength_target_lat_idx, ", ", amoc_strength_target_depth_idx
      Write(*,*) "AMOC strength : ", amoc_strength
      Write(*,*) "==============================="
      Write(*,*) "  "

    endif ! master_task

    deallocate ( WORK1, WORK2, WORK1_G )

!-----------------------------------------------------------------------
!EOC
  end subroutine compute_amoc

end module tavg
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
