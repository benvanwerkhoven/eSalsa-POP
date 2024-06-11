!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module io_netcdf_stoch

!BOP
! !MODULE: io_netcdf
! !DESCRIPTION:
!  This module provides a generic input/output interface
!  for writing arrays in netCDF format.
!
! !REVISION HISTORY:
!  SVN:$Id: io_netcdf.F90 17212 2009-07-20 23:01:42Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod

   use kinds_mod
   use domain_size
   use domain
   use blocks
   use constants
   use communicate
   use broadcast
   use gather_scatter
   use exit_mod
   use io_types
   use io_tools
   use netcdf

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: check_status, &
             read_AR_netcdf_file_header, &
             read_AR_netcdf_file_data, &
             read_EOF_netcdf_file

!EOP
!BOC


!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE:  check_status
! !INTERFACE:

 subroutine check_status(status)

! !DESCRIPTION:
!  This exception handler subroutine can be used to check error status
!  after a netcdf call.  It prints out a text message assigned to
!  an error code but does not exit because this routine is typically
!  only called from a single process.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (i4), intent (in) ::  &
      status                     ! status returned by netCDF call

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  call netCDF routine to return error message
!
!-----------------------------------------------------------------------

   if (status /= nf90_noerr) then
      write(stdout,*) trim(nf90_strerror(status))
      call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine check_status

 subroutine read_AR_netcdf_file_header(data_file, ext_name, eof_l, lag_l)

   implicit none

   type (datafile), intent (inout)  :: data_file
   character (*) :: &
      ext_name        ! The name of the group to load data from

   integer (i4), intent(out) ::  &
     eof_l, &
     lag_l

   ! local variables
   character (char_len) :: &
      path            ! filename to read

   character (80) :: &
      work_line,     &! temporary to use for parsing file lines
      att_name        ! temporary to use for attribute names

   integer (i4) ::  &
      iostat,       &! status flag
      ncid,         &! netCDF file id
      dimid,        &! netCDF dimension id
      ndim, nvar,natt,k_un, k_form, k, d,&
      nsize,        &! size parameter returned by inquire function
      n,            &! loop index
      itype,        &! netCDF data type
      att_ival,     &! netCDF data type
      num_atts       ! number of global attributes

    integer (i4), dimension(20) :: &
      ncids          ! netCDF group ids

   character (20) :: &
      dim_eof, dim_lag


   logical (log_kind) :: &
      att_lval           ! temp space for logical attribute

   real (r4) ::     &
      att_rval       ! temp space for real attribute

   real (r8) ::     &
      att_dval       ! temp space for double attribute

   logical (log_kind) :: &
      attrib_error        ! error flag for reading attributes

!-----------------------------------------------------------------------
!
!  set the readonly flag in the data file descriptor
!
!-----------------------------------------------------------------------

   data_file%readonly = .true.

!-----------------------------------------------------------------------
!
!  open the netCDF file
!
!-----------------------------------------------------------------------

   iostat = nf90_noerr
   data_file%id = 0

   if (my_task == master_task) then
      path = trim(data_file%full_name)
      iostat = nf90_open(path=trim(path), mode=nf90_nowrite, ncid=ncid)
      call check_status(iostat)
   endif

   call broadcast_scalar(iostat, master_task)
   if (iostat /= nf90_noerr) then
      write(stdout,*) 'filename = ', trim(data_file%full_name)
      call exit_POP(sigAbort,'error opening netCDF file for reading')
   endif

   call broadcast_scalar(ncid, master_task)
   data_file%id(1) = ncid

   if (my_task == master_task) then
     write(dim_eof,"(2A)") "eof_d_",trim(ext_name)
     write(dim_lag,"(2A)") "lag_d_",trim(ext_name)
     iostat = nf90_inquire(ncid, nDimensions=ndim, nVariables=nvar, nAttributes=natt)
     do k = 1, ndim
       iostat = nf90_inquire_dimension(ncid, k, att_name, att_ival)
       if (att_name == trim(dim_eof)) then
         eof_l = att_ival
       else if (att_name == trim(dim_lag)) then
         lag_l = att_ival
       endif
     enddo
   endif

   call broadcast_scalar(iostat, master_task)
   if (iostat /= nf90_noerr) &
      call exit_POP(sigAbort, &
                    'error getting dimensions in group')

   call broadcast_scalar(eof_l, master_task)
   call broadcast_scalar(lag_l, master_task)

   if (my_task == master_task) then
      iostat = nf90_close(ncid)
   end if

 end subroutine read_AR_netcdf_file_header

 subroutine read_AR_netcdf_file_data(data_file, ext_name,   &
                                     eof_l, lag_l,          &
                                     lags, sig, rho, hist, rnd)

   implicit none

   type (datafile), intent (inout)  :: data_file
   character (*) :: &
      ext_name        ! The name of the group to load data from

   integer (i4), intent(in) ::  &
     eof_l, &
     lag_l

   integer(int_kind), dimension(eof_l), intent(out) :: lags
   real(r8) , dimension(eof_l), intent(out) :: sig, rnd
   real(r8) , dimension(lag_l, eof_l), intent(out) :: rho, hist

   ! local variables
   character (char_len) :: &
      path            ! filename to read

   character (80) :: &
      work_line,     &! temporary to use for parsing file lines
      att_name        ! temporary to use for attribute names

   integer (i4) ::  &
      iostat,       &! status flag
      ncid,         &! netCDF file id
      dimid,        &! netCDF dimension id
      ndim, nvar,natt,k_un, k_form, k, xtype, &
      nsize,        &! size parameter returned by inquire function
      n,            &! loop index
      itype,        &! netCDF data type
      att_ival,     &! netCDF data type
      num_atts       ! number of global attributes

    integer (i4), allocatable, dimension(:) :: dimids

   character (20) :: &
      v_lag, v_rho, v_rnd, v_sig, v_hist

   integer (i4), dimension(20) :: &
      ncids          ! netCDF group ids

   logical (log_kind) :: &
      att_lval           ! temp space for logical attribute

   real (r4) ::     &
      att_rval       ! temp space for real attribute

   real (r8) ::     &
      att_dval       ! temp space for double attribute

   logical (log_kind) :: &
      attrib_error        ! error flag for reading attributes

!-----------------------------------------------------------------------
!
!  set the readonly flag in the data file descriptor
!
!-----------------------------------------------------------------------

   data_file%readonly = .true.

!-----------------------------------------------------------------------
!
!  open the netCDF file
!
!-----------------------------------------------------------------------

   iostat = nf90_noerr
   data_file%id = 0

   if (my_task == master_task) then
      path = trim(data_file%full_name)
      iostat = nf90_open(path=trim(path), mode=nf90_nowrite, ncid=ncid)
      call check_status(iostat)
   endif

   call broadcast_scalar(iostat, master_task)
   if (iostat /= nf90_noerr) then
      write(stdout,*) 'filename = ', trim(data_file%full_name)
      call exit_POP(sigAbort,'error opening netCDF file for reading')
   endif

   call broadcast_scalar(ncid, master_task)
   data_file%id(1) = ncid

   if (my_task == master_task) then
     write(v_rho,"(2A)") "rho_",trim(ext_name)
     write(v_sig,"(2A)") "sig_",trim(ext_name)
     write(v_rnd,"(2A)") "rnd_",trim(ext_name)
     write(v_hist,"(2A)") "hist_",trim(ext_name)
     write(v_lag,"(2A)") "lags_",trim(ext_name)

     iostat = nf90_inquire(ncid, nDimensions=ndim, nVariables=nvar, nAttributes=natt)
     if (iostat /= nf90_noerr) then
       write(*,*) "Error getting variables information"
     endif
     do k = 1, ndim
        iostat = nf90_inquire_dimension(ncid, k, att_name, att_ival )
     enddo
     allocate(dimids(1:ndim))
     do k = 1, nvar
       iostat = nf90_inquire_variable(ncid, k, att_name, xtype, ndim, dimids, natt )
       if (iostat /= nf90_noerr) then
         write(*,*) "Error getting information for variable ", trim(att_name)
       endif
       if (trim(att_name) == v_lag) then
         iostat = nf90_get_var(ncid, k, lags)
       else if (trim(att_name) == v_sig) then
         iostat = nf90_get_var(ncid, k, sig)
       else if (trim(att_name) == v_rnd) then
         iostat = nf90_get_var(ncid, k, rnd)
       else if (trim(att_name) == v_rho) then
         iostat = nf90_get_var(ncid, k, rho, (/ 1,1 /))
       else if (trim(att_name) == v_hist) then
         iostat = nf90_get_var(ncid, k, hist, (/ 1,1 /))
       endif
       if (iostat /= nf90_noerr) then
         write(*,*) "Error reading variable data ", trim(att_name), NF90_STRERROR(iostat)
       endif
     enddo
   endif

   call broadcast_scalar(iostat, master_task)
   if (iostat /= nf90_noerr) &
      call exit_POP(sigAbort, &
                    'error getting AR data file')

   call broadcast_array(lags, master_task)
   call broadcast_array(sig,  master_task)
   call broadcast_array(rnd,  master_task)
   call broadcast_array(rho,  master_task)
   call broadcast_array(hist, master_task)

   if (my_task == master_task) then
      iostat = nf90_close(ncid)
   end if

 end subroutine read_AR_netcdf_file_data

 subroutine read_EOF_netcdf_file(data_file, &
                                 eof_l, eof)

   implicit none

   type (datafile), intent (inout)  :: data_file
   integer (i4), intent(in) ::  &
     eof_l

   real(r8), dimension(nx_block,ny_block,max_blocks_clinic,eof_l), intent(out) :: eof

   ! local variables
   real(r8), dimension(:,:), allocatable :: eof_gbl
   real(r8), dimension(:,:), allocatable :: eof_gbl_flipped

   character (char_len) :: &
      path            ! filename to read

   character (80) :: &
      work_line,     &! temporary to use for parsing file lines
      att_name        ! temporary to use for attribute names

   integer (i4) ::  &
      iostat,       &! status flag
      ncid,         &! netCDF file id
      dimid,        &! netCDF dimension id
      ndim, nvar,natt,k_un, k_form, k, xtype, e, i,j, &
      nsize,        &! size parameter returned by inquire function
      n,            &! loop index
      itype,        &! netCDF data type
      att_ival,     &! netCDF data type
      num_atts       ! number of global attributes

    integer (i4), allocatable, dimension(:) :: dimids

    integer (i4), dimension(20) :: &
      ncids          ! netCDF group ids


   logical (log_kind) :: &
      att_lval           ! temp space for logical attribute

   real (r4) ::     &
      att_rval       ! temp space for real attribute

   real (r8) ::     &
      att_dval       ! temp space for double attribute

   logical (log_kind) :: &
      attrib_error        ! error flag for reading attributes

   if (my_task == master_task) then
     allocate(eof_gbl(nx_global,ny_global))
     allocate(eof_gbl_flipped(ny_global,nx_global))
   endif

!-----------------------------------------------------------------------
!
!  set the readonly flag in the data file descriptor
!
!-----------------------------------------------------------------------

   data_file%readonly = .true.

!-----------------------------------------------------------------------
!
!  open the netCDF file
!
!-----------------------------------------------------------------------

   iostat = nf90_noerr
   data_file%id = 0

   if (my_task == master_task) then
      path = trim(data_file%full_name)
      iostat = nf90_open(path=trim(path), mode=nf90_nowrite, ncid=ncid)
      call check_status(iostat)
   endif

   call broadcast_scalar(iostat, master_task)
   if (iostat /= nf90_noerr) then
      write(stdout,*) 'filename = ', trim(data_file%full_name)
      call exit_POP(sigAbort,'error opening netCDF file for reading')
   endif

   call broadcast_scalar(ncid, master_task)
   data_file%id(1) = ncid

   ! Get and check dimension data
   if (my_task == master_task) then
     iostat = nf90_inquire(ncid, nDimensions=ndim, nVariables=nvar, nAttributes=natt)
     do k = 1, ndim
       iostat = nf90_inquire_dimension(ncid, k, att_name, att_ival )
       if (trim(att_name) == "eof_d") then
         if (att_ival /= eof_l) then
           write(*,*) "Error EOF data length ", att_ival, " do not match AR one ", eof_l
         endif
       else if (trim(att_name) == "lon_d") then
         if (att_ival /= nx_global) then
           write(*,*) "Error EOF data lon length ", att_ival, " do not match sim. one ", nx_global
         endif
       else if (trim(att_name) == "lat_d") then
         if (att_ival /= ny_global) then
           write(*,*) "Error EOF data lat length ", att_ival, " do not match sim. one ", ny_global
         endif
       endif
     enddo

     allocate(dimids(1:ndim))
   endif

   call broadcast_scalar(nvar, master_task)

   do k = 1, nvar
     if (my_task == master_task) then
       iostat = nf90_inquire_variable(ncid, k, att_name, xtype, ndim, dimids, natt )
     endif
     do e = 1, eof_l
       if (my_task == master_task) then
         iostat = nf90_get_var(ncid, k, eof_gbl_flipped, start=(/ 1,1,e /), count=(/ny_global,nx_global,1/))
         do i = 1, ny_global
           do j = 1, nx_global
             eof_gbl(j,i) = eof_gbl_flipped(i,j)
           enddo
         enddo
       endif
       call broadcast_scalar(iostat, master_task)

       if (iostat /= nf90_noerr) then
         if (my_task == master_task) then
           write(*,*) "Error while loading EOF data "
         endif
         exit
       else
         call scatter_global(eof(:,:,:,e), &
                             eof_gbl, master_task, distrb_clinic, &
                             field_loc_center, field_type_scalar)
       endif
     enddo
   enddo

   call broadcast_scalar(iostat, master_task)
   if (iostat /= nf90_noerr) &
      call exit_POP(sigAbort, &
                    'error while reading EOF data file')

   if (my_task == master_task) then
     deallocate(eof_gbl,eof_gbl_flipped)
   endif

 end subroutine read_EOF_netcdf_file

!***********************************************************************
 end module io_netcdf_stoch

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
