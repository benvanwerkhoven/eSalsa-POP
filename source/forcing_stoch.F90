module forcing_stoch

!BOP
! !MODULE: forcing_stoch
!
! !DESCRIPTION:
!  This module contains the surface stochastic forcing
!  data and methods.
!
! !REVISION HISTORY:
!
! !USES:

  use kinds_mod
  use blocks
  use distribution
  use POP_CommMod
  use domain
  use constants
  use prognostic
  use io
  use grid
  use io_netcdf_stoch
  use time_management
  use global_reductions, only:global_sum
  use exit_mod

  implicit none

  real (r8) :: sf_fwf_fraction = 0.0

  ! Data for ERA5 forcing
  ! Dimensions
  integer (int_kind) :: &
    n_eof_ep,                   &! Number of EOF for E-P forcing
    n_eof_t2m,                  &! Number of EOF for T forcing
    lag_max_ep,                 &! Largest lag for E-P forcing
    lag_max_t2m,                &! Largest lag for T forcing
    nlat_eof,                   &! Zonal length of EOFs
    nlon_eof                     ! Meridional length of EOFs

  ! Autoregression (AR) data
  integer (int_kind), allocatable, dimension(:) :: &
    lags_ep,                    &! Lags for E-P forcing (n_eof)
    lags_t2m                     ! Lags for T forcing
  real (r8), allocatable, dimension(:) :: &
    sig_ep,                     &! White noise ampl for E-P (n_eof)
    sig_t2m                      ! White noise ampl for T
  real (r8), allocatable, dimension(:,:) :: &
    rho_ep,                     &! Decay coefficient for E-P (n_eof * l_max)
    rho_t2m                      ! Decay coefficient for T

  ! History data
  ! AR requires historical data up to lag(i) for i in n_eof
  real (r8), allocatable, dimension(:,:) :: &
    hist_ep,                    &! Forcing history for E-P (n_eof * l_max)
    hist_t2m                     ! Forcing history for T

  ! EOFs data
  ! EOFs are on a cartesian grids centered on the Atlantic
  real (r8), allocatable, dimension(:,:,:,:) :: &
    eof_ep,                     &! EOFs of E-P from PC analysis of ERA5
    eof_t2m                      ! EOFs of T from PC analysis of ERA5

  ! Monthly forcing data
  ! Keep around 2 randomly generated monthly forcing fields
  ! to interpolate from linearly
  ! "old": beginning of the current month
  ! "new": beginning of next month
  real (r8), allocatable, dimension(:,:,:,:) :: &
    stoch_data_ep,             &
    stoch_data_t2m

  ! Timing data for updating stochastic fields
  real (r8), dimension(1:12) :: &
    sf_fwf_time,                &
    sf_hf_time
  character (char_len) ::       &
    sf_fwf_interp_type,         &
    sf_hf_interp_type
  character (char_len) ::       &
    sf_fwf_data_type,           &
    sf_hf_data_type
  real (r8) ::                  &
    sf_fwf_data_update,         &
    sf_hf_data_update
  real (r8) ::                  &
    sf_fwf_forcing_update,      &
    sf_hf_forcing_update
  integer (int_kind) ::         &
    sf_fwf_data_update_idx,     &
    sf_hf_data_update_idx

  real(r8), allocatable, dimension(:) :: &
    rnd_ep,                     &! Latest set of random number for E-P
    rnd_t2m                      ! Latest set of random number for T

  logical :: stochastic_forcing_fwf
  logical :: stochastic_forcing_hf

  character (char_len) :: &
    sf_fwf_formulation,         &! Formulation of the freshwater stoch forcing
    sf_hf_formulation            ! Formulation of the heat stoch forcing

  character (char_len) :: &
    ARdata_ep_filename,         &! Autoregression datafile for E-P
    ARdata_t_filename,          &! Autoregression datafile for temperature
    EOF_ep_filename,            &! EOF datafile for E-P
    EOF_t_filename               ! EOF datafile for temperature

contains

  subroutine read_stoch_forcing_namelist
    implicit none
    integer(int_kind) :: nml_error ! namelist error flag
    namelist /forcing_stoch/ stochastic_forcing_fwf, stochastic_forcing_hf,  &
                              sf_fwf_formulation, sf_hf_formulation,            &
                              ARdata_ep_filename, ARdata_t_filename,            &
                              EOF_ep_filename, EOF_t_filename,                  &
                              sf_fwf_fraction

    ! Set some defaults
    stochastic_forcing_fwf   = .false.
    stochastic_forcing_hf    = .false.
    sf_fwf_formulation        = 'baseline-frac'
    sf_hf_formulation         = 'baseline-frac'
    ARdata_ep_filename        = 'unknown'
    ARdata_t_filename         = 'unknown'
    EOF_ep_filename           = 'unknown'
    EOF_t_filename            = 'unknown'
    sf_fwf_fraction           = 0.0

    if (my_task == master_task) then
       open (nml_in, file=nml_filename, status='old',iostat=nml_error)
       if (nml_error /= 0) then
          nml_error = -1
       else
          nml_error =  1
       endif
       do while (nml_error > 0)
          read(nml_in, nml=forcing_stoch,iostat=nml_error)
       end do
       if (nml_error == 0) close(nml_in)
    endif

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading forcing_stoch')
    endif

    call broadcast_scalar(stochastic_forcing_fwf, master_task)
    call broadcast_scalar(stochastic_forcing_hf,  master_task)
    call broadcast_scalar(sf_fwf_formulation,     master_task)
    call broadcast_scalar(sf_hf_formulation,      master_task)
    call broadcast_scalar(ARdata_ep_filename,     master_task)
    call broadcast_scalar(ARdata_t_filename,      master_task)
    call broadcast_scalar(EOF_ep_filename,        master_task)
    call broadcast_scalar(EOF_t_filename,         master_task)
    call broadcast_scalar(sf_fwf_fraction,        master_task)

  end subroutine read_stoch_forcing_namelist

  subroutine init_stoch_forcing(STF_stoch)

    implicit none

    ! !INOUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(inout) :: &
       STF_stoch! surface tracer fluxes for all tracers

    ! LOCAL
    real(r8), dimension(:), allocatable :: rnd_u1, rnd_u2
    integer (int_kind) :: seed_l
    integer (int_kind), dimension(:), allocatable:: rnd_seed
    type (datafile) :: ARfile, EOFfile

    ! Set forcing to zero
    STF_stoch(:,:,:,:) = c0

    ! Return if stoch sfwf nor shf activated
    if (.not. stochastic_forcing_fwf .and. &
        .not. stochastic_forcing_hf ) then
      return
    endif

    if (stochastic_forcing_fwf) then
      select case (sf_fwf_formulation)
      case ('ERA5-Data')
        if (my_task == master_task) then
          write(*,*) "Reading AR E-P data from ", ARdata_ep_filename
        endif

        ARfile = construct_file("nc", &
                                full_name=trim(ARdata_ep_filename))
        call read_AR_netcdf_file_header(ARfile,"ep",n_eof_ep,lag_max_ep)

        allocate(lags_ep(1:n_eof_ep))
        allocate(sig_ep(1:n_eof_ep))
        allocate(rho_ep(1:lag_max_ep,1:n_eof_ep))
        allocate(hist_ep(1:lag_max_ep,1:n_eof_ep))

        allocate(rnd_ep(1:n_eof_ep))

        call read_AR_netcdf_file_data(ARfile,"ep",n_eof_ep,lag_max_ep,&
                                      lags_ep, sig_ep, rho_ep, hist_ep, rnd_ep)

        allocate(eof_ep(nx_block,ny_block,max_blocks_clinic,n_eof_ep))

        EOFfile = construct_file("nc", &
                                full_name=trim(EOF_ep_filename))

        call read_EOF_netcdf_file(EOFfile, n_eof_ep, eof_ep)

        allocate(stoch_data_ep(nx_block,ny_block,max_blocks_clinic,2))
        stoch_data_ep(:,:,:,:) = 0.0

        ! Find time to update the stochastic fields
        ! Hard-coded types for now.
        sf_fwf_time(1:12) = 0.0                 !
        sf_fwf_interp_type = "linear"           ! Interpolation type (between 2 EOFs-based fields)
        sf_fwf_data_type = "monthly-calendar"   ! Working with monthly data here
        sf_fwf_data_update = 0.0                ! Next time a new EOFs-based field is needed
        sf_fwf_forcing_update = 0.0             ! Next time forcing need update

        ! Light version of the content of the find_forcing_times subroutine
        call find_stoch_forcing_times(sf_fwf_time, sf_fwf_forcing_update, sf_fwf_data_update, &
                                      sf_fwf_data_update_idx)

        ! Generate a forcing for current month (old)
        call generate_new_forcing(n_eof_ep, lag_max_ep, 1, rnd_ep, hist_ep, &
                                  lags_ep, rho_ep, sig_ep, eof_ep, stoch_data_ep)

        ! Generate a forcing for next month (new), using new set of random numbers
        call generate_new_forcing(n_eof_ep, lag_max_ep, 0, rnd_ep, hist_ep, &
                                  lags_ep, rho_ep, sig_ep, eof_ep, stoch_data_ep)
      end select
    endif

    if (stochastic_forcing_hf) then
      select case (sf_hf_formulation)
      case ('ERA5-Data')
        if (my_task == master_task) then
          write(*,*) "Reading AR t2m data from ", ARdata_t_filename
        endif

        ARfile = construct_file("nc", &
                                full_name=trim(ARdata_t_filename))
        call read_AR_netcdf_file_header(ARfile,"t2m",n_eof_t2m,lag_max_t2m)

        if (my_task == master_task) then
          write(*,*) n_eof_t2m, lag_max_t2m
        endif

        allocate(lags_t2m(1:n_eof_t2m))
        allocate(sig_t2m(1:n_eof_t2m))
        allocate(rho_t2m(1:lag_max_t2m,1:n_eof_t2m))
        allocate(hist_t2m(1:lag_max_t2m,1:n_eof_t2m))
        allocate(rnd_t2m(1:n_eof_t2m))

        call read_AR_netcdf_file_data(ARfile,"t2m",n_eof_t2m,lag_max_t2m,&
                                      lags_t2m, sig_t2m, rho_t2m, hist_t2m, rnd_t2m)

        allocate(eof_t2m(nx_block,ny_block,max_blocks_clinic,n_eof_t2m))

        EOFfile = construct_file("nc", &
                                full_name=trim(EOF_t_filename))

        call read_EOF_netcdf_file(EOFfile, n_eof_t2m, eof_t2m)

        allocate(stoch_data_t2m(nx_block,ny_block,max_blocks_clinic,2))
        stoch_data_t2m(:,:,:,:) = 0.0

        ! Find time to update the stochastic fields
        ! Hard-coded types for now.
        sf_hf_time(1:12) = 0.0                 !
        sf_hf_interp_type = "linear"           ! Interpolation type (between 2 EOFs-based fields)
        sf_hf_data_type = "monthly-calendar"   ! Working with monthly data here
        sf_hf_data_update = 0.0                ! Next time a new EOFs-based field is needed
        sf_hf_forcing_update = 0.0             ! Next time forcing need update

        ! Light version of the content of the find_forcing_times subroutine
        call find_stoch_forcing_times(sf_hf_time, sf_hf_forcing_update, sf_hf_data_update, &
                                      sf_hf_data_update_idx)

        ! Generate a forcing for current month (old)
        call generate_new_forcing(n_eof_t2m, lag_max_t2m, 1, rnd_t2m, hist_t2m, &
                                  lags_t2m, rho_t2m, sig_t2m, eof_t2m, stoch_data_t2m)

        ! Generate a forcing for next month (new), using new set of random numbers
        call generate_new_forcing(n_eof_t2m, lag_max_t2m, 0, rnd_t2m, hist_t2m, &
                                  lags_t2m, rho_ep, sig_t2m, eof_t2m, stoch_data_t2m)
      end select
    endif

  end subroutine init_stoch_forcing

  subroutine append_stoch_forcing_sfwf(STF)

    use forcing_fields, only : STF_stoch

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(inout) :: &
       STF     ! surface tracer fluxes for all tracers

    ! LOCAL
    integer (int_kind) :: &
       iblock              ! block loop index

    ! Return if sfwf not activated
    if (.not. stochastic_forcing_fwf) then
      return
    endif

    ! Compute the stoicastic freshwater forcing
    call calc_stochastic_sfwf(STF, STF_stoch)

    ! Add stochastic forcing to baseline forcing
    !$OMP PARALLEL DO PRIVATE(iblock)
    do iblock = 1, nblocks_clinic
       STF(:,:,2,iblock) = STF(:,:,2,iblock) + STF_stoch(:,:,2,iblock)
    enddo
    !$OMP END PARALLEL DO

  end subroutine append_stoch_forcing_sfwf

  subroutine append_stoch_forcing_shf(STF)

    use forcing_fields, only : STF_stoch

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(inout) :: &
       STF     ! surface tracer fluxes for all tracers

    ! LOCAL
    integer (int_kind) :: &
       iblock              ! block loop index

    ! Return if stoch shf not activated
    if (.not. stochastic_forcing_hf) then
      return
    endif

    ! The ERA5-Data return a stochastic temperature
    ! while baseline-frac, return a flux -> exit with the former
    if (sf_fwf_formulation == 'ERA5-Data') then
      return
    endif

   ! Compute the stoicastic heat forcing
   call calc_stochastic_shf(STF, STF_stoch)

   ! Add stochastic forcing to baseline forcing
   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
      STF(:,:,1,iblock) = STF(:,:,1,iblock) + STF_stoch(:,:,1,iblock)
   enddo
   !$OMP END PARALLEL DO

  end subroutine append_stoch_forcing_shf

  subroutine append_stoch_temp_shf(STF, AST)

    use forcing_fields, only : STF_stoch

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
       intent(in) :: &
       STF     ! Air temperature
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
       intent(inout) :: &
       AST     ! Air temperature

    ! LOCAL
    integer (int_kind) :: &
       iblock              ! block loop index

    ! Return if stoch shf not activated
    if (.not. stochastic_forcing_hf) then
      return
    endif

    ! The ERA5-Data return a stochastic temperature
    ! while baseline-frac, return a flux -> exit with the latest
    if (sf_fwf_formulation == 'baseline-frac') then
      return
    endif

    ! Compute the stocastic heat forcing
    call calc_stochastic_shf(STF, STF_stoch)

    ! Add stochastic air temperature to baseline value
    !$OMP PARALLEL DO PRIVATE(iblock)
    do iblock = 1, nblocks_clinic
       AST(:,:,iblock) = AST(:,:,iblock) + STF_stoch(:,:,1,iblock)
    enddo
    !$OMP END PARALLEL DO

  end subroutine append_stoch_temp_shf

  subroutine calc_stochastic_sfwf(STF, STF_stoch)
    ! !DESCRIPTION:
    ! Used to add stochastic forcing in the northern atlantic ocean
    ! as a fraction of the baseline forcing value in that region. The
    ! net addition of forcing should be zero, such that the mean forcing
    ! is removed.
    !
    !  Notes:

    implicit none

    ! !INPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(in) :: &
       STF     ! surface tracer fluxes for all tracers

    ! !OUTPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(out) :: &
       STF_stoch! surface tracer fluxes for all tracers

    ! LOCAL
    type(block) :: &
       this_block  ! block info for current block

    integer (int_kind) :: &
       ierr,&              ! MPI calls error flag
       i ,j,&              ! cell loop indices
       iblock              ! block loop index

    real (r8) :: &
       stf_int_lcl,&       ! Integral of FWF (local)
       stf_int,&           ! Integral of FWF (global)
       surf_lcl,&          ! Covered surface (local)
       surf                ! Covered surface (global)

    real (r8) :: pold, pnew! Interpolation coefficients for ERA5-Data forcing

    select case (sf_fwf_formulation)
    case ('baseline-frac')
      ! Stochastic forcing is based on a fraction of the
      ! baseline forcing computed previously.
      ! 1. Mask STF on northern atlantic
      ! 2. Compute mean of masked STF
      ! 3. STF_stoch = sf_fwf_fraction * (STF_m - \int(STF_m))
      ! 4. STF += STF_stoch

      ! Compute the integral of the STF term and surface
      ! of the controlled area (area converted from cm^2 to km^2).
      stf_int_lcl = 0.0
      surf_lcl = 0.0

      !$OMP PARALLEL DO PRIVATE(iblock, this_block)
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               if (KMT(i,j,iblock) > 0 .and. &
                  TLAT(i,j,iblock) >= 0.26179939 .and. &
                  TLAT(i,j,iblock) <= 1.22173048 .and. &
                  TLON(i,j,iblock) >= 4.53785606)  then
                  STF_stoch(i,j,2,iblock) = STF(i,j,2,iblock)
                  stf_int_lcl = stf_int_lcl + TAREA(i,j,iblock) * 1.0e-10_r8 * STF_stoch(i,j,2,iblock)
                  surf_lcl = surf_lcl + TAREA(i,j,iblock) * 1.0e-10_r8
               else
                  STF_stoch(i,j,2,iblock) = c0
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      stf_int = global_sum(stf_int_lcl, distrb_clinic)
      !write(*,*) " >> Integral of stoch SFWF before correction: ", stf_int
      surf = global_sum(surf_lcl, distrb_clinic)
      !write(*,*) " >> Surface covered by stoch SFWF: ", surf
      stf_int = stf_int / surf
      !write(*,*) " >> Mean stoch SFWF: ", stf_int

      stf_int_lcl = 0.0
      !$OMP PARALLEL DO PRIVATE(iblock, this_block)
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               if (KMT(i,j,iblock) > 0 .and. &
                  TLAT(i,j,iblock) >= 0.26179939 .and. &
                  TLAT(i,j,iblock) <= 1.22173048 .and. &
                  TLON(i,j,iblock) >= 4.53785606)  then
                  STF_stoch(i,j,2,iblock) = sf_fwf_fraction * (STF_stoch(i,j,2,iblock) - stf_int)
                  stf_int_lcl = stf_int_lcl + TAREA(i,j,iblock) * 1.0e-10_r8 * STF_stoch(i,j,2,iblock)
               else
                  STF_stoch(i,j,2,iblock) = c0
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      stf_int = global_sum(stf_int_lcl, distrb_clinic)

    case ('ERA5-Data')
      ! Check if a new random field has to be generated
      if (thour00 >= sf_fwf_data_update) then
        ! Generate a forcing field from current data
        call generate_new_forcing(n_eof_ep, lag_max_ep, 0, rnd_ep, hist_ep, &
                                  lags_ep, rho_ep, sig_ep, eof_ep, stoch_data_ep)

        ! Udpate the time management
        sf_fwf_data_update_idx = mod(imonth+1,12)
        if (sf_fwf_data_update_idx == 0) sf_fwf_data_update_idx = 12

        if (imonth == 12 .AND. sf_fwf_time(1) < thour00) then
          sf_fwf_time(:) = sf_fwf_time(:) + hours_in_year
        endif

        sf_fwf_data_update = sf_fwf_time(sf_fwf_data_update_idx)
      endif

      ! Linear interpolation between an 'old' and 'new' monthly random fields.
      call get_linear_coeff(sf_fwf_time, pold, pnew)

      ! Conversion: stoch data from EOFs in m/s. -> convert to kg/s/m^2
      ! then to salf flux (CGS) msu*cm/s
      ! using water density rho_fw
      do iblock = 1, nblocks_clinic
        STF_stoch(:,:,2,iblock) =  (pold * stoch_data_ep(:,:,iblock,2) + &
                                    pnew * stoch_data_ep(:,:,iblock,1)) * &
                                    1.0e3_r8 / rho_fw  * & ! kg/s/m^2
                                    salinity_factor        ! msu*cm/s
      end do
    end select

  end subroutine calc_stochastic_sfwf

  subroutine calc_stochastic_shf(STF, STF_stoch)
    ! !DESCRIPTION:
    ! Used to add stochastic forcing in the northern atlantic ocean
    ! as a fraction of the baseline forcing value in that region. The
    ! net addition of forcing should be zero, such that the mean forcing
    ! is removed.
    !
    !  Notes:

    implicit none

    ! !INPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(in) :: &
       STF        ! surface tracer fluxes for all tracers

    ! !OUTPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(out) :: &
       STF_stoch ! surface tracer fluxes for all tracers

    ! LOCAL
    integer (int_kind) :: &
       iblock              ! block loop index

    real (r8) :: pold, pnew! Interpolation coefficients for ERA5-Data forcing

    select case (sf_fwf_formulation)
    case ('baseline-frac')      ! TODO
      do iblock = 1, nblocks_clinic
        STF_stoch(:,:,1,iblock) = c0
      end do
    case ('ERA5-Data')
      ! Check if a new random field has to be generated
      if (thour00 >= sf_hf_data_update) then
        ! Generate a forcing field from current data
        call generate_new_forcing(n_eof_t2m, lag_max_t2m, 0, rnd_t2m, hist_t2m, &
                                  lags_t2m, rho_t2m, sig_t2m, eof_t2m, stoch_data_t2m)

        ! Udpate the time management
        sf_hf_data_update_idx = mod(imonth+1,12)
        if (sf_hf_data_update_idx == 0) sf_hf_data_update_idx = 12

        if (imonth == 12 .AND. sf_hf_time(1) < thour00) then
          sf_hf_time(:) = sf_hf_time(:) + hours_in_year
        endif

        sf_hf_data_update = sf_hf_time(sf_hf_data_update_idx)
      endif

      ! Linear interpolation between an 'old' and 'new' monthly random fields.
      call get_linear_coeff(sf_hf_time, pold, pnew)

      ! Conversion: stoch data from EOFs in K to flux.
      do iblock = 1, nblocks_clinic
        STF_stoch(:,:,1,iblock) =  (pold * stoch_data_t2m(:,:,iblock,2) + &
                                    pnew * stoch_data_t2m(:,:,iblock,1))
      end do
    end select

  end subroutine calc_stochastic_shf

  subroutine find_stoch_forcing_times(forcing_time, sf_forcing_time_next, sf_data_update, &
                                      sf_data_update_idx)

    ! Always assume stochastic forcing work with calendar month + linear interpolation

    implicit none

    ! INOUT
    real (r8), dimension(1:12), intent(inout) :: &
      forcing_time

    ! OUT
    real (r8), intent(out) :: &
      sf_forcing_time_next, sf_data_update
    integer (int_kind), intent(out) :: &
      sf_data_update_idx

    ! LOCAL
    integer (int_kind) :: &
      n, month_minloc, second

    forcing_time(:) = thour00_begin_this_year + thour00_begmonth_calendar(:)
    do n = 1, 12
      if (forcing_time(n) > thour00) exit
    enddo

    month_minloc = n - 1
    month_minloc = mod(month_minloc,12)
    if (month_minloc <= 0 ) then
      month_minloc = month_minloc + 12
    endif

    second = mod(month_minloc+1,12)
    if (second == 0) second = 12

    sf_forcing_time_next = forcing_time(second+1)

    sf_data_update = forcing_time(second)

  end subroutine find_stoch_forcing_times

  subroutine generate_new_forcing(eof_l, lag_l, offset, &
                                  new_rnd_set, rnd_hist, lags, &
                                  rho, sigma, eofs, stoch_fields)

    implicit none

    ! IN
    integer (int_kind), intent(in) :: &
      eof_l,                          & ! Number of EOFs to assemble
      lag_l,                          & ! Lag of each EOF
      offset                            ! Offset

    real (r8), dimension(eof_l), intent(in) :: &
      new_rnd_set                       ! New set of random number (N(0,1)) to use and add to history

    integer (int_kind), dimension(eof_l), intent(in) :: &
      lags                              ! Lag of each EOF
    real (r8), dimension(lag_l, eof_l), intent(in) :: &
      rho                               ! Weight of each EOFs for each lag time in AR
    real (r8), dimension(eof_l), intent(in) :: &
      sigma                             ! Scaling factor applied on random process for each EOFs
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic,eof_l), intent(in) :: &
      eofs                              ! EOFs obtained from PCA, ...

    ! INOUT
    real (r8), dimension(lag_l, eof_l), intent(inout) :: &
      rnd_hist                          ! Random number history, included in AR
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic,2), intent(inout) :: &
      stoch_fields                      ! Final stochastic fields

    ! LOCAL
    integer (int_kind) :: &
      i_eof, i_lag

    ! stoch_fields contains an 'old (2)' and a 'new (1)' fields.
    ! Move new into old
    stoch_fields(:,:,:,2) = stoch_fields(:,:,:,1)

    ! Reinit new field
    stoch_fields(:,:,:,1) = 0.0

    ! If we have an offset, no need to regenerate a set of data,
    ! we are using an older set in the history
    if (offset == 0) then
      ! Add to history using AR and the new set of random
      ! Placing the new data at the end to easily access the beginning
      ! with a sum
      do i_eof = 1, eof_l
        rnd_hist(lag_l,i_eof) =  sum(rnd_hist(1:lags(i_eof),i_eof)*rho(1:lags(i_eof),i_eof)) &
                               + sigma(i_eof) * new_rnd_set(i_eof)
      enddo

      ! Then roll the history data in the lag direction so that our newly computed
      ! data is now at the front
      rnd_hist = Cshift(rnd_hist,-1,1)
    endif

    ! Assemble new field from EOFs and weights stored in history
    do i_eof = 1, eof_l
      stoch_fields(:,:,:,1) = stoch_fields(:,:,:,1) + eofs(:,:,:,i_eof) * rnd_hist(1+offset,i_eof)
    enddo

  end subroutine generate_new_forcing

  subroutine get_linear_coeff(sf_time, pold, pnew)

    implicit none

    ! IN
    real (r8), dimension(1:12), intent(in) :: sf_time

    ! OUT
    real (r8), intent(out) :: pold, pnew

    ! LOCAL
    real (r8) :: t_old, t_new           ! Old and new time in hour
    integer (int_kind) :: m_old, m_new  ! Old and new month idx

    m_old = imonth
    m_new = mod(m_old + 1, 12)
    if (m_new <= 0) m_new = m_new + 12

    t_old = sf_time(m_old)
    t_new = sf_time(m_new)

    ! sf_time was already update in Dec.
    ! Flip back to a year earlier.
    if (imonth == 12) then
      t_old = t_old - hours_in_year
    endif

    pold = (t_new - thour00) / (t_new - t_old)
    pnew = max(0.0, 1.0 - pold)
  end subroutine get_linear_coeff

  subroutine get_rnd_stoch(rnd, rnd_l, ascii_file)
    implicit none

    ! IN
    integer (int_kind), intent(in) :: rnd_l
    character(*), intent(in) :: ascii_file
    ! OUT
    real (r8), dimension(rnd_l), intent(out) :: rnd
    ! LOCAL
    integer (int_kind) :: funit = 21
    logical :: fcheck

    if (my_task == master_task) then
      inquire( file=trim(ascii_file), exist=fcheck )
      if (.not. fcheck) then
        write(*,*) "Unable to open random number file: ", trim(ascii_file)
      else
        open(funit, file=trim(ascii_file), status="old", action="read")
        read(funit, fmt=*) rnd(:)
        close(funit)
      endif
    endif

    call broadcast_scalar(fcheck, master_task)
    if (.not. fcheck) then
      call exit_POP(sigAbort,'error while reading random numbers data file')
    endif

    call broadcast_array(rnd, master_task)
  end subroutine get_rnd_stoch

end module forcing_stoch
