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
  real (r8), allocatable, dimension(:,:,:,:) :: &
    stoch_data_ep,             &
    stoch_data_t2m

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
                                      lags_ep, sig_ep, rho_ep, hist_ep)

        allocate(eof_ep(nx_block,ny_block,max_blocks_clinic,n_eof_ep))

        EOFfile = construct_file("nc", &
                                full_name=trim(EOF_ep_filename))

        call read_EOF_netcdf_file(EOFfile, n_eof_ep, eof_ep)

        allocate(stoch_data_ep(nx_block,ny_block,max_blocks_clinic,2))
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

        allocate(lags_t2m(1:n_eof_t2m))
        allocate(sig_t2m(1:n_eof_t2m))
        allocate(rho_t2m(1:lag_max_t2m,1:n_eof_t2m))
        allocate(hist_t2m(1:lag_max_t2m,1:n_eof_t2m))

        allocate(rnd_t2m(1:n_eof_t2m))

        call read_AR_netcdf_file_data(ARfile,"t2m",n_eof_t2m,lag_max_t2m,&
                                      lags_t2m, sig_t2m, rho_t2m, hist_t2m)

        allocate(eof_t2m(nx_block,ny_block,max_blocks_clinic,n_eof_t2m))

        EOFfile = construct_file("nc", &
                                full_name=trim(EOF_t_filename))

        call read_EOF_netcdf_file(EOFfile, n_eof_t2m, eof_t2m)
        
        allocate(stoch_data_t2m(nx_block,ny_block,max_blocks_clinic,2))

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

    ! Compute a new random monthly field if needed
    call update_stoch_forcing_data_fwf()

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

   ! Compute the stoicastic heat forcing
   call calc_stochastic_shf(STF, STF_stoch)

   ! Add stochastic forcing to baseline forcing
   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
      STF(:,:,1,iblock) = STF(:,:,1,iblock) + STF_stoch(:,:,1,iblock)
   enddo
   !$OMP END PARALLEL DO

  end subroutine append_stoch_forcing_shf

  subroutine update_stoch_forcing_data_fwf()
    implicit none
  end subroutine update_stoch_forcing_data_fwf

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
      !write(*,*) " >> Integral of stoch SFWF after correction and scaling: ", stf_int

    case ('ERA5-Data')

       do iblock = 1, nblocks_clinic
         STF_stoch(:,:,2,iblock) = eof_ep(:,:,iblock,1)*salinity_factor ! kg/s/m^2
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

    !$OMP PARALLEL DO PRIVATE(iblock)
    do iblock = 1, nblocks_clinic
       STF_stoch(:,:,1,iblock) = c0
    end do
    !$OMP END PARALLEL DO

  end subroutine calc_stochastic_shf

end module forcing_stoch
