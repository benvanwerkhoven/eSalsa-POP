module forcing_stoich

!BOP
! !MODULE: forcing_stoich
!
! !DESCRIPTION:
!  This module contains the surface stoichastic forcing
!  data and methods.
!
! !REVISION HISTORY:
!
! !USES:

  use kinds_mod
  use blocks
  use distribution
  use domain
  use constants
  use prognostic
  use io
  use grid
  use io_netcdf
  use global_reductions, only:global_sum
  use exit_mod

  implicit none

  real (r8), public :: sf_fwf_fraction = 0.0

  logical :: stoichastic_forcing_fwf
  logical :: stoichastic_forcing_hf

  character (char_len) :: &
    sf_fwf_formulation,         &! Formulation of the freshwater stoich forcing
    sf_hf_formulation            ! Formulation of the heat stoich forcing

  character (char_len) :: &
    ARdata_ep_filename,         &! Autoregression datafile for E-P
    ARdata_t_filename,          &! Autoregression datafile for temperature
    EOF_ep_filename,            &! EOF datafile for E-P
    EOF_t_filename               ! EOF datafile for temperature

contains

  subroutine read_stoich_forcing_namelist
    implicit none
    integer(int_kind) :: nml_error ! namelist error flag
    namelist /stoich_forcing/ stoichastic_forcing_fwf, stoichastic_forcing_hf,  &
                              sf_fwf_formulation, sf_hf_formulation,            &
                              ARdata_ep_filename, ARdata_t_filename,            &
                              EOF_ep_filename, EOF_t_filename,                  &
                              sf_fwf_fraction

    ! Set some defaults
    stoichastic_forcing_fwf   = .false.
    stoichastic_forcing_hf    = .false.
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
          read(nml_in, nml=stoich_forcing,iostat=nml_error)
       end do
       if (nml_error == 0) close(nml_in)
    endif

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading forcing_sfwf_nml')
    endif

    call broadcast_scalar(stoichastic_forcing_fwf, master_task)
    call broadcast_scalar(stoichastic_forcing_hf,  master_task)
    call broadcast_scalar(sf_fwf_formulation,      master_task)
    call broadcast_scalar(sf_hf_formulation,       master_task)
    call broadcast_scalar(ARdata_ep_filename,      master_task)
    call broadcast_scalar(ARdata_t_filename,       master_task)
    call broadcast_scalar(sf_fwf_fraction,         master_task)

  end subroutine read_stoich_forcing_namelist

  subroutine init_stoich_forcing(STF_stoich)

    implicit none

    ! !INOUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(inout) :: &
       STF_stoich! surface tracer fluxes for all tracers

    STF_stoich(:,:,:,:) = c0

  end subroutine init_stoich_forcing

  subroutine read_netCDF_ARdata()

    implicit none

  end subroutine read_netCDF_ARdata

  subroutine append_stoich_forcing_sfwf(STF)

    use forcing_fields, only : STF_stoich

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(inout) :: &
       STF     ! surface tracer fluxes for all tracers

    ! LOCAL
    integer (int_kind) :: &
       iblock              ! block loop index

    ! Return if sfwf not activated
    if (.not. stoichastic_forcing_fwf) then
      return
    endif

   ! Compute the stoicastic freshwater forcing
   call calc_stoichastic_sfwf(STF, STF_stoich)

   ! Add stoichastic forcing to baseline forcing
   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
      STF(:,:,2,iblock) = STF(:,:,2,iblock) + STF_stoich(:,:,2,iblock)
   enddo
   !$OMP END PARALLEL DO

  end subroutine append_stoich_forcing_sfwf

  subroutine append_stoich_forcing_shf(STF)

    use forcing_fields, only : STF_stoich

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
       intent(inout) :: &
       STF     ! surface tracer fluxes for all tracers

    ! LOCAL
    integer (int_kind) :: &
       iblock              ! block loop index

    ! Return if stoich shf not activated
    if (.not. stoichastic_forcing_hf) then
      return
    endif

   ! Compute the stoicastic heat forcing
   call calc_stoichastic_shf(STF, STF_stoich)

   ! Add stoichastic forcing to baseline forcing
   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
      STF(:,:,1,iblock) = STF(:,:,1,iblock) + STF_stoich(:,:,1,iblock)
   enddo
   !$OMP END PARALLEL DO

  end subroutine append_stoich_forcing_shf

  subroutine calc_stoichastic_sfwf(STF, STF_stoich)
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
       STF_stoich! surface tracer fluxes for all tracers

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
      ! Stoichastic forcing is based on a fraction of the
      ! baseline forcing computed previously.
      ! 1. Mask STF on northern atlantic
      ! 2. Compute mean of masked STF
      ! 3. STF_stoich = sf_fwf_fraction * (STF_m - \int(STF_m))
      ! 4. STF += STF_stoich

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
                  STF_stoich(i,j,2,iblock) = STF(i,j,2,iblock)
                  stf_int_lcl = stf_int_lcl + TAREA(i,j,iblock) * 1.0e-10_r8 * STF_stoich(i,j,2,iblock)
                  surf_lcl = surf_lcl + TAREA(i,j,iblock) * 1.0e-10_r8
               else
                  STF_stoich(i,j,2,iblock) = c0
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      stf_int = global_sum(stf_int_lcl, distrb_clinic)
      !write(*,*) " >> Integral of stoich SFWF before correction: ", stf_int
      surf = global_sum(surf_lcl, distrb_clinic)
      !write(*,*) " >> Surface covered by stoich SFWF: ", surf
      stf_int = stf_int / surf
      !write(*,*) " >> Mean stoich SFWF: ", stf_int

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
                  STF_stoich(i,j,2,iblock) = sf_fwf_fraction * (STF_stoich(i,j,2,iblock) - stf_int)
                  stf_int_lcl = stf_int_lcl + TAREA(i,j,iblock) * 1.0e-10_r8 * STF_stoich(i,j,2,iblock)
               else
                  STF_stoich(i,j,2,iblock) = c0
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      stf_int = global_sum(stf_int_lcl, distrb_clinic)
      !write(*,*) " >> Integral of stoich SFWF after correction and scaling: ", stf_int

    case ('ERA5-Data')
      call abort()
    end select

  end subroutine calc_stoichastic_sfwf

  subroutine calc_stoichastic_shf(STF, STF_stoich)
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
       STF_stoich ! surface tracer fluxes for all tracers

    ! LOCAL
    integer (int_kind) :: &
       iblock              ! block loop index

    !$OMP PARALLEL DO PRIVATE(iblock)
    do iblock = 1, nblocks_clinic
       STF_stoich(:,:,1,iblock) = c0
    end do
    !$OMP END PARALLEL DO

  end subroutine calc_stoichastic_shf

end module forcing_stoich
