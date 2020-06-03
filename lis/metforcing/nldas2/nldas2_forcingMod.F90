!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nldas2_forcingMod
!BOP
! !MODULE: nldas2_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data used in the North American
!  Land Data Assimilation System Phase II.  The variables are produced 
!  at 0.125 degree spatial resolution, and at hourly intervals.  For more
!  details please view the forcing files manual available at the 
!  following URL:
!
!  http://ldas.gsfc.nasa.gov//nldas/NLDAS2forcing.php
! 
!  The implemenatation in LIS has the derived data type {\tt nldas2\_struc}
!  that includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nldas2time1]
!    The nearest, previous hourly instance of the incoming 
!    data (as a real time). 
!  \item[nldas2time2]
!    The nearest, next hourly instance of the incoming 
!    data (as a real time).
!  \item[nldas2dir]
!    Directory containing the input data
!  \item[nldas2\_filesrc]
!    Center(GES-DISC|NCEP)-based NLDAS-2 filename source option
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for nearest neighbor interpolation. 
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !REVISION HISTORY: 
! 02 Feb 2004: Sujay Kumar; Initial Specification
! 24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
! 14 Mar 2014: David Mocko: Added CAPE and PET forcing from NLDAS-2
! 
! !USES: 
      use ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_NLDAS2      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nldas2_struc
  public :: num_nldas2_fields
  public :: list_nldas2_fields
!EOP


  type, public ::  nldas2_type_dec 
     real          :: ts
     integer       :: ncold, nrold   ! AWIPS 212 dimensions
     character*50  :: nldas2_filesrc
     character*80  :: nldas2dir      ! NLDAS-2 Forcing Directory
     real*8        :: nldas2time1,nldas2time2
     integer       :: model_level_data 
     integer       :: model_level_press 
     integer       :: model_pcp_data 
     integer       :: model_dswrf_data 
     
     real,  allocatable     :: orig_ediff(:)

     real                   :: gridDesc(50)
     integer                :: mi
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)

     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)

     integer, allocatable   :: n113(:)
     integer                :: findtime1, findtime2

     integer           :: nIter, st_iterid,en_iterid  ! Forecast parameters

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

     ! For ESMF regridding
     type(ESMF_FieldBundle)       :: forcing_bundle
     type(ESMF_RouteHandle)       :: routehandle
     type(ESMF_DynamicMask)       :: dynamicMask
     type(ESMF_TypeKind_Flag)     :: type_kind = ESMF_TYPEKIND_R4
     type(ESMF_STAGGERLOC)        :: staggerloc
     type(ESMF_RegridMethod_Flag) :: regridMethod
     real                         :: undefined_value ! for missing value
     
  end type nldas2_type_dec

  integer            :: num_nldas2_fields       ! number of available fields
  character(len=100) :: list_nldas2_fields(30)  ! list of name of fields

  type(nldas2_type_dec), allocatable :: nldas2_struc(:)
  character(len=30), parameter :: nldas2_bundle_bname = "nldas2_bundle_"
!EOP
contains
  
!BOP
!
! !ROUTINE: init_NLDAS2
! \label{init_NLDAS2}
!
! !INTERFACE:
  subroutine init_NLDAS2(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain, LIS_vm, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit,LIS_endrun, LIS_verify
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use map_utils,      only : proj_latlon
    use LIS_forecastMod

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for NLDAS-2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_nldas2](\ref{readcrd_nldas2}) \newline
!     reads the runtime options specified for NLDAS-2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_nldas2\_elev](\ref{read_nldas2_elev}) \newline
!    reads the native elevation of the NLDAS-2 data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP
    
    integer          :: n
    integer          :: ic, rc
    real             :: forcing_gridDesc(10)
    
    allocate(nldas2_struc(LIS_rc%nnest))
    call readcrd_nldas2()

    do n=1, LIS_rc%nnest
       nldas2_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, nldas2_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 13

  ! Set NLDAS-2 grid dimensions and extent information:
    nldas2_struc(:)%ncold = 464
    nldas2_struc(:)%nrold = 224

    IF (LIS_rc%do_esmfRegridding) THEN
       CALL set_list_nldas2_fields()
    ENDIF

    do n=1,LIS_rc%nnest


       ! Forecast mode:
       if(LIS_rc%forecastMode.eq.1) then

          if(mod(LIS_rc%nensem(n),&
               LIS_forecast_struc(1)%niterations).ne.0) then
             write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
             write(LIS_logunit,*) '[ERR] of the number of iterations '
             write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
             write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
             call LIS_endrun()
          endif

          allocate(nldas2_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
                   LIS_rc%met_nf(findex),&
                   LIS_rc%ngrid(n)))
          allocate(nldas2_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
                   LIS_rc%met_nf(findex),&
                   LIS_rc%ngrid(n)))

          nldas2_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          nldas2_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          nldas2_struc(n)%nIter = LIS_forecast_struc(1)%niterations

       ! Regular retrospective or non-forecast mode:
       else
          allocate(nldas2_struc(n)%metdata1(1,&
                   LIS_rc%met_nf(findex),&
                   LIS_rc%ngrid(n)))
          allocate(nldas2_struc(n)%metdata2(1,&
                   LIS_rc%met_nf(findex),&
                   LIS_rc%ngrid(n)))

          nldas2_struc(n)%st_iterid = 1
          nldas2_struc(n)%en_iterId = 1
          nldas2_struc(n)%nIter = 1

       endif

       nldas2_struc(n)%metdata1 = 0
       nldas2_struc(n)%metdata2 = 0

       nldas2_struc(n)%gridDesc = 0        
       nldas2_struc(n)%findtime1 = 0 
       nldas2_struc(n)%findtime2 = 0 

       nldas2_struc(n)%gridDesc(1) = 0
       nldas2_struc(n)%gridDesc(2) = nldas2_struc(n)%ncold
       nldas2_struc(n)%gridDesc(3) = nldas2_struc(n)%nrold
       nldas2_struc(n)%gridDesc(4) = 25.0625
       nldas2_struc(n)%gridDesc(5) = -124.9375
       nldas2_struc(n)%gridDesc(6) = 128
       nldas2_struc(n)%gridDesc(7) = 52.9375
       nldas2_struc(n)%gridDesc(8) = -67.0625
       nldas2_struc(n)%gridDesc(9) = 0.125
       nldas2_struc(n)%gridDesc(10) = 0.125
       nldas2_struc(n)%gridDesc(20) = 64

     ! Check for grid and interp option selected:
       if( nldas2_struc(n)%gridDesc(9)  == LIS_rc%gridDesc(n,9) .and. &
           nldas2_struc(n)%gridDesc(10) == LIS_rc%gridDesc(n,10).and. &
           LIS_rc%gridDesc(n,1) == proj_latlon .and. &
           LIS_rc%met_interp(findex) .ne. "neighbor" ) then
         write(LIS_logunit,*) "[ERR] The NLDAS2 grid was selected for the"
         write(LIS_logunit,*) "[ERR] LIS run domain; however, 'bilinear', 'budget-bilinear',"
         write(LIS_logunit,*) "[ERR] or some other unknown option was selected to spatially"
         write(LIS_logunit,*) "[ERR] downscale the grid, which will cause errors during runtime."
         write(LIS_logunit,*) "[ERR] Program stopping ..."
         call LIS_endrun()
       endif

       nldas2_struc(n)%mi = nldas2_struc(n)%ncold*nldas2_struc(n)%nrold

       IF (LIS_rc%do_esmfRegridding) THEN
          !nldas2_struc(n)%type_kind       = ESMF_TYPEKIND_R4
          nldas2_struc(n)%undefined_value = 9999.0

          nldas2_struc(n)%regridMethod = ESMF_REGRIDMETHOD_BILINEAR ! ESMF_REGRIDMETHOD_NEAREST_STOD
          nldas2_struc(n)%staggerloc = ESMF_STAGGERLOC_CENTER
          LIS_domain(n)%staggerloc   = ESMF_STAGGERLOC_CENTER

          forcing_gridDesc(:) = 0
          forcing_gridDesc(2) = nldas2_struc(n)%gridDesc(2) ! num points along x
          forcing_gridDesc(3) = nldas2_struc(n)%gridDesc(3) ! num points along y
          forcing_gridDesc(4) = nldas2_struc(n)%gridDesc(4) ! lower lat
          forcing_gridDesc(5) = nldas2_struc(n)%gridDesc(5) ! lower lon
          forcing_gridDesc(7) = nldas2_struc(n)%gridDesc(7) ! upper lat
          forcing_gridDesc(8) = nldas2_struc(n)%gridDesc(8) ! upper lon
          forcing_gridDesc(9) = nldas2_struc(n)%gridDesc(9) ! x-grid size
          forcing_gridDesc(10)= nldas2_struc(n)%gridDesc(10) ! y-grid size

          CALL create_nldas2_Forcing_ESMFbundle(n, forcing_gridDesc(1:10))
          CALL create_nldas2_Model_ESMFbundle(n)
          CALL create_nldas2_ESMFroutehandle(n)
       ELSE
          ! Setting up weights for spatial interpolation:
          select case( LIS_rc%met_interp(findex) )
           case( "bilinear" )
             allocate(nldas2_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

             call bilinear_interp_input(n,&
                  nldas2_struc(n)%gridDesc(:),&
                  nldas2_struc(n)%n111,nldas2_struc(n)%n121,nldas2_struc(n)%n211,&
                  nldas2_struc(n)%n221,nldas2_struc(n)%w111,nldas2_struc(n)%w121,&
                  nldas2_struc(n)%w211,nldas2_struc(n)%w221)
   
           case( "budget-bilinear" )

             allocate(nldas2_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(nldas2_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

             call bilinear_interp_input(n,nldas2_struc(n)%gridDesc(:),&
                  nldas2_struc(n)%n111,nldas2_struc(n)%n121,&
                  nldas2_struc(n)%n211,nldas2_struc(n)%n221,&
                  nldas2_struc(n)%w111,nldas2_struc(n)%w121,&
                  nldas2_struc(n)%w211,nldas2_struc(n)%w221)
   
             allocate(nldas2_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(nldas2_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(nldas2_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(nldas2_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(nldas2_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(nldas2_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(nldas2_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(nldas2_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

             call conserv_interp_input(n,nldas2_struc(n)%gridDesc(:),&
                  nldas2_struc(n)%n112,nldas2_struc(n)%n122,&
                  nldas2_struc(n)%n212,nldas2_struc(n)%n222,&
                  nldas2_struc(n)%w112,nldas2_struc(n)%w122,&
                  nldas2_struc(n)%w212,nldas2_struc(n)%w222)

           case( "neighbor" )

             allocate(nldas2_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
             call neighbor_interp_input(n,nldas2_struc(n)%gridDesc(:),&
                  nldas2_struc(n)%n113)
           case default
             write(LIS_logunit,*) "[ERR] Interpolation option not specified for NLDAS2"
             write(LIS_logunit,*) "[ERR] Program stopping ..."
             call LIS_endrun()
          end select
       ENDIF

     ! Read in elevation difference and NLDAS2 elevation maps:
       if( LIS_rc%met_ecor(findex).ne."none" ) then 

          allocate(nldas2_struc(n)%orig_ediff(&
               nldas2_struc(n)%ncold*nldas2_struc(n)%nrold))

          call read_orig_nldas2_elevdiff(n)
          call read_nldas2_elev(n,findex)
       endif

     ! Set up precipitation climate downscaling:
       if(LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,&
               nint(nldas2_struc(n)%gridDesc(2)),&
               nint(nldas2_struc(n)%gridDesc(3)))
          
       endif

    enddo

  end subroutine init_NLDAS2

!------------------------------------------------------------------------------
!BOP
      subroutine set_list_nldas2_fields()
!
! !USES:
      use LIS_FORC_AttributesMod
!
! !DESCRIPTION:
! Determine the number of available fields that will be regridded
!
! !LOCAL VARIABLES:
      integer :: ic
!EOP
!------------------------------------------------------------------------------
!BOC
      ic = 0
      if (LIS_FORC_Tair%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_Tair%varname(1))
      endif
      if (LIS_FORC_Qair%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_Qair%varname(1))
      endif
      if (LIS_FORC_SWdown%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_SWdown%varname(1))
      endif
      if (LIS_FORC_LWdown%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_LWdown%varname(1))
      endif
      if (LIS_FORC_Wind_E%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_Wind_E%varname(1))
      endif
      if (LIS_FORC_Wind_N%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_Wind_N%varname(1))
      endif
      if (LIS_FORC_Psurf%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_Psurf%varname(1))
      endif
      if (LIS_FORC_Rainf%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_Rainf%varname(1))
      endif
      if (LIS_FORC_CRainf%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_CRainf%varname(1))
      endif
      if (LIS_FORC_PET%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_PET%varname(1))
      endif
      if (LIS_FORC_CAPE%selectOpt .eq. 1) then
         ic = ic + 1
         list_nldas2_fields(ic) = TRIM(LIS_FORC_CAPE%varname(1))
      endif


      num_nldas2_fields = ic

      end subroutine set_list_nldas2_fields
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_nldas2_Forcing_ESMFbundle(n, forcing_gridDesc)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod,      only : create_regular_grid
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
      real,    intent(in) :: forcing_gridDesc(10)
!
! !DESCRIPTION:
! Create the NLDAS2 forcing ESMF grid and bundle.
! This subroutine might be called several times depending on the integration date.
! Howver, it will be called once for any period when the NLDAS2 forcing resolution
! does not change.
!
! !LOCAL VARIABLES:
      character(len=2) :: num_st
      integer          :: rc, ic
      type(ESMF_Grid)  :: nldas2_grid
      real             :: dummy_array(1,1) = 0.0
!EOP
!------------------------------------------------------------------------------
!BOC
     write(LIS_logunit,*) '[INFO] Initialize NLDAS2 ESMF object.'
      ! Create the bundles
      write(num_st, '(i2.2)') n

      ! ---> Bundle for nldas2
      nldas2_struc(n)%forcing_bundle = ESMF_FieldBundleCreate(name = TRIM(nldas2_bundle_bname)//num_st, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleCreate failed for nldas2 Forcing Data')

      nldas2_grid =   create_regular_grid(forcing_gridDesc, "nldas2 Grid", &
                                        LIS_rc%npesx, LIS_rc%npesy, &
                                        staggerloc = nldas2_struc(n)%staggerloc)
      ! Add fields to the bundle
      DO ic = 1, num_nldas2_fields
         if (LIS_masterproc) PRINT*,"--->Adding-Forcing: ", ic, TRIM(list_nldas2_fields(ic))
         call addTracerToBundle(nldas2_struc(n)%forcing_bundle, nldas2_grid, &
                   TRIM(list_nldas2_fields(ic)), nldas2_struc(n)%type_kind, dummy_array)
      ENDDO

      end subroutine create_nldas2_Forcing_ESMFbundle
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_nldas2_Model_ESMFbundle(n)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod,       only : create_regular_grid
!
! !INPUT PRAMETERS:
      integer, intent(in) :: n
!
! !DESCRIPTION:
! Create the model ESMF grid and bundle.
! This subroutine is called once as the model resolution does not change.
!
! !LOCAL VARIABLES:
      character(len=2) :: num_st
      integer          :: rc, ic
      real             :: model_gridDesc(10)
      type(ESMF_Grid)  :: model_grid
      real             :: dummy_array(1,1) = 0.0
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize model ESMF objects.'

      ! Create the bundles
      write(num_st, '(i2.2)') n

      ! ---> Bundle for the model
      LIS_domain(n)%nldas2_bundle = ESMF_FieldBundleCreate(name = TRIM(nldas2_bundle_bname)//num_st, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleCreate failed for model Data')

      ! ---> ESMF grid for for the model
      model_gridDesc(:) = 0
      model_gridDesc(2) = LIS_rc%gnc(n)         ! LIS_rc%gridDesc(n,2) ! Global num points along x
      model_gridDesc(3) = LIS_rc%gnr(n)         ! LIS_rc%gridDesc(n,3) ! Global num points along y
      model_gridDesc(4) = LIS_rc%gridDesc(n,34) ! lower lat
      model_gridDesc(5) = LIS_rc%gridDesc(n,35) ! lower lon
      model_gridDesc(7) = LIS_rc%gridDesc(n,37) ! 7) ! upper lat
      model_gridDesc(8) = LIS_rc%gridDesc(n,38) ! 8) ! upper lon
      model_gridDesc(9) = LIS_rc%gridDesc(n,39) ! x-grid size
      model_gridDesc(10)= LIS_rc%gridDesc(n,40) ! y-grid size

      model_grid =  create_regular_grid(model_gridDesc, "Model Grid", &
                                        LIS_rc%npesx, LIS_rc%npesy, &
                                        staggerloc = LIS_domain(n)%staggerloc)

      ! Add fields to the bundle
      DO ic = 1, num_nldas2_fields
         if (LIS_masterproc) PRINT*,"<---Adding-Model: ", ic, TRIM(list_nldas2_fields(ic))
         call addTracerToBundle(LIS_domain(n)%nldas2_bundle, model_grid, &
                   TRIM(list_nldas2_fields(ic)), nldas2_struc(n)%type_kind, dummy_array)
      ENDDO

      end subroutine create_nldas2_Model_ESMFbundle
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_nldas2_ESMFroutehandle(n)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun

      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_ESMF_Regrid_Utils, only : createESMF_RouteHandle
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
! 
! !DESCRIPTION:
! Determine the ESMF routehandle needed for thr regridding between 
! the NLDAS2 forcing and the model.
!
! !LOCAL VARIABLES:
      type(ESMF_FIELD) :: model_field, nldas2_field
      integer          :: ftc(2), ftlb(2), ftub(2)
      real, pointer    :: farray2dd(:,:)
      integer          :: rc
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Determine the ESMF routehandle.'

      ! Get one field from the NLDAS2 focing bundle
      call ESMF_FieldBundleGet(nldas2_struc(n)%forcing_bundle, TRIM(list_nldas2_fields(1)), &
                              field = nldas2_field, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for nldas2 field'//TRIM(list_nldas2_fields(1)))

!    call ESMF_FieldGet(nldas2_field, localDe=0, farrayPtr=farray2dd, &
!        totalLBound=ftlb, totalUBound=ftub, totalCount=ftc, rc=rc)
!        print"(a14,7i7)", "ForcingGrid: ", LIS_localPet,ftlb,ftub,ftc
!      call LIS_verify(rc, 'ESMF_FieldGet failed for NLDAS2  field'//TRIM(list_nldas2_fields(1)))


      ! Get the corresponding field from the model bundle
      call ESMF_FieldBundleGet(LIS_domain(n)%nldas2_bundle, TRIM(list_nldas2_fields(1)), &
                               field = model_field, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for model field'//TRIM(list_nldas2_fields(1)))


!    call ESMF_FieldGet(model_field, localDe=0, farrayPtr=farray2dd, &
!        totalLBound=ftlb, totalUBound=ftub, totalCount=ftc, rc=rc)
!        print"(a14,7i7)", "ModelGrid: ", LIS_localPet,ftlb,ftub,ftc
!      call LIS_verify(rc, 'ESMF_FieldGet failed for model field'//TRIM(list_nldas2_fields(1)))


      ! Compute the ESMF routehandle
      call createESMF_RouteHandle(nldas2_field, model_field,  nldas2_struc(n)%regridMethod, &
                                  nldas2_struc(n)%undefined_value, nldas2_struc(n)%routehandle, &
                                  nldas2_struc(n)%dynamicMask, rc)
      call LIS_verify(rc, 'createESMF_RouteHandle failed')

      end subroutine create_nldas2_ESMFroutehandle
!EOC
!------------------------------------------------------------------------------
end module nldas2_forcingMod
