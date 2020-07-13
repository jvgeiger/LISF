!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module merra2_forcingMod
!BOP
! !MODULE: merra2_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the MERRA2 forcing data.
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt merra2\_struc}
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the ECMWF data
!  \item[merra2time1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[merra2time2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[merra2dir]
!    Directory containing the input data
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
!    for each grid point in LIS, for n. neighbor interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \end{description}
!
! !USES:
      USE ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_merra2      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: merra2_struc
  public :: num_merra2_fields
  public :: list_merra2_fields

!EOP
  type, public ::  merra2_type_dec
     real         :: ts
     integer      :: ncold, nrold
     character*40 :: merra2dir   !MERRA2 Forcing Directory
     real*8       :: merra2time1,merra2time2
     logical      :: reset_flag

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
     logical                :: startFlag, dayFlag
     real, allocatable      :: merraforc1(:,:,:,:), merraforc2(:,:,:,:)

     integer            :: nvars
     integer            :: uselml
     integer            :: usecorr

     real*8             :: ringtime
     
     integer            :: nIter, st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

     integer                 :: usescalef
     integer                 :: usepcpsampling
     integer                 :: pcpscal_cmo
     integer                 :: use2mwind
     character*100           :: scaleffile
     integer                 :: nbins
     real, allocatable       :: refxrange(:,:,:,:)
     real, allocatable       :: refcdf(:,:,:,:)
     real, allocatable       :: refmean(:,:,:)
     real, allocatable       :: refmean_ip(:)
     real, allocatable       :: refstdev(:,:,:)
     real, allocatable       :: refstdev_ip(:)
     real, allocatable       :: merraxrange(:,:,:,:)
     real, allocatable       :: merracdf(:,:,:,:)
     integer, allocatable    :: rseed(:,:)

     ! For ESMF regridding
     type(ESMF_FieldBundle)       :: forcing_bundle
     type(ESMF_RouteHandle)       :: routehandle
     type(ESMF_DynamicMask)       :: dynamicMask
     type(ESMF_TypeKind_Flag)     :: type_kind = ESMF_TYPEKIND_R4
     type(ESMF_STAGGERLOC)        :: staggerloc
     type(ESMF_RegridMethod_Flag) :: regridMethod
     type(ESMF_CoordSys_Flag)     :: coordSys
     real                         :: undefined_value ! for missing value

  end type merra2_type_dec

  type(merra2_type_dec), allocatable :: merra2_struc(:)

  integer            :: num_merra2_fields       ! number of available fields
  character(len=100) :: list_merra2_fields(30)  ! list of name of fields
  character(len=30), parameter :: merra2_bundle_bname = "merra2_bundle_"

contains

!BOP
!
! !ROUTINE: init_merra2
! \label{init_merra2}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
!
! !INTERFACE:
  subroutine init_merra2(findex)

! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use LIS_forecastMod

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for MERRA2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_merra2](\ref{readcrd_merra2}) \newline
!     reads the runtime options specified for MERRA2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real :: gridDesci(LIS_rc%nnest,50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n

    integer          :: ic, rc
    real             :: forcing_gridDesc(10)

    allocate(merra2_struc(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       merra2_struc(n)%ncold = 576
       merra2_struc(n)%nrold = 361
    enddo

    call readcrd_merra2()
    LIS_rc%met_nf(findex) = 14

    merra2_struc%reset_flag = .false.

    do n=1, LIS_rc%nnest
       merra2_struc(n)%ts = 3600  !check
       call LIS_update_timestep(LIS_rc, n, merra2_struc(n)%ts)
    enddo

    IF (LIS_rc%do_esmfRegridding) THEN
       CALL set_list_merra2_fields(findex)
    ENDIF


    gridDesci = 0

    do n=1,LIS_rc%nnest
       gridDesci(n,1) = 0
       gridDesci(n,2) = merra2_struc(n)%ncold
       gridDesci(n,3) = merra2_struc(n)%nrold
       gridDesci(n,4) = -90.000
       gridDesci(n,5) = -180.000
       gridDesci(n,6) = 128
       gridDesci(n,7) = 90.000
       gridDesci(n,8) = 179.375
       gridDesci(n,9) = 0.625
       gridDesci(n,10) = 0.5
       gridDesci(n,20) = 0

       merra2_struc(n)%mi = merra2_struc(n)%ncold*merra2_struc(n)%nrold

       IF (LIS_rc%do_esmfRegridding) THEN
          !merra2_struc(n)%type_kind       = ESMF_TYPEKIND_R4
          merra2_struc(n)%undefined_value = LIS_rc%udef

          if ((LIS_rc%met_interp(findex)) .eq. "bilinear") then
             merra2_struc(n)%regridMethod = ESMF_REGRIDMETHOD_BILINEAR
          elseif(trim(LIS_rc%met_interp(findex)) .eq. "neighbor") then
             merra2_struc(n)%regridMethod = ESMF_REGRIDMETHOD_NEAREST_STOD     
          elseif(trim(LIS_rc%met_interp(findex)) .eq. "conservative") then
             merra2_struc(n)%regridMethod =  ESMF_REGRIDMETHOD_CONSERVE
          endif

          merra2_struc(n)%staggerloc = ESMF_STAGGERLOC_CENTER
          LIS_domain(n)%staggerloc   = ESMF_STAGGERLOC_CENTER
          merra2_struc(n)%coordSys   = ESMF_COORDSYS_SPH_DEG
          LIS_domain(n)%coordSys     = ESMF_COORDSYS_SPH_DEG

          forcing_gridDesc(:) = 0
          forcing_gridDesc(2) = gridDesci(n, 2) ! num points along x
          forcing_gridDesc(3) = gridDesci(n, 3) ! num points along y
          forcing_gridDesc(4) = gridDesci(n, 4) ! lower lat
          forcing_gridDesc(5) = gridDesci(n, 5) ! lower lon
          forcing_gridDesc(7) = gridDesci(n, 7) ! upper lat
          forcing_gridDesc(8) = gridDesci(n, 8) ! upper lon
          forcing_gridDesc(9) = gridDesci(n, 9) ! x-grid size
          forcing_gridDesc(10)= gridDesci(n, 10) ! y-grid size

          CALL create_merra2_Forcing_ESMFbundle(n, forcing_gridDesc(1:10))
          CALL create_merra2_Model_ESMFbundle(n)
          CALL create_merra2_ESMFroutehandle(n)
       ELSE
          ! Setting up weights for Interpolation
          if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
             allocate(merra2_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call bilinear_interp_input(n, gridDesci(n,:),&
                  merra2_struc(n)%n111,merra2_struc(n)%n121,&
                  merra2_struc(n)%n211,merra2_struc(n)%n221,&
                  merra2_struc(n)%w111,merra2_struc(n)%w121,&
                  merra2_struc(n)%w211,merra2_struc(n)%w221)

          elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
             allocate(merra2_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call bilinear_interp_input(n, gridDesci(n,:),&
                  merra2_struc(n)%n111,merra2_struc(n)%n121,&
                  merra2_struc(n)%n211,merra2_struc(n)%n221,&
                  merra2_struc(n)%w111,merra2_struc(n)%w121,&
                  merra2_struc(n)%w211,merra2_struc(n)%w221)

             allocate(merra2_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             call conserv_interp_input(n, gridDesci(n,:),&
                  merra2_struc(n)%n112,merra2_struc(n)%n122,&
                  merra2_struc(n)%n212,merra2_struc(n)%n222,&
                  merra2_struc(n)%w112,merra2_struc(n)%w122,&
                  merra2_struc(n)%w212,merra2_struc(n)%w222)

          elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
             allocate(merra2_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call neighbor_interp_input(n, gridDesci(n,:),&
                  merra2_struc(n)%n113)
   
          else
             write(LIS_logunit,*) '[ERR] Interpolation option '// &
                  trim(LIS_rc%met_interp(findex))//&
                  ' for MERRA2 forcing is not supported'
             call LIS_endrun()
          endif
       ENDIF

       call LIS_registerAlarm("MERRA2 forcing alarm", 86400.0, 86400.0)
       merra2_struc(n)%startFlag = .true.
       merra2_struc(n)%dayFlag = .true.

       merra2_struc(n)%nvars = 14

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

          allocate(merra2_struc(n)%merraforc1(&
               LIS_forecast_struc(1)%niterations,&
               merra2_struc(n)%nvars, 24, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2_struc(n)%merraforc2(&
               LIS_forecast_struc(1)%niterations,&
               merra2_struc(n)%nvars, 24, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          merra2_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          merra2_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          merra2_struc(n)%nIter = LIS_forecast_struc(1)%niterations
       
          allocate(merra2_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(merra2_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       ! Regular retrospective or non-forecast mode:
       else
          allocate(merra2_struc(n)%merraforc1(1,&
               merra2_struc(n)%nvars, 24, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2_struc(n)%merraforc2(1,&
               merra2_struc(n)%nvars, 24, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          merra2_struc(n)%st_iterid = 1
          merra2_struc(n)%en_iterId = 1
          merra2_struc(n)%nIter = 1
       
          allocate(merra2_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(merra2_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
       
       endif

       merra2_struc(n)%metdata1 = 0
       merra2_struc(n)%metdata2 = 0

       merra2_struc(n)%merraforc1 = LIS_rc%udef
       merra2_struc(n)%merraforc2 = LIS_rc%udef

       if ( LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
            LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" ) then
          call read_merra2_elev(n,findex)
       endif

       ! Set up precipitation climate downscaling:
       if(LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,&
               merra2_struc(n)%ncold,&
               merra2_struc(n)%nrold)
       endif

    enddo   ! End nest loop

  end subroutine init_merra2

!------------------------------------------------------------------------------
!BOP
      subroutine set_list_merra2_fields(findex)
!
! !USES:
      use LIS_FORC_AttributesMod
      use LIS_coreMod,    only : LIS_rc
!
! !INPUT PARAMETERS: 
      integer, intent(in) :: findex
!
! !DESCRIPTION:
! Determine the number of available fields that will be regridded
!
! !LOCAL VARIABLES:
      integer :: ic
      character(len=2) :: num_spc
!EOP
!------------------------------------------------------------------------------
!BOC
      ic = 0
      if (LIS_FORC_Tair%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Tair%varname(1))
      endif
      if (LIS_FORC_Qair%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Qair%varname(1))
      endif
      if (LIS_FORC_SWdown%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_SWdown%varname(1))
      endif
      if (LIS_FORC_LWdown%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_LWdown%varname(1))
      endif
      if (LIS_FORC_Wind_E%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Wind_E%varname(1))
      endif
      if (LIS_FORC_Wind_N%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Wind_N%varname(1))
      endif
      if (LIS_FORC_Psurf%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Psurf%varname(1))
      endif
      if (LIS_FORC_Rainf%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Rainf%varname(1))
      endif
      if (LIS_FORC_CRainf%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_CRainf%varname(1))
      endif
      if (LIS_FORC_Snowf%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Snowf%varname(1))
      endif
      if (LIS_FORC_Pardr%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Pardr%varname(1))
      endif
      if (LIS_FORC_Pardf%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Pardf%varname(1))
      endif
      if (LIS_FORC_Swnet%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Swnet%varname(1))
      endif
      if (LIS_FORC_Forc_Hgt%selectOpt .eq. 1) then
         ic = ic + 1
         list_merra2_fields(ic) = TRIM(LIS_FORC_Forc_Hgt%varname(1))
      endif

      num_merra2_fields = ic

      num_merra2_fields = LIS_rc%met_nf(findex)
      DO ic = 1, num_merra2_fields
         write(num_spc, '(i2.2)') ic
         list_merra2_fields(ic) = "merra2_var_"//num_spc
      ENDDO

      end subroutine set_list_merra2_fields
!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
      subroutine create_merra2_Forcing_ESMFbundle(n, forcing_gridDesc)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod !,      only : create_regular_grid
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
      real,    intent(in) :: forcing_gridDesc(10)
!
! !DESCRIPTION:
! Create the MERRA2 forcing ESMF grid and bundle.
! This subroutine might be called several times depending on the integration date.
! Howver, it will be called once for any period when the MERRA2 forcing resolution
! does not change.
!
! !LOCAL VARIABLES:
      character(len=2) :: num_st
      integer          :: rc, ic
      type(ESMF_Grid)  :: merra2_grid
      real             :: dummy_array(1,1) = 0.0
      real, pointer    :: lat_points(:)
      real, pointer    :: lon_points(:)
      real             :: min_lon, min_lat, dx, dy
      integer          :: num_lons, num_lats
!EOP
!------------------------------------------------------------------------------
!BOC
     write(LIS_logunit,*) '[INFO] Initialize MERRA2 ESMF object.'
      ! Create the bundles
      write(num_st, '(i2.2)') n

      ! ---> Bundle for merra2
      merra2_struc(n)%forcing_bundle = ESMF_FieldBundleCreate(name = TRIM(merra2_bundle_bname)//num_st, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleCreate failed for merra2 Forcing Data')

      num_lons = forcing_gridDesc(2)
      min_lon  = forcing_gridDesc(5)
      dx       = forcing_gridDesc(9)

      num_lats = forcing_gridDesc(3)
      min_lat  = forcing_gridDesc(4)
      dy       = forcing_gridDesc(10)

      ALLOCATE(lon_points(num_lons))
      do ic = 1, num_lons
         lon_points(ic) = (ic-1)*dx + min_lon
      enddo

      ALLOCATE(lat_points(num_lats))
      do ic = 1, num_lats
         lat_points(ic) = (ic-1)*dy + min_lat
      enddo

      merra2_grid = create_rectilinear_grid(lon_points, lat_points, &
                               "MERRA2 Grid", LIS_rc%npesx, LIS_rc%npesy,  &
                                merra2_struc(n)%coordSys, &
                                periodic   = .TRUE., &
                                staggerloc = merra2_struc(n)%staggerloc)

      DEALLOCATE(lon_points, lat_points)

      ! Add fields to the bundle
      DO ic = 1, num_merra2_fields
         call addTracerToBundle(merra2_struc(n)%forcing_bundle, merra2_grid, &
                   TRIM(list_merra2_fields(ic)), merra2_struc(n)%type_kind, dummy_array)
      ENDDO

      end subroutine create_merra2_Forcing_ESMFbundle
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_merra2_Model_ESMFbundle(n)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod !,       only : create_regular_grid
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
      integer          :: rc, ic, num_lats, num_lons
      real             :: model_gridDesc(10)
      type(ESMF_Grid)  :: model_grid
      real             :: dummy_array(1,1) = 0.0
      real, pointer    :: lat_points(:)
      real, pointer    :: lon_points(:)
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize model ESMF objects.'

      ! Create the bundles
      write(num_st, '(i2.2)') n

      ! ---> Bundle for the model
      LIS_domain(n)%merra2_bundle = ESMF_FieldBundleCreate(name = TRIM(merra2_bundle_bname)//num_st, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleCreate failed for model Data')

      ! ---> ESMF grid for for the model
      num_lons = LIS_rc%gnc(n)
      num_lats = LIS_rc%gnr(n)

      ALLOCATE(lon_points(num_lons))
      lon_points(:) = LIS_domain(n)%glon(1:num_lons)

      ALLOCATE(lat_points(num_lats))
      lat_points(:) = LIS_domain(n)%glat(1:(num_lats-1)*num_lons:num_lons)

      model_grid = create_rectilinear_grid(lon_points, lat_points, &
                               "Model Grid", LIS_rc%npesx, LIS_rc%npesy, &
                                LIS_domain(n)%coordSys, &
                                staggerloc = LIS_domain(n)%staggerloc)

      DEALLOCATE(lon_points, lat_points)

      ! Add fields to the bundle
      DO ic = 1, num_merra2_fields
         call addTracerToBundle(LIS_domain(n)%merra2_bundle, model_grid, &
                   TRIM(list_merra2_fields(ic)), merra2_struc(n)%type_kind, dummy_array)
      ENDDO

      end subroutine create_merra2_Model_ESMFbundle
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_merra2_ESMFroutehandle(n)
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
! the MERRA2 forcing and the model.
!
! !LOCAL VARIABLES:
      type(ESMF_FIELD) :: model_field, merra2_field
      integer          :: ftc(2), ftlb(2), ftub(2)
      real, pointer    :: farray2dd(:,:)
      integer          :: rc
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Determine the ESMF routehandle.'

      ! Get one field from the MERRA2 focing bundle
      call ESMF_FieldBundleGet(merra2_struc(n)%forcing_bundle, TRIM(list_merra2_fields(1)), &
                              field = merra2_field, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for merra2 field'//TRIM(list_merra2_fields(1)))

      ! Get the corresponding field from the model bundle
      call ESMF_FieldBundleGet(LIS_domain(n)%merra2_bundle, TRIM(list_merra2_fields(1)), &
                               field = model_field, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for model field'//TRIM(list_merra2_fields(1)))


      ! Compute the ESMF routehandle
      call createESMF_RouteHandle(merra2_field, model_field,  merra2_struc(n)%regridMethod, &
                                  merra2_struc(n)%undefined_value, merra2_struc(n)%routehandle, &
                                  merra2_struc(n)%dynamicMask, rc)
      call LIS_verify(rc, 'createESMF_RouteHandle failed')

      end subroutine create_merra2_ESMFroutehandle
!EOC
!------------------------------------------------------------------------------
end module merra2_forcingMod
