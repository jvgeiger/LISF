!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdasT1534_forcingMod
!BOP
! !MODULE: gdasT1534_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Data
!  Assimilation System (GDAST1534) developed at the Environmental Modeling
!  Center (EMC) of NOAA/NCEP. GDAST1534 forcing variables are produced
!  on a quadratic gaussian grid. 
!  The implementation in LIS has the derived data type {\tt gdasT1534\_struc} that
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the GDAST1534 data
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \item[gdasT1534dir]
!    Directory containing the input data
!  \item[elevfile]
!    File with the elevation definition for the input data. 
!  \item[mi]
!    Number of points in the input grid
!  \item[n11,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES: 
      use ESMF

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDAST1534      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdasT1534_struc
  public :: num_gdasT1534_fields
  public :: list_gdasT1534_fields
!EOP

  type, public :: gdasT1534_type_dec
     real          :: ts
     integer       :: ncold, nrold   !AWIPS 212 dimensions
     integer       :: nmif
     character*100 :: gdasT1534dir   !GDAST1534 Forcing Directory
     real*8        :: gdasT1534time1,gdasT1534time2
     integer       :: mi
     integer       :: findtime1, findtime2

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
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

     ! For ESMF regridding
     type(ESMF_FieldBundle)       :: forcing_bundle
     type(ESMF_RouteHandle)       :: routehandle
     type(ESMF_DynamicMask)       :: dynamicMask
     type(ESMF_TypeKind_Flag)     :: type_kind = ESMF_TYPEKIND_R4
     type(ESMF_STAGGERLOC)        :: staggerloc
     type(ESMF_RegridMethod_Flag) :: regridMethod
     real                         :: undefined_value ! for missing value

  end type gdasT1534_type_dec
  
  type(gdasT1534_type_dec), allocatable :: gdasT1534_struc(:)

  integer            :: num_gdasT1534_fields       ! number of available fields
  character(len=100) :: list_gdasT1534_fields(30)  ! list of name of fields
  character(len=30), parameter :: gdasT1534_bundle_bname = "gdasT1534_bundle_"

contains
  
!BOP
!
! !ROUTINE: init_GDAST1534
!  \label{init_GDAST1534}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_GDAST1534(findex)
! !USES: 
    use LIS_coreMod !,    only : LIS_rc, LIS_domain
    use LIS_logMod !,     only : LIS_logunit, LIS_endrun
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep

      use LIS_FORC_AttributesMod
    
    implicit none
! !ARGUMENTS: 
    integer,  intent(in)    :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GDAST1534
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes \ref{interp}. Based on the GDAST1534 data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights. The dates of the GDAST1534 resolution switches are also 
!  defined in this routine. 
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gdasT1534](\ref{readcrd_gdasT1534}) \newline
!     reads the runtime options specified for GDAST1534 data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n

    integer          :: ic, rc

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GDAS-T1534 forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(gdasT1534_struc(LIS_rc%nnest))
    
    do n=1, LIS_rc%nnest
       gdasT1534_struc(n)%ts = 21600
       call LIS_update_timestep(LIS_rc, n, gdasT1534_struc(n)%ts)
    enddo

    call readcrd_gdasT1534()
    
    gdasT1534_struc(:)%nmif    = 16 
    LIS_rc%met_nf(findex) = 16 !number of met variables in GDAST1534 forcing

    IF (LIS_rc%do_esmfRegridding) THEN
       CALL set_list_gdasT1534_fields()
    ENDIF
    
    do n=1,LIS_rc%nnest

       allocate(gdasT1534_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gdasT1534_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gdasT1534_struc(n)%metdata1 = 0
       gdasT1534_struc(n)%metdata2 = 0

       gdasT1534_struc(n)%findtime1 = 0 
       gdasT1534_struc(n)%findtime2 = 0 

       ! This grid is good for some time in the 1980's.
       ! Look up the exact dates.
       gridDesci = 0 
       gridDesci(1) = 4
       gridDesci(2) = 3072
       gridDesci(3) = 1536
       gridDesci(4) = 89.91153
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) =  -89.91153
       gridDesci(8) = -0.117187
       gridDesci(9) = 0.117187
       gridDesci(10) = 768.0
       gridDesci(20) = 0
       gdasT1534_struc(n)%ncold = 3072
       gdasT1534_struc(n)%nrold = 1536
       gdasT1534_struc(n)%mi = gdasT1534_struc(n)%ncold*gdasT1534_struc(n)%nrold
       
       IF (LIS_rc%do_esmfRegridding) THEN
          !gdasT1534_struc(n)%type_kind       = ESMF_TYPEKIND_R4
          gdasT1534_struc(n)%undefined_value = 9999.0

          gdasT1534_struc(n)%regridMethod = ESMF_REGRIDMETHOD_NEAREST_STOD ! ESMF_REGRIDMETHOD_BILINEAR 
          gdasT1534_struc(n)%staggerloc = ESMF_STAGGERLOC_CENTER
          LIS_domain(n)%staggerloc      = ESMF_STAGGERLOC_CENTER

          CALL create_gdasT1534Forcing_ESMFbundle(n, gridDesci(1:10))
          CALL create_gdasT1534Model_ESMFbundle(n)
          CALL create_gdasT1534_ESMFroutehandle(n)
       ELSE
          !Setting up weights for Interpolation
          if(LIS_rc%met_interp(findex).eq."bilinear") then 

             allocate(gdasT1534_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call bilinear_interp_input(n,gridDesci,&
                  gdasT1534_struc(n)%n111,gdasT1534_struc(n)%n121,&
                  gdasT1534_struc(n)%n211,gdasT1534_struc(n)%n221,&
                  gdasT1534_struc(n)%w111,gdasT1534_struc(n)%w121,&
                  gdasT1534_struc(n)%w211,gdasT1534_struc(n)%w221)
          elseif(LIS_rc%met_interp(findex).eq."budget-bilinear") then 
             allocate(gdasT1534_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call bilinear_interp_input(n,gridDesci,&
                  gdasT1534_struc(n)%n111,gdasT1534_struc(n)%n121,&
                  gdasT1534_struc(n)%n211,gdasT1534_struc(n)%n221,&
                  gdasT1534_struc(n)%w111,gdasT1534_struc(n)%w121,&
                  gdasT1534_struc(n)%w211,gdasT1534_struc(n)%w221)
             allocate(gdasT1534_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             call conserv_interp_input(n,gridDesci,&
                  gdasT1534_struc(n)%n112,gdasT1534_struc(n)%n122,&
                  gdasT1534_struc(n)%n212,gdasT1534_struc(n)%n222,&
                  gdasT1534_struc(n)%w112,gdasT1534_struc(n)%w122,&
                  gdasT1534_struc(n)%w212,gdasT1534_struc(n)%w222)
          elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
             allocate(gdasT1534_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call neighbor_interp_input(n,gridDesci,&
                  gdasT1534_struc(n)%n113)
          else
             write(LIS_logunit,*) 'The specified spatial interpolation option not'
             write(LIS_logunit,*) 'supported for GDAST1534..'
             write(LIS_logunit,*) 'program stopping..'
             call LIS_endrun()
          endif

       ENDIF

    enddo
  end subroutine init_GDAST1534

!------------------------------------------------------------------------------
!BOP
      subroutine set_list_gdasT1534_fields()
!
! !USES:
      use LIS_FORC_AttributesMod
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
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
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Tair%varname(1))
      endif
      if (LIS_FORC_Qair%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Qair%varname(1))
      endif
      if (LIS_FORC_SWdown%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_SWdown%varname(1))
      endif
      if (LIS_FORC_LWdown%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_LWdown%varname(1))
      endif
      if (LIS_FORC_Wind_E%selectOpt .eq. 1) then
         ic = ic + 1
          list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Wind_E%varname(1))
        endif
      if (LIS_FORC_Wind_N%selectOpt .eq. 1) then
        ic = ic + 1
        list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Wind_N%varname(1))
      endif
      if (LIS_FORC_Psurf%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Psurf%varname(1))
      endif
      if (LIS_FORC_Rainf%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Rainf%varname(1))
      endif

      if (LIS_FORC_CRainf%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_CRainf%varname(1))
      endif

      if (LIS_FORC_Forc_Hgt%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Forc_Hgt%varname(1))
       endif
      if (LIS_FORC_Ch%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Ch%varname(1))
       endif
      if (LIS_FORC_Z0%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Z0%varname(1))
      endif
      if (LIS_FORC_GVF%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_GVF%varname(1))
      endif
      if (LIS_FORC_Alb%selectOpt .eq. 1) then
         ic = ic + 1
         list_gdasT1534_fields(ic) = TRIM(LIS_FORC_Alb%varname(1))
      endif

      num_gdasT1534_fields = ic

      end subroutine set_list_gdasT1534_fields
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_gdasT1534Forcing_ESMFbundle(n, forcing_gridDesc)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod,      only :create_gaussian_grid
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
      real,    intent(in) :: forcing_gridDesc(10)
!
! !DESCRIPTION:
! Create the GDAS forcing ESMF grid and bundle.
! This subroutine might be called several times depending on the integration date.
! Howver, it will be called once for any period when the GDAS forcing resolution
! does not change.
!
! !LOCAL VARIABLES:
      character(len=2) :: num_st
      integer          :: rc, ic
      type(ESMF_Grid)  :: gdasT1534_grid
      real             :: dummy_array(1,1) = 0.0
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize GDAS T1534 ESMF object.'
      ! Create the bundles
      write(num_st, '(i2.2)') n

      ! ---> Bundle for gdasT1534
      gdasT1534_struc(n)%forcing_bundle = ESMF_FieldBundleCreate(name = TRIM(gdasT1534_bundle_bname)//num_st, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleCreate failed for gdasT1534 Forcing Data')

      gdasT1534_grid =  create_gaussian_grid(forcing_gridDesc, "gdasT1534 Grid", &
                                        LIS_rc%npesx, LIS_rc%npesy, &
                                        staggerloc = gdasT1534_struc(n)%staggerloc, &
                                        global_domain = .TRUE.)
      ! Add fields to the bundle
      DO ic = 1, num_gdasT1534_fields
         if (LIS_masterproc) PRINT*,"--->Adding-Forcing: ", ic, TRIM(list_gdasT1534_fields(ic))
         call addTracerToBundle(gdasT1534_struc(n)%forcing_bundle, gdasT1534_grid, &
                   TRIM(list_gdasT1534_fields(ic)), gdasT1534_struc(n)%type_kind, dummy_array)
      ENDDO

      end subroutine create_gdasT1534Forcing_ESMFbundle
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_gdasT1534Model_ESMFbundle(n)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod,       only : create_gaussian_grid
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
      LIS_domain(n)%gdasT1534_bundle = ESMF_FieldBundleCreate(name = TRIM(gdasT1534_bundle_bname)//num_st, rc=rc)
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

     !PRINT"(i4,8f10.4)",LIS_localPet, &
     !                    model_gridDesc(2),model_gridDesc(3),model_gridDesc(4), &
     !                    model_gridDesc(5),model_gridDesc(7),model_gridDesc(8), &
     !                    model_gridDesc(9),model_gridDesc(10)

      model_grid =  create_gaussian_grid(model_gridDesc, "Model Grid", &
                                        LIS_rc%npesx, LIS_rc%npesy, &
                                        staggerloc = LIS_domain(n)%staggerloc, &
                                        global_domain = .TRUE.)

      ! Add fields to the bundle
      DO ic = 1, num_gdasT1534_fields
         if (LIS_masterproc) PRINT*,"<---Adding-Model: ", ic, TRIM(list_gdasT1534_fields(ic))
         call addTracerToBundle(LIS_domain(n)%gdasT1534_bundle, model_grid, &
                   TRIM(list_gdasT1534_fields(ic)), gdasT1534_struc(n)%type_kind, dummy_array)
      ENDDO

      end subroutine create_gdasT1534Model_ESMFbundle
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_gdasT1534_ESMFroutehandle(n)
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
! the GDAS forcing and the model.
!
! !LOCAL VARIABLES:
      type(ESMF_FIELD) :: model_field, gdasT1534_field
      integer          :: ftc(2), ftlb(2), ftub(2)
      real, pointer    :: farray2dd(:,:)
      integer          :: rc
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Determine the ESMF routehandle.'

      ! Get one field from the GDAS focing bundle
      call ESMF_FieldBundleGet(gdasT1534_struc(n)%forcing_bundle, TRIM(list_gdasT1534_fields(1)), &
                              field = gdasT1534_field, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for gdasT1534 field'//TRIM(list_gdasT1534_fields(1)))

!    call ESMF_FieldGet(gdasT1534_field, localDe=0, farrayPtr=farray2dd, &
!        totalLBound=ftlb, totalUBound=ftub, totalCount=ftc, rc=rc)
!        print"(a14,7i7)", "ForcingGrid: ", LIS_localPet,ftlb,ftub,ftc
!      call LIS_verify(rc, 'ESMF_FieldGet failed for GDAS T1534 field'//TRIM(list_gdasT1534_fields(1)))


      ! Get the corresponding field from the model bundle
      call ESMF_FieldBundleGet(LIS_domain(n)%gdasT1534_bundle, TRIM(list_gdasT1534_fields(1)), &
                               field = model_field, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for model field'//TRIM(list_gdasT1534_fields(1)))


!    call ESMF_FieldGet(model_field, localDe=0, farrayPtr=farray2dd, &
!        totalLBound=ftlb, totalUBound=ftub, totalCount=ftc, rc=rc)
!        print"(a14,7i7)", "ModelGrid: ", LIS_localPet,ftlb,ftub,ftc
!      call LIS_verify(rc, 'ESMF_FieldGet failed for model field'//TRIM(list_gdasT1534_fields(1)))


      ! Compute the ESMF routehandle
      call createESMF_RouteHandle(gdasT1534_field, model_field,  gdasT1534_struc(n)%regridMethod, &
                                  gdasT1534_struc(n)%undefined_value, gdasT1534_struc(n)%routehandle, &
                                  gdasT1534_struc(n)%dynamicMask, rc)
      call LIS_verify(rc, 'createESMF_RouteHandle failed')

      end subroutine create_gdasT1534_ESMFroutehandle
!EOC
!------------------------------------------------------------------------------

end module gdasT1534_forcingMod

