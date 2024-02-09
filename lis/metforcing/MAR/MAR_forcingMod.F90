!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module MAR_forcingMod
!BOP
! !MODULE: MAR_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data extracted from MAR output
!  files.  Here MAR output files are consider input data.
!
!  The implementation in LIS has the derived data type {\tt MAR\_struc} that
!  includes the variables that specify the runtime options.
!  They are described below: 
!  \begin{description}
!  \item[MARdir]
!    Directory containing the input data
!  \item[ts]
!    Frequency in seconds of the forcing data
!  \item[MARtime1]
!    The nearest, previous 1 hour instance of the incoming 
!    data (as a real time). 
!  \item[MARtime2]
!    The nearest, next 1 hour instance of the incoming 
!    data (as a real time).
!  \item[findtime1, findtime2]
!    boolean flags to indicate which time is to be read for 
!    temporal interpolation.
!  \end{description}
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
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

! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_MAR      !defines the native resolution of 
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: MAR_struc
!EOP

  type, public    :: MAR_type_dec
     real         :: ts
     integer      :: nc, nr
     integer      :: nest_id
     character(len=LIS_CONST_PATH_LEN) :: MARdir
     real*8       :: MARtime1,MARtime2
     integer      :: findtime1,findtime2
     integer      :: mi

     real         :: gridDesci(50)

     integer, allocatable  :: n111(:)
     integer, allocatable  :: n121(:)
     integer, allocatable  :: n211(:)
     integer, allocatable  :: n221(:)
     real, allocatable     :: w111(:),w121(:)
     real, allocatable     :: w211(:),w221(:)
     integer, allocatable  :: n112(:,:)
     integer, allocatable  :: n122(:,:)
     integer, allocatable  :: n212(:,:)
     integer, allocatable  :: n222(:,:)
     real, allocatable     :: w112(:,:),w122(:,:)
     real, allocatable     :: w212(:,:),w222(:,:)

     ! Neighbor
     integer, allocatable  :: n113(:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
     !integer      :: nIter, st_iterid, en_iterid

  end type MAR_type_dec

  type(MAR_type_dec), allocatable :: MAR_struc(:)
!EOP
contains
  
!BOP
!
! !ROUTINE: init_MAR
! \label{init_MAR}
!
! !REVISION HISTORY: 
! 18 Aug 2023: Mahdi Navari; Initial Specification
! 
! !INTERFACE:
  subroutine init_MAR(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain, LIS_localPet,&
        LIS_ews_ind, LIS_ewe_ind,&
        LIS_nss_ind, LIS_nse_ind,&
        LIS_ews_halo_ind,LIS_ewe_halo_ind,&
        LIS_nss_halo_ind, LIS_nse_halo_ind
  use LIS_timeMgrMod, only : LIS_update_timestep, LIS_calendar
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, LIS_verify
  use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
  use LIS_forecastMod
  use map_utils

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for LIS output
!  data. 
!  Note that no interpolation is performed for
!  this forcing. The data is expected to be in the same map projection
!  and resolution as that of the current LIS run.
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_MAR](\ref{readcrd_MAR}) \newline
!     reads the runtime options specified for MAR output data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    
    integer :: n
    !real  :: gridDesci(LIS_rc%nnest, 50)

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the MAR forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    
    allocate(MAR_struc(LIS_rc%nnest))

    call readcrd_MAR()

    do n=1, LIS_rc%nnest
       MAR_struc(n)%ts = 60*60
       call LIS_update_timestep(LIS_rc, n, MAR_struc(n)%ts)
    enddo

    !gridDesci = 0
    LIS_rc%met_nf(findex) = 9
 
! Number of x- and y-dir points:
    MAR_struc(:)%nc = 76 !60
    MAR_struc(:)%nr =136 ! 110

    do n=1, LIS_rc%nnest

       allocate(MAR_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(MAR_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       MAR_struc(n)%metdata1 = 0
       MAR_struc(n)%metdata2 = 0


      ! See below link for more information on polar stereographic grids:
      ! https://epsg.io/3413
      ! PROJECTION["Polar_Stereographic"],
      ! PARAMETER["latitude_of_origin",70],
      ! PARAMETER["central_meridian",-45],
      ! PARAMETER["false_easting",0],
      ! PARAMETER["false_northing",0],

       MAR_struc(n)%gridDesci = 0
       MAR_struc(n)%gridDesci(1) = 5         ! Polar-stereographic
       MAR_struc(n)%gridDesci(2) = MAR_struc(n)%nc !number of columns in the domain
       MAR_struc(n)%gridDesci(3) = MAR_struc(n)%nr !number of rows in the domain
       MAR_struc(n)%gridDesci(4) = 60        ! latitude of origin -- LL Lat
       MAR_struc(n)%gridDesci(5) = -180      ! longitude of origin -- LL Lon
       MAR_struc(n)%gridDesci(6) = 8         ! not used
       MAR_struc(n)%gridDesci(7) = -45       ! orientation of the grid
       MAR_struc(n)%gridDesci(8) = 20        ! grid spacing in km
       MAR_struc(n)%gridDesci(9) = 20        ! grid spacing in km
       MAR_struc(n)%gridDesci(10) = 70       ! true lat1
       MAR_struc(n)%gridDesci(11) = -45      ! standard long
       !MAR_struc(n)%gridDesci(13) =          ! regional domain
       MAR_struc(n)%gridDesci(20) = 255      !used to specify the ordering of data (non-divisible by 32 indicates E-W ordering else N-S ordering)

       MAR_struc(n)%mi = MAR_struc(n)%nc*MAR_struc(n)%nr

! The data is expected to be in the same map projection
! and resolution as the current LIS run.

! TODO: We can test and activate the following options in the future.
#if 0
! NOTE: Bilinear interpolation (bilinear_interp_input.F90), 
! conservative interpolation (conserv_interp_input.F90), 
! and neighbor interpolation (neighbor_interp_input.F90) issue calls to get_fieldpos.F90. However, 
! please be aware that this subroutine does not support polar stereographic projection.


     ! Setting up weights for Interpolation
       select case( LIS_rc%met_interp(findex) )

         case( "bilinear" )
          allocate(MAR_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,MAR_struc(n)%gridDesci,&
               MAR_struc(n)%n111,MAR_struc(n)%n121,&
               MAR_struc(n)%n211,MAR_struc(n)%n221,&
               MAR_struc(n)%w111,MAR_struc(n)%w121,&
               MAR_struc(n)%w211,MAR_struc(n)%w221)

         case( "budget-bilinear" )
          allocate(MAR_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(MAR_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,MAR_struc(n)%gridDesci,&
               MAR_struc(n)%n111,MAR_struc(n)%n121,&
               MAR_struc(n)%n211,MAR_struc(n)%n221,&
               MAR_struc(n)%w111,MAR_struc(n)%w121,&
               MAR_struc(n)%w211,MAR_struc(n)%w221)

          allocate(MAR_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(MAR_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(MAR_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(MAR_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(MAR_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(MAR_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(MAR_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(MAR_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,MAR_struc(n)%gridDesci,&
               MAR_struc(n)%n112,MAR_struc(n)%n122,&
               MAR_struc(n)%n212,MAR_struc(n)%n222,&
               MAR_struc(n)%w112,MAR_struc(n)%w122,&
               MAR_struc(n)%w212,MAR_struc(n)%w222)

         case( "neighbor" )

          allocate(MAR_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,MAR_struc(n)%gridDesci,&
                        MAR_struc(n)%n113)

       case default
         write(LIS_logunit,*) "[ERR] User-input issue with MAR forcing ..."
         write(LIS_logunit,*) " -- Currently only supported interpolation options include:"
         write(LIS_logunit,*) "  - bilinear, budget-bilinear, or neighbor - "
         call LIS_endrun
       end select
#endif 
    enddo
    ! MAR
    ! TTZ (C)      tair     <--> LIS_forc%metdata1(1,:)
    ! QQZ (g/kg)   qai      <--> LIS_forc%metdata1(2,:)
    ! SWD (W/m2)  swd     <--> LIS_forc%metdata1(3,:)
    ! LWD (W/m2)  lwd     <--> LIS_forc%metdata1(4,:)
    ! UUZ  (m/s)   uwind   <--> LIS_forc%metdata1(5,:)
    ! VVZ  (m/s)   vwind   <--> LIS_forc%metdata1(6,:)
    ! SP  (hpa)   psurf   <--> LIS_forc%metdata1(7,:)
    ! RF  (mmWE/day) (railf)  <--> LIS_forc%metdata1(8,:)    
    ! SF  (mmWE/day) (snowf)  <--> LIS_forc%metdata1(9,:)

  end subroutine init_MAR
end module MAR_forcingMod
