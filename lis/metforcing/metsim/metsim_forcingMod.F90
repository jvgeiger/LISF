!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module metsim_forcingMod
!BOP
! !MODULE: metsim_forcingMod
! 
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the 150-yr MetSim data in netcdf format.
!  The data is in latlon projection, and at 3 hourly intervals. 
!  The implemenatation in LIS has the derived data type {\tt metsim\_struc} 
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[fmodeltime1]
!    The nearest, previous 6 hour instance of the incoming 
!    data (as a real time). 
!  \item[fmodeltime2]
!    The nearest, next 6 hour instance of the incoming 
!    data (as a real time).
!  \item[metsimdir]
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
!  \item[n12,n122,n212,n222]
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
! 26 Jan 2007: Hiroko Beaudoing; Initial Specification adopted from 
!                                LIS/retberg.F90
! 16 Feb 2016: Hiroko Beaudoing; Fixed indexes for precip and pressure fields 
!                                so budget-bilinear interpolation applies to 
!                                correct individual fields.
! 15 May 2017: Bailing Li; Added changes for reading in version 2.2 data
!                          that is in 3D array (4D in version 1 & 2).
! 22 Oct 2018: Daniel Sarmiento; Added changes to support version 3 data
! 29 Apr 2020: Zhuo Wang: Added MetSim forcing for LIS-SUMMA 
!
! !USES: 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_MetSim      !defines the native resolution of 
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: metsim_struc

!EOP

  type, public :: metsim_type_dec
     real                   :: ts
     integer                :: ncold, nrold   
     character*100          :: metsimdir    ! MetSim Forcing Directory 
     integer                :: mi
     real*8                 :: metsimtime1,metsimtime2

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

     integer           :: nIter, st_iterid,en_iterid  ! Forecast parameters  ??

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type metsim_type_dec

  type(metsim_type_dec), allocatable :: metsim_struc(:) 

contains
  
!BOP
!
! !ROUTINE: init_MetSim
! \label{init_MetSim}
!
! !REVISION HISTORY: 
! 26 Jan 2007: Hiroko Kato; Initial Specification
! 15 May 2017: Bailing Li; Added changes for reading in version 2.2 data
! 22 Oct 2018: Daniel Sarmiento; Added changes to support version 3 data
! 29 Apr 2020: Zhuo Wang: Added MetSim forcing for LIS-SUMMA.
! 
! !INTERFACE:
  subroutine init_MetSim(findex)

! !USES: 
   use LIS_coreMod
   use LIS_timeMgrMod
   use LIS_logMod
   use LIS_forecastMod

   implicit none
   integer, intent(in)  :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for MetSim 
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_metsim](\ref{readcrd_metsim}) \newline
!     reads the runtime options specified for MetSim data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \end{description}
!EOP
    real     :: gridDesci(50)
    integer  :: n 

    ! ??? Added according to merra-land ???
    integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt

    gridDesci = 0.0 
    allocate(metsim_struc(LIS_rc%nnest))
    call readcrd_metsim()

    do n=1, LIS_rc%nnest
       metsim_struc(n)%ts = 3*3600
       call LIS_update_timestep(LIS_rc, n, metsim_struc(n)%ts)
    enddo

    ! ??? 7 forcing variables ??
    LIS_rc%met_nf(findex) = 7 

    ! Set MetSim grid dimensions and extent information:
    ! ???
    metsim_struc(:)%ncold = 464
    metsim_struc(:)%nrold = 224

    do n=1,LIS_rc%nnest

       ! Regular retrospective or non-forecast mode:
       allocate(metsim_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(metsim_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       metsim_struc(n)%st_iterid = 1
       metsim_struc(n)%en_iterId = 1
       metsim_struc(n)%nIter = 1


       metsim_struc(n)%metdata1 = 0
       metsim_struc(n)%metdata2 = 0

       ! Added based on NLDAS2 data
      !metsim_struc(n)%gridDesci = 0
      !metsim_struc(n)%findtime1 = 0
      !metsim_struc(n)%findtime2 = 0

      gridDesci(1) = 0
      gridDesci(2) = metsim_struc(n)%ncold
      gridDesci(3) = metsim_struc(n)%nrold
      !Define driver data domains
      gridDesci(4) = 25.0625      ! min(lat)
      gridDesci(5) = -124.9375    ! min(lon)
      gridDesci(6) = 128
      gridDesci(7) = 52.9375      ! max(lat) 
      gridDesci(8) = -67.0625     ! max(lon) 
      gridDesci(9) = 0.125        ! dlat
      gridDesci(10) = 0.125       ! dlon

      !???? For merra-land and merra, gridDesci(n,20) = 0
      gridDesci(20) = 64

!     ! Check for grid and interp option selected:
!     if( gridDesci(9)  == LIS_rc%gridDesc(n,9) .and. &
!         gridDesci(10) == LIS_rc%gridDesc(n,10).and. &
!         LIS_rc%gridDesc(n,1) == proj_latlon .and. &
!         LIS_rc%met_interp(findex) .ne. "neighbor" ) then
!       write(LIS_logunit,*) "[ERR] The MetSim grid was selected for the"
!       write(LIS_logunit,*) "[ERR] LIS run domain; however, 'bilinear', 'budget-bilinear',"
!       write(LIS_logunit,*) "[ERR] or some other unknown option was selected to spatially"
!       write(LIS_logunit,*) "[ERR] downscale the grid, which will cause errors during runtime."
!       write(LIS_logunit,*) "[ERR] Program stopping ..."
!       call LIS_endrun()
!     endif

      metsim_struc(n)%mi = metsim_struc(n)%ncold*metsim_struc(n)%nrold

      ! Setting up weights for spatial Interpolation
      if(trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
         allocate(metsim_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,&
              metsim_struc(n)%n111,metsim_struc(n)%n121,&
              metsim_struc(n)%n211,metsim_struc(n)%n221,&
              metsim_struc(n)%w111,metsim_struc(n)%w121,&
              metsim_struc(n)%w211,metsim_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") then
         allocate(metsim_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(metsim_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,&
              metsim_struc(n)%n111,metsim_struc(n)%n121,&
              metsim_struc(n)%n211,metsim_struc(n)%n221,&
              metsim_struc(n)%w111,metsim_struc(n)%w121,&
              metsim_struc(n)%w211,metsim_struc(n)%w221)

         allocate(metsim_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(metsim_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(metsim_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(metsim_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(metsim_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(metsim_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(metsim_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(metsim_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               metsim_struc(n)%n112,metsim_struc(n)%n122,&
               metsim_struc(n)%n212,metsim_struc(n)%n222,&
               metsim_struc(n)%w112,metsim_struc(n)%w122,&
               metsim_struc(n)%w212,metsim_struc(n)%w222)
        elseif(trim(LIS_rc%met_interp(findex)) .eq. "neighbor") then 
          allocate(metsim_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
        
          call neighbor_interp_input(n,gridDesci,&
               metsim_struc(n)%n113)
       else
         write(LIS_logunit,*) '[ERR] Interpolation option '// &
               trim(LIS_rc%met_interp(findex))//&
               ' for MetSim forcing is not supported'
         call LIS_endrun
       endif
!-------------------------------------------------------------
    enddo

  end subroutine init_MetSim
end module metsim_forcingMod
