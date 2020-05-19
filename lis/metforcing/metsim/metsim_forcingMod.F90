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
!  \item[nmif]
!    Number of forcing variables in the MetSim data
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
     integer                :: nmif
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
     integer                :: nIter, st_iterid, en_iterid ! Forecast parameters

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
! 29 Apr 2020: Zhuo Wang: Added MetSim forcing for LIS-SUMMA.
! 
! !INTERFACE:
  subroutine init_MetSim(findex)

! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_logMod,     only : LIS_logunit, LIS_endrun
   use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
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

    integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt

    allocate(metsim_struc(LIS_rc%nnest))

    do n=1, LIS_rc%nnest
       metsim_struc(n)%ts = 3*3600
       call LIS_update_timestep(LIS_rc, n, metsim_struc(n)%ts)
    enddo

    call readcrd_metsim()

    metsim_struc(:)%nmif  = 7
    LIS_rc%met_nf(findex) = 7 

    metsim_struc(:)%ncold = 464
    metsim_struc(:)%nrold = 224

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

          metsim_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          metsim_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          metsim_struc(n)%nIter = LIS_forecast_struc(1)%niterations

          allocate(metsim_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
             LIS_rc%met_nf(findex),&
             LIS_rc%ngrid(n)))
          allocate(metsim_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
             LIS_rc%met_nf(findex),&
             LIS_rc%ngrid(n)))

       ! Regular retrospective or non-forecast mode:
       else
          metsim_struc(n)%st_iterid = 1
          metsim_struc(n)%en_iterId = 1
          metsim_struc(n)%nIter = 1

          allocate(metsim_struc(n)%metdata1(1, LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(metsim_struc(n)%metdata2(1, LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

       endif

      metsim_struc(n)%metdata1 = 0
      metsim_struc(n)%metdata2 = 0
      metsim_struc(n)%findtime1 = 0
      metsim_struc(n)%findtime2 = 0

      gridDesci = 0.0
      gridDesci(1) = 0
      gridDesci(2) = metsim_struc(n)%ncold
      gridDesci(3) = metsim_struc(n)%nrold
      gridDesci(4) = 25.0625      ! min(lat)
      gridDesci(5) = -124.9375    ! min(lon)
      gridDesci(6) = 128
      gridDesci(7) = 52.9375      ! max(lat) 
      gridDesci(8) = -67.0625     ! max(lon) 
      gridDesci(9) = 0.125        ! dlat
      gridDesci(10) = 0.125       ! dlon
      gridDesci(20) = 64

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
