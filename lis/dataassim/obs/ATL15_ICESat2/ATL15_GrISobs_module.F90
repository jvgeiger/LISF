!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
#include "LIS_NetCDF_inc.h"
!BOP
!
! !MODULE: ATL15_GrISobs_module
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle ICESat2-ATL15 observations for Greenland (GrIS)
!   ATL15 provides coarser resolution (1 km, 10 km, 20 km, and 40 km) height-change maps at 3 month intervals.
!   Temporal Coverage: start at 29 March 2019 to present
!   Temporal Resolution: 3 month
!   Spatial Reference System(s): WGS 84 / NSIDC Sea Ice Polar Stereographic NorthEPSG:3413 (Greenland)
!   Spatial Coverage: N:90S:59E:180W:-180 (Greenland)  

! Based on data, the temporal coverage are as follows:
!   First and second observations for fo dh are @ 2018-10-01 22:30 and 2019-01-01 06:00
!   First observation for dhdt is @ 2018-11-16 14:15 ( in the middle of the 1st and 2st dh obervation times)   
!   Time interval is 91.3125 days or 7889400 seconds

! !REVISION HISTORY: 
!  13 Dec 2023    Mahdi Navari;   Initial Specification
! 
module ATL15_GrISobs_module
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ATL15_GrISobs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ATL15_GrIS_struc
!EOP
  type, public :: ATL15_GrIS_dec
     type(ESMF_TimeInterval) :: ts
     !integer           :: mo
     !real              :: alarmStartdate
     integer           :: useDistErr
     logical           :: startMode
     integer           :: nc, nr, num_obs
     real,allocatable  :: obs_time 
  end type ATL15_GrIS_dec
  
  type(ATL15_GrIS_dec),allocatable :: ATL15_GrIS_struc(:)
  
contains
!BOP
! 
! !ROUTINE: ATL15_GrISobs_setup
! \label{ATL15_GrISobs_setup}
! 
! !INTERFACE: 
  subroutine ATL15_GrISobs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_DAobservationsMod
    use LIS_perturbMod

    implicit none 

! !ARGUMENTS: 
    integer                ::  k 
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for ATL15_GrIS (elevation change)
!   assimilation
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n 
    integer                ::  ftn
    integer                ::  i
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    type(ESMF_Time)        ::  startTime
    character(len=LIS_CONST_PATH_LEN) ::  ATL15_GrISobsdir
    character*100          ::  temp
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    real,      allocatable     ::  ssdev(:)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)

    allocate(ATL15_GrIS_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: ATL15_GrISobs_setup')

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: ATL15_GrISobs_setup')

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: ATL15_GrISobs_setup')
    
    call ESMF_ConfigFindLabel(LIS_config,"ATL15_GrIS data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ATL15_GrISobsdir,&
            rc=status)
       call LIS_verify(status,'ATL15_GrIS data directory: not defined')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            ATL15_GrISobsdir, rc=status)
       call LIS_verify(status)
    enddo

!    call ESMF_ConfigFindLabel(LIS_config,"ATL15_GrIS use reported measurement error values:",&
!         rc=status)
!    do n=1,LIS_rc%nnest
!       call ESMF_ConfigGetAttribute(LIS_config,ATL15_GrIS_struc(n)%useDistErr,&
!            rc=status)
!       call LIS_verify(status,'ATL15_GrIS use reported measurement error values: not defined')
!    enddo


    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(OBS_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Update Time",&
            -99.0, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
            LIS_rc%obs_ngrid(k),rc=status)
       call LIS_verify(status)
       
    enddo

    !do n=1,LIS_rc%nnest
    !   call ESMF_TimeIntervalSet(ATL15_GrIS_struc(n)%ts,s=7889400,rc=status)
    !   call LIS_verify(status, 'Error in ESMF_TimeIntervalSet in ATL15 GrIS')
    !enddo

    write(LIS_logunit,*)'[INFO] read ATL15_GrIS data specifications'

!----------------------------------------------------------------------------
!   Create the array containers that will contain the observations and
!   the perturbations. For the ATL15_GrIS case, it is assumed that the 
!   observations are in the grid space. 
!   
!----------------------------------------------------------------------------

    do n=1,LIS_rc%nnest

       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid
       
       obsField(n) = ESMF_FieldCreate(&
            arrayspec=realarrspec,grid=LIS_obsvecGrid(n,k),&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status, 'Error ESMF_FieldCreate: ATL15_GrISobs_setup')

!Perturbations State
       write(LIS_logunit,*) '[INFO] Opening attributes for observations ',&
            trim(LIS_rc%obsattribfile(k))
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
       read(ftn,*)
       read(ftn,*) LIS_rc%nobtypes(k)
       read(ftn,*)
    
       allocate(vname(LIS_rc%nobtypes(k)))
       allocate(varmax(LIS_rc%nobtypes(k)))
       allocate(varmin(LIS_rc%nobtypes(k)))
       
       do i=1,LIS_rc%nobtypes(k)
          read(ftn,fmt='(a40)') vname(i)
          read(ftn,*) varmin(i),varmax(i)
          write(LIS_logunit,*) '[INFO] ', vname(i),varmin(i),varmax(i)
       enddo
       call LIS_releaseUnitNumber(ftn)  
       
       allocate(ssdev(LIS_rc%obs_ngrid(k)))

       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 

          allocate(obs_pert%vname(1))
          allocate(obs_pert%perttype(1))
          allocate(obs_pert%ssdev(1))
          allocate(obs_pert%stdmax(1))
          allocate(obs_pert%zeromean(1))
          allocate(obs_pert%tcorr(1))
          allocate(obs_pert%xcorr(1))
          allocate(obs_pert%ycorr(1))
          allocate(obs_pert%ccorr(1,1))

          call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
               obs_pert)

          ssdev = obs_pert%ssdev(1)

          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsEnsOnGrid(n,k),&
               name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status, 'Error ESMF_FieldCreate: ATL15_GrISobs_setup')
          
! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status, 'Error ESMF_FieldGet: ATL15_GrISobs_setup')
          obs_temp(:,:) = 0 

          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1), rc=status)
          call LIS_verify(status,'Error ESMF_AttributeSet: ATL15_GrISobs_setup')
          
          if(LIS_rc%obs_ngrid(k).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%obs_ngrid(k), rc=status)
             call LIS_verify(status,'Error ESMF_AttributeSet: ATL15_GrISobs_setup')
          endif

          call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
               obs_pert%stdmax(1), rc=status)
          call LIS_verify(status,'Error ESMF_AttributeSet: ATL15_GrISobs_setup')
          
          call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
               obs_pert%zeromean(1),rc=status)
          call LIS_verify(status,'Error ESMF_AttributeSet: ATL15_GrISobs_setup')
          
          call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
               obs_pert%tcorr(1),rc=status)
          call LIS_verify(status,'Error ESMF_AttributeSet: ATL15_GrISobs_setup')
          
          call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
               obs_pert%xcorr(1),rc=status)
          
          call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
               obs_pert%ycorr(1),rc=status)

          call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
               obs_pert%ccorr(1,:),itemCount=1,rc=status)

       endif
       
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)

    enddo
    write(LIS_logunit,*) '[INFO] Created the States to hold the observations data'

!--------------------------------------------------------------------------------
! The data will be read and kept in memory at the beginning of the simulation
!--------------------------------------------------------------------------------    
! ATL15 time "days since 2018-01-01". 
! First observation for dhdt is @ 2018-11-16 14:15 ( in the middle of the 1st and 2st dh obervation times)   


     !call ESMF_TimeSet(startTime, yy=2018,&
     !     mm = 11, dd=16, h =0, &
     !     m = 0, s = 0, calendar = LIS_calendar,&
     !     rc=status)
     !call LIS_verify(status, 'ESMF_TimeSet failed in ATL15_GrISobs_setup') 

   do n=1,LIS_rc%nnest

       !ATL15_GrIS_struc(n)%alarmStartdate = startTime
       !!if(LIS_rc%syr.ge.2018.and. LIS_rc%smo.eq.&
       !call LIS_registerAlarm("ATL15 GrIS read alarm",&
       !     7889400.0, 7889400.0)   ! dhdt interval is 91.3125 days or 7889400 s
       ATL15_GrIS_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status,'Error ESMF_StateAdd: ATL15_GrISobs_setup')

       call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
       call LIS_verify(status,'Error ESMF_StateAdd: ATL15_GrISobs_setup')  

    enddo
    
    
  end subroutine ATL15_GrISobs_setup
  
end module ATL15_GrISobs_module
