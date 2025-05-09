!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISdssatSM_obsMod
!BOP
! 
! !MODULE: LISdssatSM_obsMod
! 
! !DESCRIPTION: 
!  This module handles the use of a LIS-DSSAT model simulation 
!  output as "observations" for data assimilation. This
!  plugin is typically used to handle the computations of 
!  scaling factors such as cumulative distribution function
!  (CDF) for use in DA
! 
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
!  03 Oct 2023    Pang-Wei Liu Revise to DSSAT SM

  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISdssatSM_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: dssatsmobs
!
!EOP
  
  type, public :: dssatsmobsdec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr
     integer       :: soillyer !P-W. Liu
     real          :: datares
     real          :: run_dd(8)
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character(len=LDT_CONST_PATH_LEN) :: odir
     character*20  :: security_class
     character*20  :: distribution_class
     character*20  :: data_category
     character*20  :: area_of_data
     character*20  :: write_interval
     character*20  :: lsmordssat !P-W. Liu
!--------------------------------------------------------
!  interpolation/upscaling weights
!--------------------------------------------------------
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real,  allocatable     :: w11(:)
     real,  allocatable     :: w12(:)
     real,  allocatable     :: w21(:)
     real,  allocatable     :: w22(:)

  end type dssatsmobsdec

  type(dssatsmobsdec)  :: dssatsmobs

contains

!BOP
! !ROUTINE: LISdssatSM_obsInit
! \label{LISdssatSM_obsInit}
! 
! !INTERFACE: 
  subroutine LISdssatSM_obsInit()
! !USES: 
    use ESMF
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_logMod

    implicit none
! 
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of a 
! land surface model and dssat outputs (from a LIS-DSSAT simulation) as observations.  
! 
!EOP
    integer                 :: n
    integer                 :: rc
    real                    :: gridDesci(20)

    n = 1

    dssatsmobs%run_dd             = LDT_rc%udef

    dssatsmobs%security_class     = ''
    dssatsmobs%distribution_class = ''
    dssatsmobs%data_category      = ''
    dssatsmobs%area_of_data       = ''
    dssatsmobs%write_interval     = ''

    call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%format, &
         label="LIS-DSSAT soil moisture output format:",rc=rc)
    call LDT_verify(rc,'LIS-DSSAT soil moisture output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%wopt, &
         label="LIS-DSSAT soil moisture output methodology:",rc=rc)
    call LDT_verify(rc,'LIS-DSSAT soil moisture output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%wstyle, &
         label="LIS-DSSAT soil moisture output naming style:",rc=rc)
    call LDT_verify(rc,'LIS-DSSAT soil moisture output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%map_proj, &
         label="LIS-DSSAT soil moisture output map projection:",rc=rc)
    call LDT_verify(rc,'LIS-DSSAT soil moisture output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%nest, &
         label="LIS-DSSAT soil moisture output nest index:",rc=rc)
    call LDT_verify(rc,'LIS-DSSAT soil moisture output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%odir, &
         label="LIS-DSSAT soil moisture output directory:",rc=rc)
    call LDT_verify(rc,'LIS-DSSAT soil moisture output directory: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%lsmordssat, &
         label="LIS-DSSAT soil moisture from LSM/DSSAT:",rc=rc)
    call LDT_verify(rc,'LIS-DSSAT soil moisture from LSM/DSSAT: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%soillyer, &
         label="LIS-DSSAT soil moisture layer (1-4):",rc=rc)
    call LDT_verify(rc,'LIS-DSSAT soil moisture layer: not defined')

    ! WMO-convention specific identifiers
    if ( dssatsmobs%wstyle == "WMO convention") then 
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%security_class, &
       label="LIS-DSSAT soil moisture security class:",rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture security class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%distribution_class, &
       label="LIS-DSSAT soil moisture distribution class:",rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture distribution class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%data_category, &
       label="LIS-DSSAT soil moisture data category:",rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture data category: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%area_of_data, &
       label="LIS-DSSAT soil moisture area of data:",rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture area of data: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%write_interval, &
       label="LIS-DSSAT soil moisture write interval:",rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture write interval: not defined')
    endif

    if(dssatsmobs%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(6),rc=rc)       
       
       dssatsmobs%datares = min(dssatsmobs%run_dd(5),dssatsmobs%run_dd(6))

       dssatsmobs%nc    = (nint((dssatsmobs%run_dd(4)-dssatsmobs%run_dd(2))/&
            dssatsmobs%run_dd(5))) + 1
       dssatsmobs%nr    = (nint((dssatsmobs%run_dd(3)-dssatsmobs%run_dd(1))/&
            dssatsmobs%run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = dssatsmobs%nc
       gridDesci(3) = dssatsmobs%nr
       gridDesci(4) = dssatsmobs%run_dd(1)
       gridDesci(5) = dssatsmobs%run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = dssatsmobs%run_dd(3)
       gridDesci(8) = dssatsmobs%run_dd(4)
       gridDesci(9) = dssatsmobs%run_dd(5)
       gridDesci(10) = dssatsmobs%run_dd(6)
       gridDesci(20) = 64

    elseif(dssatsmobs%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS-DSSAT soil moisture domain y-dimension size: not defined')

       dssatsmobs%datares = dssatsmobs%run_dd(6)/100.0

       dssatsmobs%nc    = dssatsmobs%run_dd(7)
       dssatsmobs%nr    = dssatsmobs%run_dd(8)

       gridDesci = 0
       gridDesci(1) = 3 
       gridDesci(2) = dssatsmobs%nc
       gridDesci(3) = dssatsmobs%nr
       gridDesci(4) = dssatsmobs%run_dd(1)
       gridDesci(5) = dssatsmobs%run_dd(2)
       gridDesci(6) = 8
       gridDesci(7) = dssatsmobs%run_dd(4)
       gridDesci(8) = dssatsmobs%run_dd(6)
       gridDesci(9) = dssatsmobs%run_dd(6)
       gridDesci(10) = dssatsmobs%run_dd(3)
       gridDesci(11) = dssatsmobs%run_dd(5)
       gridDesci(20) = 64

    elseif(dssatsmobs%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(1),rc=rc)
 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS-DSSAT soil moisture domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,dssatsmobs%run_dd(8),rc=rc)

    endif

!-------------------------------------------------------------------
!  if the LIS-DSSAT output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,dssatsmobs%datares)) then 

       allocate(dssatsmobs%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(dssatsmobs%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(dssatsmobs%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(dssatsmobs%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(dssatsmobs%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(dssatsmobs%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(dssatsmobs%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(dssatsmobs%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            dssatsmobs%n11, &
            dssatsmobs%n12, dssatsmobs%n21, &
            dssatsmobs%n22, dssatsmobs%w11, &
            dssatsmobs%w12, dssatsmobs%w21, &
            dssatsmobs%w22)

    else

       allocate(dssatsmobs%n11(&
            dssatsmobs%nc*&
            dssatsmobs%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            dssatsmobs%nc*dssatsmobs%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            dssatsmobs%n11)
       
    endif
    
!  which variable we want in the DA obs computations. 
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%soilmoist_obs, "m3/m3",1,1)
    LDT_DAobsData(1)%soilmoist_obs%selectStats = 1

  end subroutine LISdssatSM_obsInit
  
end module LISdssatSM_obsMod
