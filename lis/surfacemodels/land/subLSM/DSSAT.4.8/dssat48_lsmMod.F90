!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: dssat48_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of DSSAT48 model. It also provides the entry method
! for the initialization of DSSAT48-specific variables. The derived
! data type {\tt dssat48\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the dssat48 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for Crocus81
! \item[ts]
!   DSSAT48 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[crocus81]
!  Crocus81 model specific variables
! \item[forc\_count]
!   counter of forcing data


! !REVISION HISTORY:
!  10 May 2020::  Pang-Wei Liu; Created for DSSAT Model
!

module dssat48_lsmMod

! !USES:
  use dssat48_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  USE ModuleDefs, only: ControlType,SwitchType,SoilType !From DSSAT
  implicit none
  
  PRIVATE
  !-------------------------------------------------------------------------
  ! PUBLIC MEMBER FUNCTIONS
  !-------------------------------------------------------------------------
  public :: dssat48_init
  !-------------------------------------------------------------------------
  ! PUBLIC TYPES
  !-------------------------------------------------------------------------
  public :: dssat48_struc
  !EOP
  type, public :: dssat48_type_dec
     !-------------------------------------------------------------------------
     !  LIS related parameters:
     !-------------------------------------------------------------------------
     character(len=LIS_CONST_PATH_LEN) :: rfile
     character*256      :: rformat
     !-------------------------------------------------------------------------
     ! Parameter file names
     !-------------------------------------------------------------------------
     character*128      :: LDT_ncvar_CDL
     character*128      :: LDT_ncvar_MUKEY
     !-------------------------------------------------------------------------
     ! ts, Count, rstInterval, outInterval
     !-------------------------------------------------------------------------
     real               :: ts
     integer            :: count
     real               :: rstInterval
     integer            :: outInterval
     integer            :: forc_count
     !-------------------------------------------------------------------------
     ! Initial Model State for cold start
     !-------------------------------------------------------------------------
     !REAL, pointer      :: init_SNOWSWE(:)
     !REAL               :: init_SNOWALB
     !-------------------------------------------------------------------------
     ! Constant Parameter
     !-------------------------------------------------------------------------
     !integer            :: nsnow
     
     ! for each grid cell in the driver using PERMSNOWFRAC    
     !LOGICAL            :: PRODSNOWMAK_BOOL 
     type(dssat48dec), pointer :: dssat48(:)
     TYPE (ControlType),pointer ::  CONTROL(:) !Use from DSSAT ModuleDefs
     TYPE (SwitchType), pointer ::  ISWITCH(:) !Use from DSSAT ModuleDefs
     TYPE (SoilType), pointer :: SOILPROP(:)   !Use from DSSAT ModuleDefs
  end type dssat48_type_dec
  
  type(dssat48_type_dec), allocatable :: dssat48_struc(:)
contains 
  
!BOP
!
! !ROUTINE: dssat48_init
! \label{dssat48_init}
!
! !INTERFACE:
  subroutine dssat48_init(eks)
! !USES:
    use ESMF
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod
    use LIS_surfaceModelDataMod
    use LIS_lsmMod
    USE ModuleDefs, only: ControlType,SwitchType,SoilType, ModelVerTxt !From DSSAT
    USE CSMVersion
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for dssat48-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for dssat48 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[dssat48\_readcrd](\ref{dssat48_readcrd})\\
!    reads the runtime options for dssat48 model
!  \end{description}
!EOP
    implicit none        
    integer, intent(in) :: eks    
    integer  :: n, t     
    character*3             :: fnest
    integer  :: status   
    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: smField
    ! DEFINITION FOR DSSAT INIT
    CHARACTER*1   :: RNMODE
    CHARACTER*8   :: MODELARG
    CHARACTER*12  :: FILEX
    CHARACTER*30  :: FILEIO
    CHARACTER*80  :: PATHEX
    CHARACTER*120 :: FILECTL
    INTEGER       :: ROTNUM, TRTNUM, YRSIM, YRDOY, MULTI, YRDIF
    INTEGER       :: ERRCODE, RUNINIT, YREND, EXPNO, TRTALL, NYRS, ENDYRS
    INTEGER       :: RUN, YRDOY_END, NREPS, REPNO, TRTREP, YRPLT, MDATE
    INTEGER       :: ERRNUM
    LOGICAL       :: FEXIST
    !CHARACTER(LEN=3)  ModelVerTxt

    !The variable "CONTROL" is of type "ControlType".
    !TYPE (ControlType) CONTROL

    !The variable "ISWITCH" is of type "SwitchType".
    !TYPE (SwitchType) ISWITCH

    ! allocate memory for nest 
    allocate(dssat48_struc(LIS_rc%nnest))
    
    ! read configuation information from lis.config file
    call dssat48_readcrd()
    PRINT*,'Im in dssat48_init'
        do n=1, LIS_rc%nnest
            !Get  Dssat version # which is used to read dssat prms files
            WRITE(ModelVerTxt,'(I2.2,I1)') Version%Major, Version%Minor !Obtain Dssat version #

            ! allocate memory for all tiles in current nest 
            allocate(dssat48_struc(n)%dssat48(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%CONTROL(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%ISWITCH(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%SOILPROP(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !PRINT*, 'LIS_rc%lsm_index is the number of tiles'
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            
            ! allocate memory for state variables
            !do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                !allocate(dssat48_struc(n)%dssat48(t)%SNOWSWE(dssat48_struc(n)%XXXX))
                
            !enddo
            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                dssat48_struc(n)%dssat48(t)%tair = 0.0
                dssat48_struc(n)%dssat48(t)%tmax = 0.0
                dssat48_struc(n)%dssat48(t)%tmin= 0.0          
                dssat48_struc(n)%dssat48(t)%qair = 0.0
                dssat48_struc(n)%dssat48(t)%uwind = 0.0
                dssat48_struc(n)%dssat48(t)%vwind = 0.0
                dssat48_struc(n)%dssat48(t)%rainf = 0.0
                dssat48_struc(n)%dssat48(t)%snowf = 0.0
                dssat48_struc(n)%dssat48(t)%lwdown = 0.0
                dssat48_struc(n)%dssat48(t)%swdown = 0.0
                dssat48_struc(n)%dssat48(t)%psurf = 0.0
            enddo ! end of tile (t) loop

            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                 PRINT*, 'LIS_rc%npatch(n, LIS_rc%lsm_index): ', LIS_rc%npatch(n, LIS_rc%lsm_index)
                 !These System Initial Values May Be Loaded From Config File
                 ! Start CSM 
                 !DONE =  .FALSE.  ! We don't use DONE to control run or not run
                 !YRDOY_END = 9999999 !Can be removed
                 RNMODE = 'Q'
                 FILEIO = 'DSSAT48.INP'
                 ROTNUM = 1
                 TRTNUM = 1
                 FILEX = 'NASA2019.SQX'
               !----------------------------------------------------------------------------------------
               ! Check If There is FILEIO and Delete It  
                 INQUIRE (FILE = FILEIO,EXIST = FEXIST)
                 IF (FEXIST) THEN
                     OPEN (21, FILE = FILEIO,STATUS = 'UNKNOWN',IOSTAT=ERRNUM)
                     CLOSE (21,STATUS = 'DELETE')
                 ENDIF
               !----------------------------------------------------------------------------------------
               ! Initialization CSM
                 YREND = -99
                 YRPLT = 0 !Initialization, Pang
                 MDATE = 0 !Initialization, Pang
                 RUN = 1
                 REPNO = 1
                 !PATHEX = ""
                 dssat48_struc(n)%CONTROL(t)%repno = REPNO
                 dssat48_struc(n)%CONTROL(t)%run = RUN
                 dssat48_struc(n)%CONTROL(t)%yrdoy = 0 !YRDOY is initialized
                 dssat48_struc(n)%CONTROL(t)%fileio = FILEIO
                 dssat48_struc(n)%CONTROL(t)%filex  = FILEX
                 dssat48_struc(n)%CONTROL(t)%rnmode = RNMODE
                 dssat48_struc(n)%CONTROL(t)%rotnum = ROTNUM
                 dssat48_struc(n)%CONTROL(t)%trtnum = TRTNUM
                 dssat48_struc(n)%CONTROL(t)%errcode = 0
             !------------------------------------------------------------------------------------------
             !------- This Section Is Needed When We Still Need To Use .SQX and .INP Files -------------
                !Input Module Reads Experimental File (.SQX) and Write to Temporary IO File (.INP) 
                CALL INPUT_SUB( &
                        FILECTL, FILEIO, FILEX, MODELARG, PATHEX,       &         !Input
                        RNMODE , ROTNUM, RUN, TRTNUM,                   &         !Input
                        dssat48_struc(n)%ISWITCH(t), dssat48_struc(n)%CONTROL(t)) !Output
                !Check to see if the temporary file exists
                !Needed when we still use .INP file
                 INQUIRE (FILE = FILEIO,EXIST = FEXIST)
                 IF (.NOT. FEXIST) THEN
                     !CALL ERROR(ERRKEY,2,FILEIO,LUNIO)
                     PRINT*, FILEIO, ' Does Not Exist!'
                     CALL EXIT
                 ENDIF
              !---------------------------------------------------------------------------------------------
                 EXPNO = 1
                 TRTALL = 999
                 NYRS = 1
                 NREPS = 1
                 YRSIM = 2019110 !Day of Simulation
                 YRDOY_END = (INT(YRSIM/1000)+NYRS-1)*1000 + YRSIM-INT(YRSIM/1000.0)*1000 - 1 
                 YRDOY = YRSIM !YRDOY is initialized as same as YRSIM
                 MULTI = 0
                 YRDIF = 0
                 ENDYRS = 0

                 dssat48_struc(n)%CONTROL(t)%filex  = FILEX
                 dssat48_struc(n)%CONTROL(t)%multi = MULTI
                 dssat48_struc(n)%CONTROL(t)%run = RUN
                 dssat48_struc(n)%CONTROL(t)%trtnum = TRTNUM
                 dssat48_struc(n)%CONTROL(t)%yrdif = YRDIF
                 dssat48_struc(n)%CONTROL(t)%nyrs = NYRS
                 dssat48_struc(n)%CONTROL(t)%yrdoy = YRDOY !For Q Model; Run=1
                 dssat48_struc(n)%CONTROL(t)%yrsim = YRSIM
                 dssat48_struc(n)%CONTROL(t)%endyrs = ENDYRS
                 dssat48_struc(n)%CONTROL(t)%dynamic = 1 !1: RUNINIT


                 dssat48_struc(n)%dssat48(t)%yrend = YREND
                 dssat48_struc(n)%dssat48(t)%expno = EXPNO
                 dssat48_struc(n)%dssat48(t)%trtall= TRTALL
                 dssat48_struc(n)%dssat48(t)%nreps = NREPS
                 dssat48_struc(n)%dssat48(t)%yrdoy_end = YRDOY_END

             !-------------------- LAND INITIALIZATION ---------------------------------------------------
                CALL LAND(dssat48_struc(n)%CONTROL(t), dssat48_struc(n)%ISWITCH(t), &
                  YRPLT, MDATE, YREND, n, t) !Pang: add n, t for ensembles and tiles
                dssat48_struc(n)%dssat48(t)%yrend = YREND
                dssat48_struc(n)%dssat48(t)%mdate = MDATE
                dssat48_struc(n)%dssat48(t)%yrplt = YRPLT
                dssat48_struc(n)%dssat48(t)%doseasinit = .TRUE.
                !PRINT*, 'SOILPROP af LAND in lsmMod: ', dssat48_struc(n)%SOILPROP(t)
                PRINT*, 'ISWITCH af LAND in lsmMod: ', dssat48_struc(n)%ISWITCH(t)
                PRINT*, 'CONTROL af LAND in lsmMod lsmMod2: ', dssat48_struc(n)%CONTROL(t)
                PRINT*, 'YRPLT, MDATE, YREND', YRPLT, MDATE, YREND
                !PRINT*, 'dssat48_struc(ens)%dssat48(t)%BD_INIT: ', dssat48_struc(n)%dssat48(t)%BD_INIT
                !PRINT*, 'SOILPROP%BD: ', dssat48_struc(n)%SOILPROP(t)%BD.
             !---------------------------------------------------------------------------------------------
            enddo

            !------------------------------------------------------------------------
            ! Model timestep Alarm
            !------------------------------------------------------------------------
            dssat48_struc(n)%forc_count = 0
            PRINT*, 'dssat48_lsm mod inthe loop'
            call LIS_update_timestep(LIS_rc, n, dssat48_struc(n)%ts)

            write(fnest,'(i3.3)') n 
            call LIS_registerAlarm("DSSAT48 model alarm "//trim(fnest),&
                                   dssat48_struc(n)%ts, &
                                   dssat48_struc(n)%ts)
            call LIS_registerAlarm("DSSAT48 restart alarm "//trim(fnest),&
                                   dssat48_struc(n)%ts,&
                                   dssat48_struc(n)%rstInterval)
            PRINT*, 'dssat48_lsm mod after restart alarm'
            LIS_sfmodel_struc(n)%ts = dssat48_struc(n)%ts
!----------------------------------------------------------------------------------
!           Create fields for LSM2SUBLSM exchanges
!---------------------------------------------------------------------------------
           ! call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
           !      rc=status)
           ! call LIS_verify(status, &
           !      "ESMF_ArraySpecSet failed in dssat48_in")

           ! smField = ESMF_FieldCreate(&
           !      grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
           !      arrayspec=arrspec1, &
           !      name="Soil Moisture",&
           !      rc=status)
           ! call LIS_verify(status,&
           !      'ESMF_FieldCreate failed for SM in dssat48_init')
           ! 
           ! call ESMF_StateAdd(LIS_LSM2SUBLSM_State(n,eks),&
           !      (/smField/),rc=status)
           ! call LIS_verify(status,&
           !      'ESMF_StateAdd failed for SM in dssat48_init')

        enddo
    end subroutine dssat48_init
end module dssat48_lsmMod
