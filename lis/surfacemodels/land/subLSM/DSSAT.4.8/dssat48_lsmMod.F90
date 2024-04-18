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
  use LIS_coreMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  USE ModuleDefs, only: ControlType,SwitchType,SoilType, &
                        MulchType, OxLayerType,ResidueType, &
                        TillType, FertType, OrgMatAppType !From DSSAT
  USE FloodModule !From DSSAT
  USE GHG_types_mod, only: CH4_type, N2O_type !From DSSAT
  use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc !JE for exchange 
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
     integer :: sm_coupling
     !-------------------------------------------------------------------------
     !  DSSAT I/O
     !-------------------------------------------------------------------------
     character*128      :: outpath !path for INP, INH, and .OUT files
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
     type(dssat48dec), pointer :: dssat48(:)
     TYPE (ControlType),pointer ::  CONTROL(:) !Use from DSSAT ModuleDefs
     TYPE (SwitchType), pointer ::  ISWITCH(:) !Use from DSSAT ModuleDefs
     TYPE (SoilType), pointer :: SOILPROP(:)   !Use from DSSAT ModuleDefs
     TYPE (MulchType), pointer :: MULCH(:)     !Use from DSSAT ModuleDefs
     TYPE (ResidueType), pointer :: HARVRES(:) !Use from DSSAT ModuleDefs
     TYPE (ResidueType), pointer :: SENESCE(:) !Use from DSSAT ModuleDefs
     TYPE (TillType), pointer :: TILLVALS(:) !Use from DSSAT ModuleDefs
     TYPE (FertType), pointer :: FERTDATA(:) !Use from DSSAT ModuleDefs
     TYPE (OrgMatAppType), pointer :: OMADATA(:) !Use from DSSAT ModuleDefs
     TYPE (FloodWatType), pointer :: FLOODWAT(:) !Use from DSSAT FloodModule
     TYPE (FloodNType), pointer :: FLOODN(:) !Use from DSSAT FloodModule
     TYPE (CH4_type), pointer :: CH4_data(:) !Use from DSSAT GHG_mod/dssat48_module
     TYPE (N2O_type), pointer :: N2O_data(:) !Use from DSSAT GHG_mod/dssat48_module
     TYPE (OxLayerType), pointer :: OXLAYR(:) !Use from DSSAT ModuleDefs
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
    use LIS_fileIOMod, only : LIS_create_output_directory !JG
    USE ModuleDefs, only: ModelVerTxt, MonthTxt !From DSSAT
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
    integer  :: n, t, l     
    character*3             :: fnest
    integer  :: status   
    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: smField
    integer              :: tmp_year, tmp_month, tmp_day, tmp_hour, tmp_minute
    ! DEFINITION FOR DSSAT INIT
    !CHARACTER*1   :: RNMODE
    !CHARACTER*8   :: MODELARG
    !CHARACTER*12  :: FILEX
    CHARACTER*120  :: FILEIO
    CHARACTER*4   :: fproc      !JE
    !CHARACTER*80  :: PATHEX
    !CHARACTER*120 :: FILECTL
    !INTEGER       :: ROTNUM, TRTNUM, YRSIM, YRDOY, MULTI, YRDIF
    !INTEGER       :: ERRCODE, RUNINIT, YREND, EXPNO, TRTALL, NYRS, ENDYRS
    !INTEGER       :: RUN, YRDOY_END, NREPS, REPNO, TRTREP, YRPLT, MDATE
    INTEGER       :: RUN, REPNO, ROTNUM, TRTNUM, ERRNUM
    !INTEGER       :: ERRNUM, JULIAN
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
    !PRINT*,'Im in dssat48_init'
    
        do n=1, LIS_rc%nnest
            !Get  Dssat version # which is used to read dssat prms files
            WRITE(ModelVerTxt,'(I2.2,I1)') Version%Major, Version%Minor !Obtain Dssat version #

            dssat48_struc(n)%outpath = trim(LIS_rc%odir)
            write(LIS_logunit,*) "Writing DSSAT INP, INH, and .OUT files to ", dssat48_struc(n)%outpath
            ! allocate memory for all tiles in current nest 
            allocate(dssat48_struc(n)%dssat48(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%CONTROL(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%ISWITCH(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%SOILPROP(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%MULCH(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%FLOODWAT(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%FLOODN(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%CH4_data(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%N2O_data(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%OXLAYR(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%HARVRES(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%SENESCE(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%TILLVALS(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%FERTDATA(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            allocate(dssat48_struc(n)%OMADATA(LIS_rc%npatch(n, LIS_rc%lsm_index)))
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
                dssat48_struc(n)%dssat48(t)%totprc = 0.0
                dssat48_struc(n)%dssat48(t)%tdew = 0.0
                if (dssat48_struc(n)%sm_coupling.eq.1) then
                   do l=0, LIS_sfmodel_struc(n)%nsm_layers
                      dssat48_struc(n)%dssat48(t)%LIS_sm(l) = 0.0
                   end do
                endif
            enddo ! end of tile (t) loop
                  !PRINT*, 'LIS_rc%npatch(n, LIS_rc%lsm_index): ', LIS_rc%npatch(n, LIS_rc%lsm_index)
          
            ! Call dssat48_setup to obtain mukey number before initialization
            !  CALL dssat48_setup !Pang 2024.02.09 (This is used to obtain soil mukey)

            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            !do t=1,3 !PL for testing code
                 !Start CSM 
                 dssat48_struc(n)%dssat48(t)%DONE =  .FALSE.  !PL: Control code to run DSSAT48 or NOT
                 dssat48_struc(n)%dssat48(t)%YRDOY_END = 9999999 !PL: Initial YRDOY_END
                 dssat48_struc(n)%CONTROL(t)%rnmode = 'Q' !run mode is fixed for LIS-DSSAT

                 RUN = 0
                 REPNO = 1
                 dssat48_struc(n)%CONTROL(t)%repno = REPNO
                 dssat48_struc(n)%CONTROL(t)%run = RUN
 
                 ROTNUM = 0
                 TRTNUM = 0
                 dssat48_struc(n)%CONTROL(t)%rotnum = ROTNUM
                 dssat48_struc(n)%CONTROL(t)%trtnum = TRTNUM
                 dssat48_struc(n)%CONTROL(t)%filex  = 'NASA2019.SQX' !This can be set in config file
                !----------------------------------------------------------------------------------------
                !JE Add one .INP file per processor
                 !print*, 'Writing INP Files to ', trim(dssat48_struc(n)%outpath)
                 !JE Add one .INP file per processor
                 write(unit=fproc,fmt='(i4.4)') LIS_localPet
                 FILEIO = trim(dssat48_struc(n)%outpath)//'/INP/'//'DSSAT48.INP.'//fproc !JG: add '/INP/'
                 !print*, 'FILEIO: ', FILEIO

                 !Check If There is FILEIO and Delete It  
                 INQUIRE (FILE = FILEIO,EXIST = FEXIST)
                 IF (FEXIST) THEN
                     OPEN (21, FILE = FILEIO,STATUS = 'UNKNOWN',IOSTAT=ERRNUM)
                     CLOSE (21,STATUS = 'DELETE')
                 ENDIF
            enddo

            !------------------------------------------------------------------------
            ! Model timestep Alarm
            !------------------------------------------------------------------------
            dssat48_struc(n)%forc_count = 0
            !PRINT*, 'dssat48_lsm mod inthe loop'
            call LIS_update_timestep(LIS_rc, n, dssat48_struc(n)%ts)

            write(fnest,'(i3.3)') n 
            call LIS_registerAlarm("DSSAT48 model alarm "//trim(fnest),&
                                   dssat48_struc(n)%ts, &
                                   dssat48_struc(n)%ts)
            call LIS_registerAlarm("DSSAT48 restart alarm "//trim(fnest),&
                                   dssat48_struc(n)%ts,&
                                   dssat48_struc(n)%rstInterval)
            !PRINT*, 'dssat48_lsm mod after restart alarm'
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
        if ( LIS_masterproc ) then
           call LIS_create_output_directory('INP')
        endif
    end subroutine dssat48_init
end module dssat48_lsmMod
