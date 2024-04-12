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
! !ROUTINE: Crocus81_main
! \label{Crocus81_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   10/18/19: Mahdi Navari, Shugong Wang; initial implementation for Crocus81 with LIS-7
!   9 Dec 2020: Mahdi Navari; edited to take into account the Crocus slope correction
!   19 Jan 2021: Mahdi Navari, edited to properly initialize precipitation 
!   21 Jan 2021: Mahdi Navari, edited to properly assign values for TG, XWG, and XWGI for the stand-alone version
!
! !INTERFACE:
subroutine dssat48_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod 
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
    use dssat48_lsmMod
    use NoahMP401_lsmMod ! will need to generalize JE
    USE ModuleDefs, only: MonthTxt !From DSSAT
   !use other modules
  
    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    integer              :: l
    real                 :: dt
    real                 :: lat, lon, elev
    integer              :: row, col
    integer              :: year, month, day, hour, minute, second
    logical              :: alarmCheck
    integer              :: nyrs_lis, yrdoy_end_lis 

    real                 :: tmp_LAT                ! Latitude in decimal degree  (latitude (degrees +North)) [degrees]
    real                 :: tmp_LON                ! Longitude in decimal year    (longitude (degrees +East)) [degrees]
    integer              :: tmp_year, tmp_month, tmp_day, tmp_hour, tmp_minute
    integer              :: year_end, month_end, day_end
    real                 :: tmp_pres, tmp_precip, tmp_tmax, tmp_tmin, tmp_tdew   ! Weather Forcing
    real                 :: tmp_swrad, tmp_wind, tmp_sm                          ! Weather Forcing
    character*3          :: fnest ! MN Bug in toolkit (added to this code)
! CONTRO
    INTEGER             :: RUN, REPNO, EXPNO, TRTALL, NREPS
    CHARACTER*1         :: RNMODE
    INTEGER             :: YREND, MDATE, YRPLT, YRDOY, JULIAN
    INTEGER             :: YRDOY_END, YRDIF, YRSIM_SAVE
    INTEGER             :: NYRS, ENDYRS, MULTI, YRSIM, YR, ISIM
    INTEGER             :: DAS, DOY, TIMDIF, INCYD, YR0, ISIM0
    LOGICAL             :: doseasinit, FEXIST, DONE

    CHARACTER*120 :: FILECTL
    CHARACTER*30  :: FILEIO
    CHARACTER*4   :: fproc      !JE
    CHARACTER*12  :: FILEX
    CHARACTER*8   :: MODELARG
    CHARACTER*80  :: PATHEX
    INTEGER       :: ROTNUM, TRTNUM
! !DESCRIPTION:
!  This is the entry point for calling the dssat48 physics.
!  This routine calls XXXXXX routine that performs the
!  land dssat.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for dssat48
    
    !allocate( tmp_SNOWSWE( dssat48_struc(n)%nsnow ) )

!The variable "CONTROL" is of type "ControlType".
    !TYPE (ControlType) CONTROL

    !The variable "ISWITCH" is of type "SwitchType".
    !TYPE (SwitchType) ISWITCH

    ! check dssat48 alarm. If alarm is ring, run model.
    write(LIS_logunit,*) '[INFO] Call to the DSSAT48 Check Main routine ...'  
    write(fnest,'(i3.3)') n

    alarmCheck = LIS_isAlarmRinging(LIS_rc, "DSSAT48 model alarm "// trim(fnest)) !MN  Bug in the toolkit 
    if (alarmCheck) Then

        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        !do t=1, 1 !PL for testing code     
            dt = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
            elev = LIS_domain(n)%tile(t)%elev

            ! Spatial Information
            !PRINT*, "lat, lon, elev: ", lat, lon, elev
            dssat48_struc(n)%dssat48(t)%lat = lat
            dssat48_struc(n)%dssat48(t)%lon = lon
            dssat48_struc(n)%dssat48(t)%elev = elev

            !!------ This Block Is Where We Obtain Weather Forcing ------------------------------!!
            ! retrieve forcing data from DSSAT48_struc(n)%dssat(t) and assign to local variables

            write(LIS_logunit,*) 'Weather forcing for tile: ',t

            ! Daily average surface pressure (Pa)
            tmp_pres      = dssat48_struc(n)%dssat48(t)%psurf / dssat48_struc(n)%forc_count
            write(LIS_logunit,*) 'p: ',tmp_pres
            dssat48_struc(n)%dssat48(t)%forc_pres = tmp_pres
 
            ! Total daily precipitation (rain+snow) (mm)
            tmp_precip    = (dssat48_struc(n)%dssat48(t)%totprc / dssat48_struc(n)%forc_count) * 3600. * 24. !Convert from kg/ms2 to mm
            write(LIS_logunit,*) 'Precip: ',tmp_precip
            dssat48_struc(n)%dssat48(t)%forc_precip = tmp_precip

            ! Tmax: maximum daily air temperature (C)
            tmp_tmax      = dssat48_struc(n)%dssat48(t)%tmax - 273.15 !Convert from K to C
            write(LIS_logunit,*) 'Tmax: ',tmp_tmax
            dssat48_struc(n)%dssat48(t)%forc_tmax = tmp_tmax

            ! Tmin: minimum daily air temperature (C)
            tmp_tmin      = dssat48_struc(n)%dssat48(t)%tmin - 273.15 !Convert from K to C 
            write(LIS_logunit,*) 'Tmin: ',tmp_tmin
            dssat48_struc(n)%dssat48(t)%forc_tmin = tmp_tmin

            ! Tdew: average daily dewpoint temperature (C)
            tmp_tdew      = (dssat48_struc(n)%dssat48(t)%tdew / dssat48_struc(n)%forc_count) - 273.15 !Convert from K to C
            write(LIS_logunit,*) 'Tdew: ',tmp_tdew
            dssat48_struc(n)%dssat48(t)%forc_tdew = tmp_tdew

            ! SW_RAD: daily total incoming solar radiation (MJ/(m2d))
            tmp_swrad     = (dssat48_struc(n)%dssat48(t)%swdown / dssat48_struc(n)%forc_count) * 0.0864 !Convert from W/m2 to MJ/(m2d)
            write(LIS_logunit,*) 'Swrad: ',tmp_swrad
            dssat48_struc(n)%dssat48(t)%forc_swrad = tmp_swrad

            ! Wind: daily average wind speed (km/d)
            tmp_wind      = (dssat48_struc(n)%dssat48(t)%wndspd / dssat48_struc(n)%forc_count) * 86.0 !Convert from m/s to km/d
            write(LIS_logunit,*) 'Wind: ',tmp_wind
            dssat48_struc(n)%dssat48(t)%forc_wind = tmp_wind

            !!---- Here Will Check Validity of Forcings ------------------!!
            ! check validity of PPS
            !if(tmp_PPS .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable PPS in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            ! check validity of SRSNOW
            !if(tmp_SRSNOW .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable SRSNOW in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            ! check validity of RRSNOW
            !if(tmp_RRSNOW .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable RRSNOW in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            !! check validity of TA
            !if(tmp_TA .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable TA in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            !! check validity of SW_RAD
            !if(tmp_SW_RAD .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable SW_RAD in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            !! check validity of QA
            !if(tmp_QA .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable QA in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            !! check validity of Wind_E
            !if(tmp_Wind_E .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable Wind_E in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            ! check validity of Wind_N
            !if(tmp_Wind_N .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable Wind_N in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            ! check validity of LW_RAD
            !if(tmp_LW_RAD .eq. LIS_rc%udef) then
            !    write(LIS_logunit, *) "undefined value found for forcing variable LW_RAD in Crocus81"
            !    write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
            !    call LIS_endrun()
            !endif
            ! 
            tmp_LAT  = lat
            tmp_LON = lon
            tmp_year   = LIS_rc%yr
            tmp_month  = LIS_rc%mo
            tmp_day    = LIS_rc%da
            tmp_hour   = LIS_rc%hr
            tmp_minute = LIS_rc%mn
            !PRINT*, 'time: ', tmp_year, tmp_month, tmp_day, tmp_hour, tmp_minute

            ! get parameters 
            !tmp_nsnow                               = CROCUS81_struc(n)%nsnow                           


            !Here Determines if Coupling LSM              
            !if(LIS_rc%lsm.ne."none") then
            !   tmp_TG                                  = CROCUS81_struc(n)%crocus81(t)%TG
            !else
            ! For the stand-alone version, we need to provide some default values for XWGI and XWG
              !tmp_XWGI = 0.0 ! MN set to zero 
              !tmp_XWG  = tmp_POROSITY * 0.8 ! MN assume volumetric soil water content of the snow covered
            !endif
            ! get state variables
            !tmp_SNOWSWE(:)    = CROCUS81_struc(n)%crocus81(t)%SNOWSWE(:)   

            ! JE Add options for tight coupling
            ! Add model flag

            IF (dssat48_struc(n)%sm_coupling.eq.1) THEN
               do l=0, LIS_sfmodel_struc(n)%nsm_layers
                  tmp_sm = (dssat48_struc(n)%dssat48(t)%LIS_sm(l) / dssat48_struc(n)%forc_count)
                  write(LIS_logunit,*) 'Soil Moisture Layer: ',l,tmp_sm
                  dssat48_struc(n)%dssat48(t)%LIS_sm(l) = tmp_sm
                  !dssat48_struc(n)%dssat48(t)%LIS_sm(l) = NOAHMP401_struc(n)%noahmp401(t)%smc(l)
               end do
            ENDIF

!-------------------------------------------------------------------------------------
! call model physics
            YREND = dssat48_struc(n)%dssat48(t)%yrend
            MDATE = dssat48_struc(n)%dssat48(t)%mdate
            YRPLT = dssat48_struc(n)%dssat48(t)%yrplt

            RNMODE = dssat48_struc(n)%CONTROL(t)%rnmode
            YRDOY= tmp_year*1000 + JULIAN (tmp_day,MonthTxt(tmp_month),tmp_year)
            DONE = dssat48_struc(n)%dssat48(t)%DONE
            RUN = dssat48_struc(n)%CONTROL(t)%run
            FILEX = dssat48_struc(n)%CONTROL(t)%filex
            ROTNUM = dssat48_struc(n)%CONTROL(t)%rotnum 
            TRTNUM = dssat48_struc(n)%CONTROL(t)%trtnum
            YRDOY_END = dssat48_struc(n)%dssat48(t)%yrdoy_end
            YRSIM = dssat48_struc(n)%CONTROL(t)%yrsim
            NYRS = dssat48_struc(n)%CONTROL(t)%nyrs
            ENDYRS = dssat48_struc(n)%CONTROL(t)%endyrs
            MULTI = dssat48_struc(n)%CONTROL(t)%multi
            !------ CSM INITIALIZATION -----------------------------------------------
            IF (.NOT.DONE) THEN
                PRINT*, 'Im in INIT'
                !PRINT*, 'RUN: ', RUN
                year_end = LIS_rc%eyr
                month_end = LIS_rc%emo
                day_end = LIS_rc%eda

                YREND = -99
                RUN = RUN + 1
                dssat48_struc(n)%CONTROL(t) % RUN = RUN
                dssat48_struc(n)%CONTROL(t) % YRDOY = 0

                !JE Add one .INP file per processor
                 write(unit=fproc,fmt='(i4.4)') LIS_localPet
                 FILEIO = 'DSSAT48.INP.'//fproc
                 !print*, 'FILEIO, RUN: ', FILEIO, RUN
                 TRTNUM = 1 !Initialization
                 ROTNUM = 1 !Initialization
                 dssat48_struc(n)%CONTROL(t)%fileio = FILEIO
                 !dssat48_struc(n)%CONTROL(t)%filex  = FILEX
                 !dssat48_struc(n)%CONTROL(t)%rnmode = RNMODE
                 dssat48_struc(n)%CONTROL(t)%rotnum = ROTNUM
                 dssat48_struc(n)%CONTROL(t)%trtnum = TRTNUM
                 dssat48_struc(n)%CONTROL(t)%errcode = 0
                !------------------------------------------------------------------------------------------
                !------- This Section Is Needed When We Still Need To Use .SQX and .INP Files -------------
                !Input Module Reads Experimental File (.SQX) and Write to Temporary IO File (.INP) 
                 CALL INPUT_SUB( n, t,                               & !Pang 2024.02.08
                        FILECTL, FILEIO, FILEX, MODELARG, PATHEX,       &         !Input
                        RNMODE , ROTNUM, RUN, TRTNUM,                   &         !Input
                        dssat48_struc(n)%ISWITCH(t), dssat48_struc(n)%CONTROL(t)) !Output
                        !PRINT*, 'What is the mukey now: ', dssat48_struc(n)%dssat48(t)%SLNO
                !Check to see if the temporary file exists
                !Needed when we still use .INP file
                 INQUIRE (FILE = FILEIO,EXIST = FEXIST)
                 IF (.NOT. FEXIST) THEN
                     !CALL ERROR(ERRKEY,2,FILEIO,LUNIO)
                     PRINT*, FILEIO, ' Does Not Exist!'
                     CALL EXIT
                 ENDIF
                !--------------------------------------------------------------------------------------------
                  EXPNO = 1 !Always 1; doesn't really need this
                  TRTALL = 999 !Always 9999
                  NYRS = 1 !Always 1 for Q mode
                  NREPS = 1 !Always 1
                  YRSIM = tmp_year*1000 + JULIAN (tmp_day,MonthTxt(tmp_month),tmp_year) +1
                  !PL: YRSIM +1 Due to 1 day shift between LIS and DSSAT
                  !PRINT*, 'YRSIM: ', YRSIM
                  !YRDOY_END = (INT(YRSIM/1000)+NYRS-1)*1000 + YRSIM-INT(YRSIM/1000.0)*1000 - 1 !Can we obtain this from config. 
                  YRDOY_END = year_end*1000 + JULIAN (day_end,MonthTxt(month_end),year_end)
                  !IF (RUN.EQ.1) THEN
                  !    YRDOY = YRSIM !YRDOY is initialized as same as YRSIM
                  !ENDIF !PL: We don't need to do this

                  MULTI = 0
                  YRDIF = 0
                  ENDYRS = 0 
                  !PRINT*, 'YRSIM, YRDOY_END, YRDOY: ', YRSIM, YRDOY_END, YRDOY
                  IF (RUN .GT. 1) THEN
                      YRSIM = INCYD(YRDOY,0)
                      CALL YR_DOY(YRSIM_SAVE, YR0, ISIM0)
                      CALL YR_DOY(YRSIM,      YR,  ISIM)
                      YRDIF = YR - YR0
                      dssat48_struc(n)%CONTROL(t)%YRDIF = YRDIF
                  ENDIF
                  !PRINT*, 'YRDIF, YRSIM: ', YRDIF, YRSIM
                  !dssat48_struc(n)%CONTROL(t)%filex  = FILEX
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
                 dssat48_struc(n)%dssat48(t)%DONE = .TRUE.
            ENDIF


            !----- SEASONAL INITIALIZATION -------------------------------------------
            IF (dssat48_struc(n)%dssat48(t)%doseasinit) THEN
                !PRINT*, 'In SEAS INIT'
                  !Input Module Reads Experimental File (.SQX) and Write to Temporary IO File (.INP) 
                !PRINT*, 'SEAS YRDOY: ', YRDOY, dssat48_struc(n)%CONTROL(t)%YRDOY
                 !JE Add one .INP file per processor
                 !write(unit=fproc,fmt='(i4.4)') LIS_localPet
                 !FILEIO = 'DSSAT48.INP.'//fproc

                 !FILEX = dssat48_struc(n)%CONTROL(t)%filex !Obtain in front of the loop
                 !ROTNUM = dssat48_struc(n)%CONTROL(t)%rotnum !Obtain in front of the loop
                 !TRTNUM = dssat48_struc(n)%CONTROL(t)%trtnum !Obtain in from of the loop
   
                 !NYRS = dssat48_struc(n)%CONTROL(t)%nyrs !Obtain in front of the loop
                 !ENDYRS = dssat48_struc(n)%CONTROL(t)%endyrs !Obtain in front of the loop
                 !MULTI = dssat48_struc(n)%CONTROL(t)%multi !Obtain in front of the loop
                 !RUN = dssat48_struc(n)%CONTROL(t)%run !Obtain in front of the loop
                 !YRSIM = dssat48_struc(n)%CONTROL(t)%yrsim !Obtain in front of the loop
                 REPNO = dssat48_struc(n)%CONTROL(t)%repno !May not really needed for Q mode
                 !dssat48_struc(n)%CONTROL(t) % YRDOY   = YRDOY !Put back later
                 
                 !CALL INPUT_SUB( n, t,                                  & !Pang 20240207
                 !       FILECTL, FILEIO, FILEX, MODELARG, PATHEX,       &         !Input
                 !       RNMODE , ROTNUM, RUN, TRTNUM,                   &         !Input
                 !       dssat48_struc(n)%ISWITCH(t), dssat48_struc(n)%CONTROL(t)) !Output

                 !!CONDITIONS
                 !!PRINT*, 'NYRS, ENDYRS, MULTI bf: ', NYRS, ENDYRS, MULTI
                 IF (NYRS .GT. 1) THEN
                     ENDYRS = ENDYRS + 1
                     IF (RNMODE .NE. 'Y') THEN
                       MULTI = MULTI + 1
                     ENDIF
                 ELSE
                     MULTI = 1
                     ENDYRS = 1
                 ENDIF
                 !Note: NYRS, ENDYRS, and ENDYRS are always 1 for Q mode
                 !IF (RUN.GT.1) THEN
                 !   YRDOY = dssat48_struc(n)%CONTROL(t)%yrsim !We have done yrsim = yrdoy for RUN>1
                 !ENDIF
                 dssat48_struc(n)%CONTROL(t) % DAS     = 0
                 dssat48_struc(n)%CONTROL(t) % RUN     = RUN
                 dssat48_struc(n)%CONTROL(t) % YRSIM   = YRSIM  !Starting Day of simulation(from config file)
                 dssat48_struc(n)%CONTROL(t) % YRDOY   = YRDOY  !The day of simulation
                 dssat48_struc(n)%CONTROL(t) % MULTI   = MULTI
                 dssat48_struc(n)%CONTROL(t) % DYNAMIC = 2   !SEASINIT
                 dssat48_struc(n)%CONTROL(t) % ENDYRS  = ENDYRS
                 dssat48_struc(n)%CONTROL(t) % REPNO   = REPNO
                 DAS = dssat48_struc(n)%CONTROL(t) % DAS
                 CALL LAND(dssat48_struc(n)%CONTROL(t), dssat48_struc(n)%ISWITCH(t), &
                    YRPLT, MDATE, YREND, n, t) !Pang: add n, t for ensembles and tiles

                 !JE Write variables to LIS_HIST file
                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD1, &
                    value=dssat48_struc(n)%dssat48(t)%SW(1),&
                    vlevel=1,unit="m^3 m-3",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD2, &
                    value=dssat48_struc(n)%dssat48(t)%SW(2),&
                    vlevel=1,unit="m^3 m-3",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD3, &
                    value=dssat48_struc(n)%dssat48(t)%SW(3),&
                    vlevel=1,unit="m^3 m-3",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD4, &
                    value=dssat48_struc(n)%dssat48(t)%SW(4),&
                    vlevel=1,unit="m^3 m-3",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_DSSAT_LAI, &
                    value=dssat48_struc(n)%dssat48(t)%XLAI,&
                    vlevel=1,unit="m2/m2",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GWAD, &
                    value=dssat48_struc(n)%dssat48(t)%SDWT*10.0,&
                    vlevel=1,unit="kg/ha",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                dssat48_struc(n)%dssat48(t)%doseasinit = .FALSE. !Pnng 2023.09.19
             ENDIF   
                !PRINT*, 'YREND in seas: ', YREND, YRDOY
                !PRINT*, 'YRDOY, YRPLT, MDATE, YREND: ', YRDOY, YRPLT, MDATE, YREND
            !----- DAILY  ------------------------------------------------------------
            IF ((.NOT.dssat48_struc(n)%dssat48(t)%doseasinit).AND.(YRDOY.GE.YRSIM)) THEN
               !PRINT*, 'Im in seas daily rate'
               !-----------------------------------------------------------------------
               !     Calculate days after simulation (DAS) 
               !-----------------------------------------------------------------------
                !PRINT*, 'DAS1: ', DAS, YRDOY, DOY
                !YRSIM = dssat48_struc(n)%CONTROL(t)%yrsim !Obtain infront of the loop
                CALL YR_DOY(YRDOY,YEAR,DOY)
                DAS   = MAX(0,TIMDIF(INCYD(YRSIM,-1),YRDOY))
                dssat48_struc(n)%CONTROL(t) % YRDOY = YRDOY
                dssat48_struc(n)%CONTROL(t) % das = DAS
               !PRINT*, 'RATES YRDOY DAS: ', YRDOY, dssat48_struc(n)%CONTROL(t) % YRDOY, DAS
               !-----------------------------------------------------------------------
               !*********************************************************************** 
               !     RATE CALCULATIONS
               !*********************************************************************** 
               dssat48_struc(n)%CONTROL(t)%dynamic = 3 !3: DAILY RATE
               CALL LAND(dssat48_struc(n)%CONTROL(t), dssat48_struc(n)%ISWITCH(t), &
                   YRPLT, MDATE, YREND, n, t) !Pang: add n, t for ensembles and tiles
                !PRINT*, 'YRDOY, YRPLT, MDATE, YREND: ', YRDOY, YRPLT, MDATE, YREND
               !*********************************************************************** 
               !     INTEGRATION 
               !*********************************************************************** 
               dssat48_struc(n)%CONTROL(t)%dynamic = 4 !4: DAILY INTEGR
               CALL LAND(dssat48_struc(n)%CONTROL(t), dssat48_struc(n)%ISWITCH(t), &
                   YRPLT, MDATE, YREND, n, t) !Pang: add n, t for ensembles and tiles
                !PRINT*, 'YRDOY, YRPLT, MDATE, YREND: ', YRDOY, YRPLT, MDATE, YREND

               !*********************************************************************** 
               !     OUTPUT
               !*********************************************************************** 
               dssat48_struc(n)%CONTROL(t)%dynamic = 5 !5: DAILY OUTPUT IF NEEDED
               CALL LAND(dssat48_struc(n)%CONTROL(t), dssat48_struc(n)%ISWITCH(t), &
                   YRPLT, MDATE, YREND, n, t) !Pang: add n, t for ensembles and tiles
               !PRINT*, 'YRDOY, YRPLT, MDATE, YREND: ', YRDOY, YRPLT, MDATE, YREND

               !JE Write variables to LIS_HIST file
               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD1, &
                  value=dssat48_struc(n)%dssat48(t)%SW(1),&
                  vlevel=1,unit="m^3 m-3",direction="-",&
                  surface_type=LIS_rc%lsm_index)

               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD2, &
                  value=dssat48_struc(n)%dssat48(t)%SW(2),&
                  vlevel=1,unit="m^3 m-3",direction="-",&
                  surface_type=LIS_rc%lsm_index)

               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD3, &
                  value=dssat48_struc(n)%dssat48(t)%SW(3),&
                  vlevel=1,unit="m^3 m-3",direction="-",&
                  surface_type=LIS_rc%lsm_index)

               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD4, &
                  value=dssat48_struc(n)%dssat48(t)%SW(4),&
                  vlevel=1,unit="m^3 m-3",direction="-",&
                  surface_type=LIS_rc%lsm_index)

               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_DSSAT_LAI, &
                  value=dssat48_struc(n)%dssat48(t)%XLAI,&
                  vlevel=1,unit="m2/m2",direction="-",&
                  surface_type=LIS_rc%lsm_index)

               call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GWAD, &
                  value=dssat48_struc(n)%dssat48(t)%SDWT*10.0,&
                  vlevel=1,unit="kg/ha",direction="-",&
                  surface_type=LIS_rc%lsm_index)

               !-----------------------------------------------------------------------

                IF (YRDOY.EQ.YREND) THEN
                !-----------------------------------------------------------------------
                !     END of DAILY SIMULATION loop
                !----------------------------------------------------------------------
                !*********************************************************************** 
                !     End of Season 
                !*********************************************************************** 
                 dssat48_struc(n)%CONTROL(t)%dynamic = 6 !6: End Season

                 CALL LAND(dssat48_struc(n)%CONTROL(t), dssat48_struc(n)%ISWITCH(t), &
                    YRPLT, MDATE, YREND, n, t) !Pang: add n, t for ensembles and tiles

                 !JE Write variables to LIS_HIST file
                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD1, &
                    value=dssat48_struc(n)%dssat48(t)%SW(1),&
                    vlevel=1,unit="m^3 m-3",direction="-",&
                    surface_type=LIS_rc%lsm_index)
               
                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD2, &
                    value=dssat48_struc(n)%dssat48(t)%SW(2),&
                    vlevel=1,unit="m^3 m-3",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD3, &
                    value=dssat48_struc(n)%dssat48(t)%SW(3),&
                    vlevel=1,unit="m^3 m-3",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMD4, &
                    value=dssat48_struc(n)%dssat48(t)%SW(4),&
                    vlevel=1,unit="m^3 m-3",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_DSSAT_LAI, &
                    value=dssat48_struc(n)%dssat48(t)%XLAI,&
                    vlevel=1,unit="m2/m2",direction="-",&
                    surface_type=LIS_rc%lsm_index)

                 call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GWAD, &
                    value=dssat48_struc(n)%dssat48(t)%SDWT*10.0,&
                    vlevel=1,unit="kg/ha",direction="-",&
                    surface_type=LIS_rc%lsm_index)
              
                 dssat48_struc(n)%dssat48(t)%DONE = .FALSE. !Return back to INIT
                !dssat48_struc(n)%dssat48(t)%doseasinit = .TRUE.
                !Pang: DO THIS at the end of the season
                !YRSIM_SAVE = YRSIM !We don't really need to use YRSIM_SAVE
                
                !--- Check If Need to Go To INIT --------------------------------------
                     IF (YRDOY .GE. YRDOY_END) THEN
                         PRINT*, 'YRDOY IS GREATER THEN YRDOY_END'
                         EXIT
                         !REPNO = dssat48_struc(n)%CONTROL(t)%repno
                         !REPNO = REPNO+1
                         !dssat48_struc(n)%CONTROL(t)%repno = REPNO
                         !IF (REPNO.GT.NREPS) THEN
                         !   DONE= .TRUE.
                         !ELSE
                         !   RUN = 0
                         !ENDIF
                     ENDIF

                ENDIF
            ENDIF

            dssat48_struc(n)%dssat48(t)%yrend = YREND
            dssat48_struc(n)%dssat48(t)%mdate = MDATE
            dssat48_struc(n)%dssat48(t)%yrplt = YRPLT
 
            ! save state variables from local variables to global variables
            !CROCUS81_struc(n)%crocus81(t)%SNOWSWE(:)    = tmp_SNOWSWE(:)   
    
            ! save output variables from local variables to global variables
            !CROCUS81_struc(n)%crocus81(t)%THRUFAL      = tmp_THRUFAL     
            
            ![ 1] output variable: SNOWSWE (unit=kg/m2). *** Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
            !do i=1, CROCUS81_struc(n)%nsnow
            !    call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWLIQPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWSWE(i), &
            !                                      vlevel=i, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            !end do

            ! Reset forcing variables to zero
            dssat48_struc(n)%dssat48(t)%tair = 0.0
            dssat48_struc(n)%dssat48(t)%tmax = 0.0
            dssat48_struc(n)%dssat48(t)%tmin = 0.0
            dssat48_struc(n)%dssat48(t)%qair = 0.0
            dssat48_struc(n)%dssat48(t)%swdown = 0.0
            dssat48_struc(n)%dssat48(t)%lwdown = 0.0
            dssat48_struc(n)%dssat48(t)%uwind = 0.0
            dssat48_struc(n)%dssat48(t)%vwind = 0.0
            dssat48_struc(n)%dssat48(t)%wndspd = 0.0
            dssat48_struc(n)%dssat48(t)%psurf = 0.0
            dssat48_struc(n)%dssat48(t)%rainf = 0.0
            dssat48_struc(n)%dssat48(t)%snowf = 0.0
            dssat48_struc(n)%dssat48(t)%totprc = 0.0
            dssat48_struc(n)%dssat48(t)%tdew = 0.0

            if (dssat48_struc(n)%sm_coupling.eq.1) then
               do l=0, LIS_sfmodel_struc(n)%nsm_layers
                  dssat48_struc(n)%dssat48(t)%LIS_sm(l) = 0.0
                end do
           endif

        enddo ! end of tile (t) loop

        ! Reset forcing counter to be zero
        dssat48_struc(n)%forc_count = 0 

    endif ! end of alarmCheck loop 

    !Deallocate Memory for Local Variable
    !deallocate( tmp_SNOWSWE )
   
end subroutine dssat48_main
