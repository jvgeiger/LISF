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
    use dssat48_lsmMod
    USE ModuleDefs, only: ControlType,SwitchType !From DSSAT
   !use other modules
  
    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    real                 :: dt
    real                 :: lat, lon
    integer              :: row, col
    integer              :: year, month, day, hour, minute, second
    logical              :: alarmCheck
    integer              :: nyrs_lis, yrdoy_end_lis 

    real                 :: tmp_LAT                ! Latitude in decimal degree  (latitude (degrees +North)) [degrees]
    real                 :: tmp_LON                ! Longitude in decimal year    (longitude (degrees +East)) [degrees]
    integer              :: tmp_year, tmp_month, tmp_day, tmp_hour, tmp_minute
    real                 :: tmp_pres, tmp_precip, tmp_tmax, tmp_tmin, tmp_tdew   ! Weather Forcing
    real                 :: tmp_swrad, tmp_wind                                  ! Weather Forcing
    character*3          :: fnest ! MN Bug in toolkit (added to this code)
! CONTRO
    integer             :: RUN, REPNO, TRTNUM, ROTNUM
    CHARACTER*30        :: FILEIO
    CHARACTER*12        :: FILEX   !,DSCSM,INPUT
    CHARACTER*1         :: RNMODE
!
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
    PRINT*, 'Im Here in dssat48_main'
     write(fnest,'(i3.3)') n
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "DSSAT48 model alarm "// trim(fnest)) !MN  Bug in the toolkit 
    if (alarmCheck) Then
        !do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        do t = 1, 1
            PRINT*, 'Am I inside the alarm on loop: Yes'
             
            dt = LIS_rc%ts
            PRINT*, "t, dt: ", t, dt
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

            !!------ This Block Is Where We Obtain Weather Forcing ------------------------------!!
            ! retrieve forcing data from DSSAT48_struc(n)%dssat(t) and assign to local variables

            write(LIS_logunit,*) 'Weather forcing for tile: ',t

            ! Daily average surface pressure (Pa)
            tmp_pres      = dssat48_struc(n)%dssat48(t)%psurf / dssat48_struc(n)%forc_count
            write(LIS_logunit,*) 'p: ',tmp_pres
 
            ! Total daily precipitation (rain+snow) (mm)
            tmp_precip    = dssat48_struc(n)%dssat48(t)%totprc * 3600. * 24. !Convert from kg/ms2 to mm
            write(LIS_logunit,*) 'Precip: ',tmp_precip

            ! Tmax: maximum daily air temperature (C)
            tmp_tmax      = dssat48_struc(n)%dssat48(t)%tmax - 273.15 !Convert from K to C
            write(LIS_logunit,*) 'Tmax: ',tmp_tmax

            ! Tmin: minimum daily air temperature (C)
            tmp_tmin      = dssat48_struc(n)%dssat48(t)%tmin - 273.15 !Convert from K to C 
            write(LIS_logunit,*) 'Tmin: ',tmp_tmin

            ! Tdew: average daily dewpoint temperature (C)
            tmp_tdew      = (dssat48_struc(n)%dssat48(t)%tdew / dssat48_struc(n)%forc_count) - 273.15 !Convert from K to C
            write(LIS_logunit,*) 'Tdew: ',tmp_tdew

            ! SW_RAD: daily total incoming solar radiation (MJ/(m2d))
            tmp_swrad     = (dssat48_struc(n)%dssat48(t)%swdown / dssat48_struc(n)%forc_count) * 0.0864 !Convert from W/m2 to MJ/(m2d)
            write(LIS_logunit,*) 'Swrad: ',tmp_swrad

            ! Wind: daily average wind speed (km/d)
            tmp_wind      = (dssat48_struc(n)%dssat48(t)%wndspd / dssat48_struc(n)%forc_count) * 86.0 !Convert from m/s to km/d
            write(LIS_logunit,*) 'Wind: ',tmp_wind

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
            PRINT*, 'time: ', tmp_year, tmp_month, tmp_day, tmp_hour, tmp_minute
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

!------------------------------------------------------------------------------------ 
! call model physics
            !INIFOR CONTROL
            !-------------------------------------------
            !RUN   = 0
            !REPNO = 1
            !CONTROL % REPNO = REPNO
            !-------------------------------------------
            !PRINT*, 'CONTROL_main_1: ', CONTROL
            !PRINT*, 'SWITCH_main_1: ', ISWITCH
            !FILEIO = 'DSSAT48.INP'
            !FILEX = 'NASA2019.SQX'
            !RNMODE = 'Q'
            !ROTNUM = 1
            !TRTNUM = 1
            !CONTROL % FILEIO = FILEIO
            !CONTROL % RNMODE = RNMODE
            !CONTROL % FILEX   = FILEX
            !CONTROL % ROTNUM  = ROTNUM !Pang: 1, SQ in run.v48
            !CONTROL % TRTNUM  = TRTNUM !Pang: 1, TRTNO in run.v48
            !CONTROL % ERRCODE = 0      ! Default
             PRINT*, 'yrend ini: ', dssat48_struc(n)%dssat48(t)%yrend
             PRINT*, 'CONTROL_main insideloop: ', dssat48_struc(n)%CONTROL(t)%run
             PRINT*, 'CONTROL_lsmMod2 fileio: ', dssat48_struc(n)%CONTROL(t)%fileio
            !PRINT*, 'CONTROL_main_2: ', CONTROL
            !PRINT*, 'SWITCH_main_2: ', ISWITCH
             nyrs_lis = dssat48_struc(n)%dssat48(t)%nyrs
             yrdoy_end_lis = dssat48_struc(n)%dssat48(t)%yrdoy_end
             PRINT*, 'nyrs_llis ', nyrs_lis, ' at t= ', t, ' and n= ', n
             PRINT*, 'nyrs_llis ', nyrs_lis, ' at t= ', t, ' and n= ', n  
            call CSM_RUN(dssat48_struc(n)%CONTROL(t), dssat48_struc(n)%ISWITCH(t), dssat48_struc(n)%dssat48(t)%yrend, dssat48_struc(n)%dssat48(t)%expno, &
                         dssat48_struc(n)%dssat48(t)%trtall,  dssat48_struc(n)%dssat48(t)%nreps, &
                         dssat48_struc(n)%dssat48(t)%endyrs, yrdoy_end_lis )
 
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
        enddo ! end of tile (t) loop

        ! Reset forcing counter to be zero
        dssat48_struc(n)%forc_count = 0 

    endif ! end of alarmCheck loop 

    !Deallocate Memory for Local Variable
    !deallocate( tmp_SNOWSWE )
   
end subroutine dssat48_main
