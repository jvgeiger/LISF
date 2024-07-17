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
! !ROUTINE: dssat48_readcrd
! \label{dssat48\_readcrd}
!
! !REVISION HISTORY:
!  11 May 2023: Pang-Wei Liu
!   
! !INTERFACE:
subroutine dssat48_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify, LIS_endrun
    use dssat48_lsmMod, only       : dssat48_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to dssat48 model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    !character*6  :: str_i
 
    write(LIS_logunit, *) "Start reading LIS configuration file for DSSAT48 model"
    
    call ESMF_ConfigFindLabel(LIS_config, "DSSAT48 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "DSSAT48 model timestep: not defined")
        call LIS_parseTimeString(time, dssat48_struc(n)%ts)
    enddo
 
    call ESMF_ConfigFindLabel(LIS_config, "DSSAT48 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"DSSAT48 restart output interval: not defined")
        call LIS_parseTimeString(time, dssat48_struc(n)%rstInterval)
    enddo

    !PL DSSAT Restart-------------------------------------------------------------------
    call ESMF_ConfigFindLabel(LIS_config, "DSSAT48 start mode:", rc = rc)
       do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, dssat48_struc(n)%dssatstartcode, rc=rc)
            call LIS_verify(rc, "DSSAT48 start mode: not defined")
            if (dssat48_struc(n)%dssatstartcode == "restart" ) then
                call ESMF_ConfigFindLabel(LIS_config,"DSSAT48 restart file:",rc=rc)
                call ESMF_ConfigGetAttribute(LIS_config,dssat48_struc(n)%rfile,rc=rc)
                call LIS_verify(rc,'DSSAT48 restart file: not defined')
            endif
       enddo

   
    !if (trim(dssat48_struc(n)%dssatstartcode) == "restart" ) then
    !    call ESMF_ConfigFindLabel(LIS_config,"DSSAT48 restart file:",rc=rc)
    !   do n=1,LIS_rc%nnest
    !      call ESMF_ConfigGetAttribute(LIS_config,dssat48_struc(n)%rfile,rc=rc)
    !      call LIS_verify(rc,'DSSAT48 restart file: not defined')
    !  enddo
    !endif

    Call ESMF_ConfigFindLabel(LIS_config, "DSSAT48 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, dssat48_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "DSSAT48 restart file format: not defined")
        enddo
     !--------------------------------------------------------------------------------
    call ESMF_ConfigFindLabel(LIS_config, "DSSAT48 experiment file:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, dssat48_struc(n)%expfile, rc=rc)
            call LIS_verify(rc, "DSSAT48 experiment file: not defined")
        enddo

    !JE Coupling Options 
    call ESMF_ConfigFindLabel(LIS_config, "DSSAT48 soil moisture coupling:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, dssat48_struc(n)%sm_coupling, rc = rc)
        call LIS_verify(rc,"DSSAT48 soil moisture coupling: not defined")
        if ((dssat48_struc(n)%sm_coupling.ne.0).AND.(dssat48_struc(n)%sm_coupling.ne.1)) then
           write(LIS_logunit,*) "[ERR] Valid options for DSSAT soil moisture coupling are 0=No or 1=Yes"
           call LIS_endrun()
        endif
        write(LIS_logunit,*) "SM Coupling Flag ", dssat48_struc(n)%sm_coupling
        if (dssat48_struc(n)%sm_coupling .eq. 1) then
          write(LIS_logunit,*) "LIS - DSSAT Soil Moisture Coupling ON"
        else
         write(LIS_logunit,*) "LIS - DSSAT Soil Moisture Coupling OFF"
        endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config, "DSSAT48 LAI coupling:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, dssat48_struc(n)%lai_coupling, rc = rc)
        call LIS_verify(rc,"DSSAT48 LAI coupling: not defined")
        if ((dssat48_struc(n)%lai_coupling.ne.0).AND.(dssat48_struc(n)%lai_coupling.ne.1)) then
           write(LIS_logunit,*) "[ERR] Valid options for DSSAT LAI coupling are 0=No or 1=Yes"
           call LIS_endrun()
        endif
        write(LIS_logunit,*) "[INFO] LAI Coupling Flag ", dssat48_struc(n)%lai_coupling
        if (dssat48_struc(n)%lai_coupling .eq. 1) then
          write(LIS_logunit,*) "[INFO] LIS - DSSAT LAI Coupling ON"
        else
         write(LIS_logunit,*) "[INFO] LIS - DSSAT LAI Coupling OFF"
        endif
    enddo 

    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number years of spin up
!    call ESMF_ConfigFindLabel(LIS_config, "DSSAT48 Spin Up:", rc = rc)
!    do n=1, LIS_rc%nnest
!        call ESMF_ConfigGetAttribute(LIS_config, dssats_struc(n)%nspin, rc=rc)
!        call LIS_verify(rc, "DSSAT48 nspin: not defined")
!    enddo

     
    write(LIS_logunit, *) "Finish reading LIS configuration file for DSSAT48 model"
     
end subroutine dssat48_readcrd
