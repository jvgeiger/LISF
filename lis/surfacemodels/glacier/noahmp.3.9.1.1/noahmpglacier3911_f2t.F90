!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmpglacier3911_f2t
! \label{noahmpglacier3911_f2t}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_f2t(n)
! !USES:
    use ESMF
    use LIS_coreMod, only       : LIS_rc , LIS_surface
    use LIS_metforcingMod, only : LIS_FORC_State
    use LIS_logMod, only        : LIS_verify
    use LIS_FORC_AttributesMod   
    use noahmpglacier3911_Mod

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n     
!
! !DESCRIPTION:
!  This routine transfers the LIS provided forcing into the noahmpglacier3911
!  model tiles.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP
    
    integer           :: t, v, status
    integer           :: tid 
 
    ! Near Surface Air Temperature [K]
    type(ESMF_Field)  :: tmpField
    real, pointer     :: tmp(:)
 
    ! Near Surface Specific Humidity [kg kg-1]
    type(ESMF_Field)  :: q2Field
    real, pointer     :: q2(:)
 
    ! Incident Shortwave Radiation [W m-2]
    type(ESMF_Field)  :: swdField
    real, pointer     :: swd(:)
 
    ! Incident Longwave Radiation [W m-2]
    type(ESMF_Field)  :: lwdField
    real, pointer     :: lwd(:)
 
    ! Eastward Wind [W m-2]
    type(ESMF_Field)  :: uField
    real, pointer     :: uwind(:)
 
    ! Northward Wind [m s-1]
    type(ESMF_Field)  :: vField
    real, pointer     :: vwind(:)
 
    ! Surface Pressure [Pa]
    type(ESMF_Field)  :: psurfField
    real, pointer     :: psurf(:)
 
    ! Rainfall Rate [kg m-2 s-1]
    type(ESMF_Field)  :: pcpField
    real, pointer     :: pcp(:)
 
    ! Snowfall Rate [kg m-2 s-1]
    type(ESMF_Field)  :: snowField
    real, pointer     :: snowf(:)

    integer, pointer   :: layer_h(:), layer_m(:)
 
    !!! GET FORCING FIELDS FROM LIS
    ! get near surface air temperature
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(1)), tmpField, rc=status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error getting Tair")

    call ESMF_FieldGet(tmpField, localDE = 0, farrayPtr = tmp, rc = status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Tair")
    
    ! get near surface specific humidity
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Qair%varname(1)), q2Field, rc=status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error getting Qair")

    call ESMF_FieldGet(q2Field, localDE = 0, farrayPtr = q2, rc = status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Qair")
    
    ! get incident shortwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_SWdown%varname(1)), swdField, rc=status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error getting Swdown")

    call ESMF_FieldGet(swdField, localDE = 0, farrayPtr = swd, rc = status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Swdown")
    
    ! get incident longwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_LWdown%varname(1)), lwdField, rc=status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error getting Lwdown")

    call ESMF_FieldGet(lwdField, localDE = 0, farrayPtr = lwd, rc = status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Lwdown")
    
    ! get eastward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_E%varname(1)), uField, rc=status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error getting Wind_E")

    call ESMF_FieldGet(uField, localDE = 0, farrayPtr = uwind, rc = status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Wind_E")
    
    ! get northward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_N%varname(1)), vField, rc=status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error getting Wind_N")

    call ESMF_FieldGet(vField, localDE = 0, farrayPtr = vwind, rc = status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Wind_N")
    
    ! get surface pressure
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Psurf%varname(1)), psurfField, rc=status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error getting Psurf")

    call ESMF_FieldGet(psurfField, localDE = 0, farrayPtr = psurf, rc = status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Psurf")
    
    ! get rainfall rate
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Rainf%varname(1)), pcpField, rc=status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error getting Rainf")

    call ESMF_FieldGet(pcpField, localDE = 0, farrayPtr = pcp, rc = status)
    call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Rainf")
    
    ! get snowfall rate
    if(LIS_Forc_Snowf%selectOpt .eq. 1) then 
        call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Snowf%varname(1)), snowField, rc=status)
        call LIS_verify(status, "noahmpglacier3911_f2t: error getting Snowf")

        call ESMF_FieldGet(snowField, localDE = 0, farrayPtr = snowf, rc = status)
        call LIS_verify(status, "noahmpglacier3911_f2t: error retrieving Snowf")
    endif 

    !!! set the forcing counter
    Noahmpgl3911_struc(n)%forc_count = Noahmpgl3911_struc(n)%forc_count + 1
 
    !!! pass forcing data to tiles
    do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
        tid = LIS_surface(n, LIS_rc%glacier_index)%tile(t)%tile_id

        ! TAIR
        Noahmpgl3911_struc(n)%noahmpgl(t)%tair = Noahmpgl3911_struc(n)%noahmpgl(t)%tair + tmp(tid)

        ! QAIR
        Noahmpgl3911_struc(n)%noahmpgl(t)%qair = Noahmpgl3911_struc(n)%noahmpgl(t)%qair + q2(tid)

        ! SWDOWN
        Noahmpgl3911_struc(n)%noahmpgl(t)%swdown = Noahmpgl3911_struc(n)%noahmpgl(t)%swdown + swd(tid)

        ! LWDOWN
        Noahmpgl3911_struc(n)%noahmpgl(t)%lwdown = Noahmpgl3911_struc(n)%noahmpgl(t)%lwdown + lwd(tid)

        ! WIND_E
        Noahmpgl3911_struc(n)%noahmpgl(t)%wind_e = Noahmpgl3911_struc(n)%noahmpgl(t)%wind_e + uwind(tid)

        ! WIND_N
        Noahmpgl3911_struc(n)%noahmpgl(t)%wind_n = Noahmpgl3911_struc(n)%noahmpgl(t)%wind_n + vwind(tid)

        ! PSURF
        Noahmpgl3911_struc(n)%noahmpgl(t)%psurf = Noahmpgl3911_struc(n)%noahmpgl(t)%psurf + psurf(tid)

        ! RAINF
        if(pcp(tid).ne.LIS_rc%udef) then
           Noahmpgl3911_struc(n)%noahmpgl(t)%prcp = Noahmpgl3911_struc(n)%noahmpgl(t)%prcp + pcp(tid)
        endif

        ! SNOWF
        ! If there is snowf add it to precipitation.  NoahMP-3.6 does not use
        ! separate rainf and snowf.  It determines what to do with precipitation.
        if(LIS_Forc_Snowf%selectOpt .eq. 1) then 
           if(snowf(tid).ne.LIS_rc%udef) then
              Noahmpgl3911_struc(n)%noahmpgl(t)%prcp = Noahmpgl3911_struc(n)%noahmpgl(t)%prcp + snowf(tid)
           endif
        endif
    enddo

end subroutine noahmpglacier3911_f2t
