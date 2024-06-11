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
! !ROUTINE: dssat48_f2t
! \label{dssat48_f2t}
!
! !REVISION HISTORY:
!  26 Jun 2023: Pang-Wei liu; Initial Code
!  29 Aug 2023: J. Erlingis;  Compute Max/Min and Averages
! !INTERFACE:

subroutine dssat48_f2t(n)

! !USES:      
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_surface
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_logMod,        only : LIS_logunit, LIS_verify, &
                                LIS_endrun
  use dssat48_lsmMod
  use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc !JE for exchange
  use NoahMP401_lsmMod !JE for soil moisture exchange

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the DSSAT48
!  model tiles. DSSAT requires solar radiation, daily maximum and
!  minimum air temperatures, daily precipitation. Daily mean dewpoint
!  temperature and wind speed are preferred in order to calculate
!  evapotranspiration. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

  integer            :: t,v,l,status
  integer            :: tid
  real               :: ee, val, td

  ! Near Surface Air Temperature [K]
  type(ESMF_Field)  :: tmpField
  real, pointer     :: tmp(:)

  ! Near Surface Specific Humidity [kg/kg]
  type(ESMF_Field)  :: q2Field
  real, pointer     :: q2(:)

  ! Eastward Wind (u) [m/s]
  type(ESMF_Field)  :: uField
  real, pointer     :: uwind(:)

  ! Northward Wind (v) [m/s]
  type(ESMF_Field)  :: vField
  real, pointer     :: vwind(:)

  ! Rainfall Rate [kg/(m2s)]
  type(ESMF_Field)  :: pcpField
  real, pointer     :: pcp(:)

  ! Snowfall Rate [kg/(m2s)]
  type(ESMF_Field)  :: snowField
  real, pointer     :: snowf(:)

  ! Incident Longwave Radiation [W/m2]
  type(ESMF_Field)  :: lwdField
  real, pointer     :: lwd(:)

  ! Incident Shortwave Radiation [W/m2]
  type(ESMF_Field)  :: swdField
  real, pointer     :: swd(:)

  ! Surface Pressure [Pa]
  type(ESMF_Field)  :: psurfField
  real, pointer     :: psurf(:)

  ! Dewpoint Temperature [K]
  real, pointer     :: tdew(:)

  ! Wind Speed [m/s]
  real, pointer     :: wndspd(:)

  integer, pointer   :: layer_windht(:), layer_relhumht(:)
! ____________________________________________________________________________________

  ! If forcing heights are not specified, then LIS will assume that forcing data
  ! corresponds to the reference heights snowmodel.par (tmprh_ht, wind_ht).

  ! Get near surface air temperature
 
  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Tair%varname(1)),tmpField,rc=status)
  call LIS_verify(status,'dssat48_f2t: error getting Tair')

  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr=tmp,rc=status)
  call LIS_verify(status, 'dssat48_f2t: error retrieving Tair')

  ! Get near surface specific humidity
  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Qair%varname(1)),q2Field,rc=status)
  call LIS_verify(status,'dssat48_f2t: error getting Qair')

  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status,'dssat48_f2t: error retrieving q2')

  ! Get eastward wind (u)
  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_E%varname(1)),uField,rc=status)
  call LIS_verify(status,'dssat48_f2t: error getting Wind_E')

  call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
  call LIS_verify(status,'dssat48_f2t: error retrieving u')

  ! Get northward wind (v)
  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_N%varname(1)),vField,rc=status)
  call LIS_verify(status,'dssat48_f2t: error getting Wind_N')

  call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
  call LIS_verify(status,'dssat48_f2t: error retrieving v')

  ! Get incident shortwave radiation
  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_SWdown%varname(1)),swdField,rc=status)
  call LIS_verify(status,'dssat48_f2t: error getting SWdown')

  call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
  call LIS_verify(status,'dssat48_f2t: error retrieving swd')

  ! Get incident longwave radiation
  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_LWdown%varname(1)),lwdField,rc=status)
  call LIS_verify(status,'dssat48_f2t: error getting LWdown')

  call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
  call LIS_verify(status,'dssat48_f2t: error retrieving lwd')

  ! Get surface pressure
  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Psurf%varname(1)),psurfField,rc=status)
  call LIS_verify(status, 'dssat48_f2t: error getting PSurf')

  call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status,'dssat48_f2t: error retrieving psurf')

  ! Get rainfall rate
  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Rainf%varname(1)),pcpField,rc=status)
  call LIS_verify(status,'dssat48_f2t: error getting Rainf')

  call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status,'dssat48_f2t: error retrieving pcp')

  ! Get snowfall rate
  if(LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Snowf%varname(1)),snowField,rc=status)
     call LIS_verify(status,'dssat48_f2t: error getting Snowf')

     call ESMF_FieldGet(snowField,localDE=0, farrayPtr=snowf,rc=status)
     call LIS_verify(status,'dssat48_f2t: error retrieving snowf')
  endif

  !Keep track of number of forcing times so that we can compute daily averages

! JE This is a test for keeping end of day soil moisture from the previous day
!  if (dssat48_struc(n)%forc_count.eq.0) then !Remember ending soil moisture from yesterday
!     write(LIS_logunit,*) "Saving LSM soil moisture "
!     do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
!         if (dssat48_struc(n)%sm_coupling.eq.1) then
!            do l=0, LIS_sfmodel_struc(n)%nsm_layers
!               dssat48_struc(n)%dssat48(t)%LIS_sm(l) = NOAHMP401_struc(n)%noahmp401(t)%smc(l)
!            end do
!         endif
!      end do
!  end if

  dssat48_struc(n)%forc_count = dssat48_struc(n)%forc_count + 1

  do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)

     ! write(LIS_logunit,*) "At tile ",t
     ! Transform tile to the patch
     tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id

     ! Soil moisture
     if (dssat48_struc(n)%sm_coupling.eq.1) then
        do l=0, LIS_sfmodel_struc(n)%nsm_layers
           dssat48_struc(n)%dssat48(t)%LIS_sm(l) = dssat48_struc(n)%dssat48(t)%LIS_sm(l) &
     &        + NOAHMP401_struc(n)%noahmp401(t)%smc(l)
         end do
     endif

     ! Air temperature (DSSAT requires daily maximum and minimum)
     dssat48_struc(n)%dssat48(t)%tair=dssat48_struc(n)%dssat48(t)%tair + tmp(tid)
 
     if (dssat48_struc(n)%forc_count.eq.1) then !First iteration set max/min 
        !write(LIS_logunit,*) "First loop temperature ", tmp(tid)
        dssat48_struc(n)%dssat48(t)%tmax = tmp(tid)
        dssat48_struc(n)%dssat48(t)%tmin = tmp(tid)
     else
        !write(LIS_logunit,*) "Loop ",dssat48_struc(n)%forc_count
        if (tmp(tid).gt.dssat48_struc(n)%dssat48(t)%tmax) then
           !write(LIS_logunit,*) "Maximum temperature replaced ",tmp(tid)
           dssat48_struc(n)%dssat48(t)%tmax=tmp(tid) !Replace maximum temperature
        endif
        if (tmp(tid).lt.dssat48_struc(n)%dssat48(t)%tmin) then
          !write(LIS_logunit,*) "Minimum temperature replaced ", tmp(tid)
          dssat48_struc(n)%dssat48(t)%tmin=tmp(tid) !Replace minimum temperature
        endif
     endif !Check loop index     

     ! Specific Humidity
     dssat48_struc(n)%dssat48(t)%qair=dssat48_struc(n)%dssat48(t)%qair + q2(tid)

     ! Shortwave Radiation
     dssat48_struc(n)%dssat48(t)%swdown=dssat48_struc(n)%dssat48(t)%swdown + swd(tid)

     ! Longwave Radiation
     dssat48_struc(n)%dssat48(t)%lwdown=dssat48_struc(n)%dssat48(t)%lwdown + lwd(tid)

     ! Wind_E (u)
     dssat48_struc(n)%dssat48(t)%uwind=dssat48_struc(n)%dssat48(t)%uwind + uwind(tid)

     ! Wind_N (v)
     dssat48_struc(n)%dssat48(t)%vwind=dssat48_struc(n)%dssat48(t)%vwind + vwind(tid)

     ! Calculate Magnitude of Wind Speed (m/s) 
     dssat48_struc(n)%dssat48(t)%wndspd = dssat48_struc(n)%dssat48(t)%wndspd + SQRT(uwind(tid)**2 + vwind(tid)**2)

     ! Surface Pressure
     dssat48_struc(n)%dssat48(t)%psurf=dssat48_struc(n)%dssat48(t)%psurf + psurf(tid)

     ! Rainfall 
     dssat48_struc(n)%dssat48(t)%rainf=dssat48_struc(n)%dssat48(t)%rainf + pcp(tid)

     if(LIS_FORC_Snowf%selectOpt.eq.1) then
        ! Snowfall
        dssat48_struc(n)%dssat48(t)%snowf= dssat48_struc(n)%dssat48(t)%snowf + snowf(tid)
        ! Total Precipitation
        dssat48_struc(n)%dssat48(t)%totprc = dssat48_struc(n)%dssat48(t)%totprc + pcp(tid) + snowf(tid)
     else
        ! Total Precipitation
        dssat48_struc(n)%dssat48(t)%totprc = dssat48_struc(n)%dssat48(t)%totprc + pcp(tid)
     endif

     ! Calculate Dewpoint

     ! Following A First Course in Atmospheric Thermodynamics, assume
     ! approximation q = epsilon*e/p

     ! Calculate vapor pressure
     ee = (q2(tid)*psurf(tid))/0.622

     ! Invert Bolton 1980 formula for saturation vapor pressure to calculate Td
     ! since es(Td) = e

     val = log(ee/611.2)
     td = (243.5 * val) / (17.67 - val) ! Dewpoint in C
     td = td + 273.15
     !write(LIS_logunit,*) 'Calculating dewpoint from T, p, and q2'
     !write(LIS_logunit,*) 'T: ',tmp(tid), 'q2: ',q2(tid),'p: ',psurf(tid),'Td: ',td
     dssat48_struc(n)%dssat48(t)%tdew = dssat48_struc(n)%dssat48(t)%tdew + td 
     

  enddo

end subroutine dssat48_f2t
