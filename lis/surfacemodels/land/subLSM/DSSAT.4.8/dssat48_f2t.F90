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
!
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

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the DSSAT48
!  model tiles. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

  integer            :: t,v,status
  integer            :: tid

  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,snowfField

  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:)

  integer, pointer   :: layer_windht(:), layer_relhumht(:)
! __________________

  ! If forcing heights are not specified, then LIS will assume that forcing data
  ! corresponds to the reference heights snowmodel.par (tmprh_ht, wind_ht).
  ! get near surface air temperature
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Tair%varname(1)),tmpField,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error getting Tair')
  !  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr=tmp,rc=status)
  !  call LIS_verify(status, 'dssat48_f2t: error retrieving Tair')

  ! get near surface specific humidity
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Qair%varname(1)),q2Field,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error getting Qair')
  !  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error retrieving q2')

  ! get eastward wind
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_E%varname(1)),uField,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error getting Wind_E')
  !  call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error retrieving u')


  ! get northward wind
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_N%varname(1)),vField,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error getting Wind_N')
  !  call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error retrieving v')

  ! get incident shortwave radiation
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_SWdown%varname(1)),swdField,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error getting SWdown')
  !  call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error retrieving swd')

  ! get incident longwave radiation
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_LWdown%varname(1)),lwdField,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error getting LWdown')
  !  call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error retrieving lwd')

  ! get surface pressure
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Psurf%varname(1)),psurfField,rc=status)
  !  call LIS_verify(status, 'dssat48_f2t: error getting PSurf')
  !  call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error retrieving psurf')

  ! get rainfall rate
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Rainf%varname(1)),pcpField,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error getting Rainf')
  !  call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error retrieving pcp')

  ! get snowfall rate
  !  call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Snowf%varname(1)),snowfField,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error getting Snowf')
  !  call ESMF_FieldGet(snowfField,localDE=0, farrayPtr=snowf,rc=status)
  !  call LIS_verify(status,'dssat48_f2t: error retrieving snowf')



   ! dssat48_struc(n)%forc_count = dssat48_struc(n)%forc_count + 1

 ! do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
 !    ! Transform tile to the patch
 !    tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
 !    ! TA
 !    dssat48_struc(n)%dssat48(t)%tair=dssat48_struc(n)%dssat48(t)%tair + tmp(tid)
 !    ! QA
 !    !dssat48_struc(n)%dssat48(t)%qair=dssat48_struc(n)%dssat48(t)%qair + q2(tid)
 !    ! SW_RAD
 !    dssat48_struc(n)%dssat48(t)%swdown=dssat48_struc(n)%dssat48(t)%swdown + swd(tid)
 !    !LW_RAD
 !    dssat48_struc(n)%dssat48(t)%lwdown=dssat48_struc(n)%dssat48(t)%lwdown + lwd(tid)
 !    !Wind_E
 !    dssat48_struc(n)%dssat48(t)%uwind=dssat48_struc(n)%dssat48(t)%uwind + uwind(tid)
 !    !WInd_N
 !    dssat48_struc(n)%dssat48(t)%vwind=dssat48_struc(n)%dssat48(t)%vwind + vwind(tid)
 !    ! PPS
 !    dssat48_struc(n)%dssat48(t)%psurf=dssat48_struc(n)%dssat48(t)%psurf + psurf(tid)
 !    ! RRSNOW
 !    dssat48_struc(n)%dssat48(t)%rainf=dssat48_struc(n)%dssat48(t)%rainf + pcp(tid)
 !    ! SRSNOW
 !    dssat48_struc(n)%dssat48(t)%snowf=dssat48_struc(n)%dssat48(t)%snowf + snowf(tid)
 ! enddo

end subroutine dssat48_f2t
