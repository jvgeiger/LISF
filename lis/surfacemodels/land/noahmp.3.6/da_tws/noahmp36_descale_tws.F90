!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp36_descale_tws
! \label{noahmp36_descale_tws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_descale_tws(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noahmp36_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Descales twsoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP

!Wanshu

  type(ESMF_Field)     :: gwField
  real, pointer        :: gws(:)
  type(ESMF_Field)     :: gwIncrField
  real, pointer        :: gwsIncr(:)
  integer              :: t
  integer              :: status

#if 0 
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp36_settws")

  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Groundwater Storage failed in noahmp36_settws")

  call ESMF_StateGet(LSM_Incr_State, "Groundwater Storage",gwIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp36_descale")

  call ESMF_FieldGet(gwIncrField, localDE=0, farrayPtr=gwsIncr,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp36_descale")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gws(t)   = gws(t)*10.0
     gwsIncr(t) = gwsIncr(t)*10.0
  enddo
#endif
 

end subroutine noahmp36_descale_tws

