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
! !ROUTINE: dssat48_getLSMexport
! \label{dssat48_getLSMexport}
!
! !REVISION HISTORY:
! 27 Jun 2023: Pang-Wei Liu
!
! !INTERFACE:
subroutine dssat48_getLSMexport(n, SubLSM2LSM_State)
! !USES:

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use dssat48_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: SubLSM2LSM_State
! 
! !DESCRIPTION:
! 
! 
!EOP
  type(ESMF_Field)   :: snwdField, sweField
  real, pointer      :: swe(:), snwd(:)
  integer            :: t
  integer            :: status
  real               :: tmp_SLOPE
 
  !call ESMF_StateGet(SubLSM2LSM_State,"Total SWE",sweField,rc=status)
  !call LIS_verify(status,'Snowmodel_getLSMexport: SM SWE state error get')
  !call ESMF_StateGet(SubLSM2LSM_State,"Total snowdepth",snwdField,rc=status)
  !call LIS_verify(status,'Snowmodel_getLSMexport: SM snowd state error get')

  !call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  !call LIS_verify(status,'Snowmodel_getLSMexport: SM swe field error get')
  !call ESMF_FieldGet(snwdField,localDE=0,farrayPtr=snwd,rc=status)
  !call LIS_verify(status,'Snowmodel_getLSMexport: SM snwd field error get')

end subroutine dssat48_getLSMexport


