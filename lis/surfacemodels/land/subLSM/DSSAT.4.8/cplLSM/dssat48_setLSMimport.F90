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
! !ROUTINE: dssat48_setLSMimport
! \label{dssat48_setLSMimport}
!
! !REVISION HISTORY:
! 26 Jun 2023: Pang-Wei Liu

! !INTERFACE:
subroutine dssat48_setLSMimport(n, LSM2SubLSM_State)
! !USES:

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use dssat48_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM2SubLSM_State
! 
! !DESCRIPTION:
! 
! 
!EOP

#if 0
  type(ESMF_Field)   :: gtField
  real, pointer      :: gt(:)
  type(ESMF_Field)   :: XWGIField
  real, pointer      :: XWGI(:)
  type(ESMF_Field)   :: XWGField
  real, pointer      :: XWG(:)
  integer            :: t
  integer            :: status

  call ESMF_StateGet(LSM2SubLSM_State,"Ground temperature",gtField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(gtField,localDE=0,farrayPtr=gt,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM2SUBLSM_State,"soil volumetric frozen water content",XWGIField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(XWGIField,localDE=0,farrayPtr=XWGI,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM2SUBLSM_State,"soil volumetric liquid water content",XWGField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(XWGField,localDE=0,farrayPtr=XWG,rc=status)
  call LIS_verify(status)


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     dssat48_struc(n)%dssat48(t)%TG = gt(t)
     dssat48_struc(n)%dssat48(t)%XWGI = XWGI(t)
     dssat48_struc(n)%dssat48(t)%XWG  = XWG(t)
  enddo

#endif

end subroutine dssat48_setLSMimport


