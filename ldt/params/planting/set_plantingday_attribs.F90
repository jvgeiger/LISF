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
! !ROUTINE: set_plantingday_attribs
!  \label{set_plantingday_attribs}
!
! !REVISION HISTORY:
!  28 Feb 2025: Pang-Wei Liu; Initial Specification!
!
! !INTERFACE:
subroutine set_plantingday_attribs()

! !USES:
  use LDT_paramDataMod
  use LDT_plantingdayMod

  implicit none

! !ARGUMENTS:

! !DESCRIPTION:
!  Set planting day attributes.
!EOP
!
  LDT_plantingday_struc(:)%plantingday%num_bins = 1
  LDT_plantingday_struc(:)%plantingday%vlevels = 1
end subroutine set_plantingday_attribs
