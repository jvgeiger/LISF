!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: reset_MAR
!  \label{reset_MAR}
!
! !REVISION HISTORY: 
! 
! !INTERFACE:
subroutine reset_MAR()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use MAR_forcingMod, only : MAR_struc
!
! !DESCRIPTION:
!  Routine to reset MAR forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n

  do n=1,LIS_rc%nnest
     MAR_struc(n)%MARtime1 = 3000.0
     MAR_struc(n)%MARtime2 = 0.0
  enddo

end subroutine reset_MAR
