!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: reset_metsim
!  \label{reset_metsim}
!
! !REVISION HISTORY: 
! 30 Apr 2020: Zhuo Wang, updated for MetSim forcing data (based on princeton) 
! 
! !INTERFACE:
subroutine reset_metsim()
! !USES:
  use LIS_coreMod,          only : LIS_rc
  use metsim_forcingMod, only : metsim_struc
!
! !DESCRIPTION:
!  Routine to reset metsim forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     metsim_struc(n)%metsimtime1 = 3000.0
     metsim_struc(n)%metsimtime2 = 0.0
  enddo

end subroutine reset_metsim
