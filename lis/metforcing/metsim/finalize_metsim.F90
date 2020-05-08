!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <LIS_misc.h>
!BOP
! !ROUTINE: finalize_metsim
! \label{finalize_metsim}
! 
! !REVISION HISTORY: 
! 29 Apr 2020: Zhuo Wang: updated for MetSim forcing data (based on princeton readers)
! 
! !INTERFACE:
subroutine finalize_metsim(findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use metsim_forcingMod, only : metsim_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for metsim forcing. 
!
!EOP
  implicit none
  integer :: findex
  integer :: n 
  
  do n=1,LIS_rc%nnest
     if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
        deallocate(metsim_struc(n)%n111)
        deallocate(metsim_struc(n)%n121)
        deallocate(metsim_struc(n)%n211)
        deallocate(metsim_struc(n)%n221)
        deallocate(metsim_struc(n)%w111)
        deallocate(metsim_struc(n)%w121)
        deallocate(metsim_struc(n)%w211)
        deallocate(metsim_struc(n)%w221)
     elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
        deallocate(metsim_struc(n)%n111)
        deallocate(metsim_struc(n)%n121)
        deallocate(metsim_struc(n)%n211)
        deallocate(metsim_struc(n)%n221)
        deallocate(metsim_struc(n)%w111)
        deallocate(metsim_struc(n)%w121)
        deallocate(metsim_struc(n)%w211)
        deallocate(metsim_struc(n)%w221)

        deallocate(metsim_struc(n)%n112)
        deallocate(metsim_struc(n)%n122)
        deallocate(metsim_struc(n)%n212)
        deallocate(metsim_struc(n)%n222)
        deallocate(metsim_struc(n)%w112)
        deallocate(metsim_struc(n)%w122)
        deallocate(metsim_struc(n)%w212)
        deallocate(metsim_struc(n)%w222)
     endif
  enddo
  deallocate(metsim_struc)
end subroutine finalize_metsim
