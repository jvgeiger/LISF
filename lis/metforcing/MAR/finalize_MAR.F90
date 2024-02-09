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
! !MODULE: finalize_MAR
!  \label{finalize_MAR}
!
! !REVISION HISTORY: 
! 
! !INTERFACE:
subroutine finalize_MAR(findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use MAR_forcingMod, only : MAR_struc
!
! !DESCRIPTION:
!  Routine to cleanup MAR forcing related memory allocations.   
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing
!  \end{description}
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
  deallocate(MAR_struc)

end subroutine finalize_MAR
