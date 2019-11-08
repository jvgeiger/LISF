!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: summa2_setup
! \label{summa2_setup}
!
! !REVISION HISTORY:
!
! !INTERFACE:
subroutine summa2_setup()
! !USES:
   use nrtype
   use summa2_lsmMod
   use LIS_coreMod,  only : LIS_rc
   use summa_util,   only : handle_err
   use summa_setup,  only : summa_paramSetup

! !ARGUMENTS: 
!
! !DESCRIPTION:
!
!  Complete the setup routines for summa2 option (forcing-only) 
! 
!EOP

   implicit none

   integer             :: n
   integer(i4b)        :: err=0                      ! error code
   character(len=1024) :: message=''                 ! error message

   do n = 1, LIS_rc%nnest
      if ( summa1_struc(n)%nGRU > 0 ) then
         call summa_paramSetup(summa1_struc(n), err, message)
         call handle_err(err, message)
      endif
   enddo

end subroutine summa2_setup
