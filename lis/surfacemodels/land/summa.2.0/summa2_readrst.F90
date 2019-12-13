!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: summa2_readrst
! \label{summa2_readrst}
! 
! !INTERFACE:
subroutine summa2_readrst()

! !USES:
   use nrtype
   use LIS_coreMod,   only : LIS_rc
   use summa2_lsmMod
   use summa_util,    only : handle_err
   use summa_restart, only : summa_readRestart

   implicit none
! !ARGUMENTS: 

!
! !DESCRIPTION:
! 
!  Calls the restart reading routines for the forcing-only option (summa2)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

   integer             :: n
   integer(i4b)        :: err=0                      ! error code
   character(len=1024) :: message=''                 ! error message

   do n = 1, LIS_rc%nnest
      if ( summa1_struc(n)%nGRU > 0 ) then
         call summa_readRestart(summa1_struc(n), err, message)
         call handle_err(err, message)
      endif
   enddo
end subroutine summa2_readrst
