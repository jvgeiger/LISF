!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: summa2_writerst
! \label{summa2_writerst}
! 
! !INTERFACE:
subroutine summa2_writerst(n)
! !USES:

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!
!  Forcing-only (summa2) option for calling the write restart routines. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

end subroutine summa2_writerst
