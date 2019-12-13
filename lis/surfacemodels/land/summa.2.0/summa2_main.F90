!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: summa2_main
! \label{summa2_main}
!
! !INTERFACE:
subroutine summa2_main(n)
! !USES:
   use LIS_coreMod,    only : LIS_rc
   use nrtype,         only : i4b
   use summa2_lsmMod
   use summa_util,     only : handle_err

   use summa_modelRun, only : summa_runPhysics

   use LIS_timeMgrMod, only : LIS_date2time, LIS_tick

  implicit none
! !ARGUMENTS: 
   integer, intent(in) :: n
!
! !DESCRIPTION:
! 
!  This is the entry point for calling the SUMMA 2.0 LSM physics.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}

!EOP

   integer(i4b)        :: err=0      ! error code
   character(len=1024) :: message='' ! error message

   ! Use LIS_rc%tscount(n) in place of modelTimeStep.
   ! modelTimeStep is an integer(i4b), 1 through end of run.
   ! It is used only to determine whether SUMMA is at the start of the run
   ! (modelTimeStep == 1).
   ! ! run the summa physics for one time step
   if ( summa1_struc(n)%nGRU > 0 ) then
      call summa_runPhysics(LIS_rc%tscount(n), summa1_struc(n), err, message)
      call handle_err(err, message)
   endif

   call summa2_output(n, summa1_struc(n))

end subroutine summa2_main
