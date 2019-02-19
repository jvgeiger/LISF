!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: summa2_lsmMod
module summa2_lsmMod
! !USES:        
   use nrtype
   use summa2_module
   use summa_type
   use summa_util, only : stop_program
   use summa_init, only : summa_initialize
!
! !DESCRIPTION:
!  Module for 1-D land model driver variable initialization
!  
! \begin{description}
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[numout]
!    number of output times 
!  \item[outInterval]
!    output writing interval
!  \item[summa2open]
!    variable to keep track of opened files
!  \item[summa2]
!   summa2 LSM specific variables
! \end{description} 
!
   implicit none
  
   PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
   public :: summa2_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
   public :: summa1_struc
!EOP

   type(summa1_type_dec), allocatable :: summa1_struc(:)

   SAVE

contains

!BOP
! 
! !ROUTINE: summa2_lsm_ini
! \label{summa2_lsm_ini}
! 
! !INTERFACE:
subroutine summa2_lsm_ini()
! !USES:
   use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
   use LIS_coreMod,             only : LIS_rc
   use LIS_timeMgrMod,          only : LIS_clock, LIS_calendar, &
                                       LIS_update_timestep, LIS_registerAlarm
   use LIS_logMod,              only : LIS_verify
   
! !DESCRIPTION:        
!
!EOP
   implicit none
   integer             :: n
   integer             :: status
   integer(i4b)        :: err=0                      ! error code
   character(len=1024) :: message=''                 ! error message

   call summa2_readcrd()
   do n = 1, LIS_rc%nnest
      allocate(summa1_struc(n), stat=err)
      if (err/=0) then
         call stop_program(1, 'problem allocating master summa structure')
      endif
      call summa_initialize(summa1_struc(n), err, message)

      !allocate(summa1_struc(n)%summa2(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      summa1_struc(n)%numout = 0

      !call LIS_update_timestep(LIS_rc, n, summa1_struc(n)%ts)

      LIS_sfmodel_struc(n)%nsm_layers = 1
      LIS_sfmodel_struc(n)%nst_layers = 1
      allocate(LIS_sfmodel_struc(n)%lyrthk(1))
      LIS_sfmodel_struc(n)%lyrthk(1) = 1
      LIS_sfmodel_struc(n)%ts = summa1_struc(n)%ts
   enddo
end subroutine summa2_lsm_ini
end module summa2_lsmMod

