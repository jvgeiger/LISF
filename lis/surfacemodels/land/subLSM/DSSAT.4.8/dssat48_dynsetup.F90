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
! !ROUTINE: dssat48_dynsetup
! \label{dssat48_dynsetup}
!
! !REVISION HISTORY:
!  26 Jun 2023: Pang-Wei Liu
!
! !INTERFACE:
subroutine dssat48_dynsetup(n)
! !USES: 
  use LIS_coreMod, only : LIS_rc,LIS_domain, LIS_surface
  use LIS_logMod,  only : LIS_logunit, LIS_endrun

  use dssat48_lsmMod, only : dssat48_struc
  use LIS_timeMgrMod, only : LIS_date2time,LIS_tick
! 
! !DESCRIPTION: 
!  This routine sets up the time-dependent variables in DSSAT48.
! 
!EOP   

  implicit none
  integer, intent(in) :: n

  write(LIS_logunit,*) '[INFO] Call to the DSSAT48 dynamic setup routine ...'

end subroutine dssat48_dynsetup
