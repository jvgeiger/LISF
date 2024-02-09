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
!
! !ROUTINE: readcrd_MAR
!  \label{readcrd_MAR}
!
! !REVISION HISTORY:
! 15 Aug 2023: Mahdi Navari; Initial Specification
!
! !INTERFACE:    
subroutine readcrd_MAR()
! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_config
  use LIS_logMod,        only : LIS_logunit, LIS_verify
  use MAR_forcingMod, only : MAR_struc

  implicit none

! !DESCRIPTION:
!
!  This routine reads the options specific to MAR output forcing from 
!  the LIS configuration file. 
!  
!EOP

  integer :: n,rc
  
  call ESMF_ConfigFindLabel(LIS_config,"MAR forcing directory:",rc=rc)
  call LIS_verify(rc, 'MAR forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,MAR_struc(n)%MARdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     MAR_struc(n)%nest_id = 1
  enddo

  write(unit=LIS_logunit,fmt=*)'[INFO] Using MAR forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] MAR forcing directory :',trim(MAR_struc(n)%MARdir)

     MAR_struc(n)%MARtime1 = 3000.0
     MAR_struc(n)%MARtime2 = 0.0
  enddo

end subroutine readcrd_MAR
