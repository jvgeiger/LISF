!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_metsim
! \label{readcrd_metsim}
!
! !REVISION HISTORY:
! 26 Jan 2007; Hiroko Kato, Initial Code adopted from readcrd_metsim.F90
! 25 Jun 2007; Hiroko Kato, upgraded to LISv5.0
! 15 May 2017: Bailing Li; Added changes for reading in version 2.2 data
! 22 Oct 2018: Daniel Sarmiento; Added changes to support version 3 data
! 29 Apr 2020: Zhuo Wang: Added MetSim forcing for LIS-SUMMA.
!
! !INTERFACE:    
subroutine readcrd_metsim()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_endrun, LIS_verify
  use metsim_forcingMod, only : metsim_struc
! !DESCRIPTION:
!
!  This routine reads the options specific to MetSim forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n, rc

  write(LIS_logunit,*)'[INFO] Using MetSim forcing'

  call ESMF_ConfigFindLabel(LIS_config,"MetSim forcing directory:",rc=rc)
  call LIS_verify(rc, 'MetSim forcing directory: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,metsim_struc(n)%metsimdir,rc=rc)
  enddo

  write(unit=LIS_logunit,fmt=*)'[INFO] Using MetSim forcing'
 
  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] MetSim forcing directory : ',metsim_struc(n)%metsimdir
 
     metsim_struc(n)%ncold = 464
     metsim_struc(n)%nrold = 224
     metsim_struc(n)%metsimtime1 = 3000.0
     metsim_struc(n)%metsimtime2 = 0.0
  enddo

end subroutine readcrd_metsim
