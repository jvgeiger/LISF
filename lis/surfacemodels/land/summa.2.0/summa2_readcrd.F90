!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: summa2_readcrd
! \label{summa2_readcrd}
!
! !REVISION HISTORY:
!
! !INTERFACE:    
subroutine summa2_readcrd()

! !USES:
   use ESMF
   use LIS_coreMod,     only : LIS_rc, LIS_config
   use LIS_timeMgrMod, only : LIS_parseTimeString
   use LIS_logMod,      only : LIS_logunit, LIS_verify
   use summa2_lsmMod, only : summa1_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to summa2 LSM 
!  option from the LIS configuration file. 
!  
!EOP
   implicit none

   integer :: rc
   integer :: n
   character*10  :: time
   real    :: ts_tmp

   call ESMF_ConfigFindLabel(LIS_config,"SUMMA.2.0 model timestep:",rc=rc)
   do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
      call LIS_verify(rc,'SUMMA.2.0 model timestep: not defined')

      call LIS_parseTimeString(time,ts_tmp)
      summa1_struc(n)%ts = ts_tmp
   enddo

   call ESMF_ConfigFindLabel(LIS_config,"SUMMA.2.0 master file:",rc=rc)
   do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,&
         summa1_struc(n)%summaFileManagerFile,rc=rc)
      call LIS_verify(rc,'SUMMA.2.0 master file: not defined')
   enddo

   write(LIS_logunit,*)'[INFO] Running SUMMA 2.0 LSM Option:'
   !do n=1,LIS_rc%nnest
   !   summa1_struc(n)%summa2open=0
   !enddo

end subroutine summa2_readcrd
