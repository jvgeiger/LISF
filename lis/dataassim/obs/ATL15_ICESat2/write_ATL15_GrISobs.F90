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
! !ROUTINE: write_ATL15_GrISobs
! \label{write_ATL15_GrISobs}
! 
! !REVISION HISTORY: 
!  19 Dec 2023    Mahdi Navari;   Initial Specification
!
! !INTERFACE: 
subroutine write_ATL15_GrISobs(n, k, OBS_State)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_historyMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none

! !ARGUMENTS: 

  integer,     intent(in)  :: n 
  integer,     intent(in)  :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION: 
! 
! writes the transformed (interpolated/upscaled/reprojected)  
! ATL15_GrIS observations to a file
! 
!EOP
  type(ESMF_Field)         :: dhdtField
  logical                  :: data_update
  real, pointer            :: dhdtobs(:)
  real                     :: dhdtobs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01",dhdtField, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(dhdtField, localDE=0, farrayPtr=dhdtobs, rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call ATL15_GrIS_obsname(obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn,file=trim(obsname), form='unformatted')
     endif

     call LIS_writevar_gridded_obs(ftn,n,k,dhdtobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_ATL15_GrISobs

!BOP
! !ROUTINE: ATL15_GrIS_obsname
! \label{ATL15_GrIS_obsname}
! 
! !INTERFACE: 
subroutine ATL15_GrIS_obsname(obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)      :: obsname
! 
! !DESCRIPTION: 
! 
!EOP

  character(len=6) :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2)') &
       LIS_rc%yr, LIS_rc%mo

  obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//'/'//cdate1//   &
            '.1gs4r'
! TODO : Why for SM we add .a01.d01 to the filename.(ie., '/LISDAOBS_'//cdate1// &
!                                                           trim(cda)//trim(cdate)//'.1gs4r')  
  
end subroutine ATL15_GrIS_obsname
