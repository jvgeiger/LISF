!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! 
! !ROUTINE: dssat48_writerst
! \label{dssat48_writerst}
!
! !REVISION HISTORY:
!  26 Jun 2023: Pang-Wei Liu
!
!
! !INTERFACE:
subroutine dssat48_writerst(n)
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod,     only : LIS_logunit, LIS_getNextUnitNumber, &
                             LIS_releaseUnitNumber, LIS_verify
  use LIS_fileIOMod,  only : LIS_create_output_directory, &
                             LIS_create_restart_filename
  use dssat48_lsmMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for dssat48.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[dssat48model\_dump\_restart](\ref{dssat48model_dump_restart}) \newline
!   writes theDSSATModel variables into the restart file
! \end{description}
!EOP
  character(len=LIS_CONST_PATH_LEN) :: filen
  logical       :: alarmCheck
  integer       :: ftn
  integer       :: status
  character*20  :: wformat
  character*3   :: fnest

  write(fnest,'(i3.3)') n
  ! set restart alarm
  alarmCheck = LIS_isAlarmRinging(LIS_rc,"DSSAT48 restart alarm "//trim(fnest))

  write(LIS_logunit,*) '[INFO] Call to the DSSAT48 write restart routine ...'

  ! set restart file format (read from LIS configration file_
  wformat = trim(dssat48_struc(n)%rformat)

  if(alarmCheck .or. (LIS_rc%endtime ==1)) then
     if (LIS_masterproc) Then
        call LIS_create_output_directory("SURFACEMODEL")
        call LIS_create_restart_filename(n, filen, "SURFACEMODEL", &
                                         "DSSAT48",wformat=wformat)
        if(wformat .eq. "binary") then
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=filen,status="unknown", form="unformatted")
        elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
           status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
           call LIS_verify(status, &
                "Error in nf90_open in dssat48_writerst")
#endif
#if (defined USE_NETCDF3)
           status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
           call LIS_verify(status, &
                "Error in nf90_open in dssat48_writerst")
#endif
        endif
     endif

     ! Call routine to write out the model states to the file:
     call dssat48_dump_restart(n, ftn, wformat) !Here is where we put variables to write out

     if (LIS_masterproc) then
        if(wformat .eq. "binary") then
           call LIS_releaseUnitNumber(ftn)
        elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           status = nf90_close(ftn)
           call LIS_verify(status, &
                "Error in nf90_close in dssat48_writerst")
#endif
        endif
        write(LIS_logunit,*)&
             "[INFO] DSSAT48 archive restart written: ",trim(filen)
     endif
  endif
end subroutine dssat48_writerst

!BOP
!
! !ROUTINE: dssat48_dump_restart
! \label{dssat48_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!  17 Apr 2024: Pang-Wei Liu
!
! !INTERFACE:
subroutine dssat48_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use dssat48_lsmMod

    implicit none

    integer, intent(in) :: ftn
    integer, intent(in) :: n
    character(len=*), intent(in) :: wformat
!
! !DESCRIPTION:
!  This routine gathers the necessary restart variables and performs
!  the actual write statements to create the restart files.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    unit number for the restart file
!   \item[wformat]
!    restart file format (binary/netcdf)
!  \end{description}
!
!  The following is the list of variables written in the SnowModel
!  restart file:
!
!  \begin{verbatim}
!    nc, nr, ntiles        - grid and tile space dimensions
!    dssat_sm_d            - Soil moisture of DSSAt (m3/m3)
!  \end{verbatim}
!
! The routines invoked are:
! \begin{description}
!   \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart})\\
!      writes the global header information
!   \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart})\\
!      writes the header information for a variable
!   \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart})\\
!      close the header
!   \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart})\\
!      writes a variable to the restart file
! \end{description}
! 
!EOP 

    integer :: l, t
    real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index))
    integer :: dimID(11)

    integer :: SMD1_d_ID, SMD2_d_ID, SMD3_d_ID, SMD4_d_ID  ! 1


!- Write the header of the restart file
   call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
                                       "DSSAT48", &
                                       dim1=1, &   ! Set as one for now
                                       dimID=dimID, &
                                       output_format = trim(wformat))

   !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max

!-  DSSAT48 states, based on preprocess_code
   ! write header for DSSAT Soil Moisture at 1st layer
   call LIS_writeHeader_restart(ftn, n, dimID, SMD1_d_ID, "DSSATSMD1", &
                                "DSSAT Soil Moisture at 1st Layer", &
                                "m3/m3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write header for DSSAT Soil Moisture at 2nd layer
   call LIS_writeHeader_restart(ftn, n, dimID, SMD2_d_ID, "DSSATSMD2", &
                                "DSSAT Soil Moisture at 2nd Layer", &
                                "m3/m3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write header for DSSAT Soil Moisture at 3rd layer
   call LIS_writeHeader_restart(ftn, n, dimID, SMD3_d_ID, "DSSATSMD3", &
                                "DSSAT Soil Moisture at 3rd Layer", &
                                "m3/m3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

   ! write header for DSSAT Soil Moisture at 4th layer
   call LIS_writeHeader_restart(ftn, n, dimID, SMD4_d_ID, "DSSATSMD4", &
                                "DSSAT Soil Moisture at 4th Layer", &
                                "m3/m3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
   ! close header of restart file
   call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, dssat48_struc(n)%rstInterval)


!- Write state variables into restart file:

   ! Soil moisture D1 (m3/m3)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             dssat48_struc(n)%dssat48%SW(1), &
                             varid=SMD1_d_ID, dim=1, wformat=wformat)

   ! Soil moisture D2 (m3/m3)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             dssat48_struc(n)%dssat48%SW(2), &
                             varid=SMD2_d_ID, dim=1, wformat=wformat)

   ! Soil moisture D3 (m3/m3)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             dssat48_struc(n)%dssat48%SW(3), &
                             varid=SMD3_d_ID, dim=1, wformat=wformat)

   ! Soil moisture D4 (m3/m3)
   call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, &
                             dssat48_struc(n)%dssat48%SW(4), &
                             varid=SMD4_d_ID, dim=1, wformat=wformat)

end subroutine dssat48_dump_restart

