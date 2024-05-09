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
! !ROUTINE: dssat48_readrst
! \label{dssat48_readrst}
!
! !REVISION HISTORY: 
!  17 Apr 2024: Pang-Wei Liu
!
! !INTERFACE:
subroutine dssat48_readrst()
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc, LIS_surface
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
                             LIS_getNextUnitNumber,   &
                             LIS_releaseUnitNumber,   &
                             LIS_verify
  use dssat48_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!  This program reads restart files for DSSAT48.  This
!  includes all relevant water/energy storages and tile information. 
!  The following is the list of variables specified in the DSSAT48
!  restart file: 
!
!  \begin{verbatim}
!    nc, nr, ntiles        - grid and tile space dimensions 
!    dssat_sm_d            - Soil moisture of DSSAt (m3/m3)
!  \end{verbatim}
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!  reads a variable from the restart file
! \item[snowmodel\_coldstart](\ref{snowmodel_coldstart}) \newline
!   initializes the SnowModel state variables
! \end{description}
!EOP

  implicit none

  integer           :: t, l, col, row
  integer           :: nc, nr, npatch
  integer           :: n
  integer           :: ftn
  integer           :: status
  logical           :: file_exists
  character*20      :: wformat
! ___________________________________

   write(LIS_logunit,*) '[INFO] Call to the DSSAT48 Read restart routine ...'

   do n=1, LIS_rc%nnest
      wformat = dssat48_struc(n)%rformat

      ! Coldstart
      if (dssat48_struc(n)%dssatstartcode == "coldstart" ) then
      !if(LIS_rc%startcode .eq. "coldstart") then
         call dssat48_coldstart(LIS_rc%lsm_index)
      ! Restart
      !elseif(LIS_rc%startcode .eq. "restart") then
      elseif (dssat48_struc(n)%dssatstartcode == "restart" ) then
         ! Check if MicroMet option is set to "SnwoModel";
         !  If so, alert user that this option currently
         !  is not supported with reading in LIS restart files:
         !if( dssat48_struc(n)%sm_micromet_opt == "DssatModel" ) then
         !   write(LIS_logunit,*) "[ERR] Dssat48 restart file read is not currently"
         !   write(LIS_logunit,*) "  supported with 'Dssat48' forcing-format option."
         !   call LIS_endrun
         !endif
         ! check the existance of restart file
         inquire(file=dssat48_struc(n)%rfile, exist=file_exists)
         if(.not. file_exists) then
            write(LIS_logunit,*) "[ERR] Dssat48 restart file ", &
                                 dssat48_struc(n)%rfile," does not exist "
            write(LIS_logunit,*) " Program stopping ..."
            call LIS_endrun
         endif
         write(LIS_logunit,*) "[INFO] DSSAT48 restart file used: ",&
                               trim(dssat48_struc(n)%rfile)

         ! open restart file
         if(wformat .eq. "binary") then
            ftn = LIS_getNextUnitNumber()
            open(ftn, file=dssat48_struc(n)%rfile, form="unformatted")
            read(ftn) nc, nr, npatch  ! time, veg class, no. tiles

            ! check for grid space conflict
            if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
               write(LIS_logunit,*) "[ERR] "//trim(dssat48_struc(n)%rfile), &
                             ":: grid space mismatch - dssat48 run halted"
               call LIS_endrun
            endif

            if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
               write(LIS_logunit,*) "[ERR] dssat48 restart tile space mismatch"
               call LIS_endrun
            endif
         elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
            status = nf90_open(path=dssat48_struc(n)%rfile, &
                               mode=NF90_NOWRITE, ncid=ftn)
            call LIS_verify(status, "Error opening file "//dssat48_struc(n)%rfile)
#endif
         endif

         ! Soil moisture D1 (m3/m3)
           call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  dssat48_struc(n)%dssat48%SW(1), &
                                  varname="DSSATSMD1", &
                                  wformat=wformat)

         ! Soil moisture D2 (m3/m3)
           call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  dssat48_struc(n)%dssat48%SW(2), &
                                  varname="DSSATSMD2", &
                                  wformat=wformat)
         ! Soil moisture D3 (m3/m3)
           call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  dssat48_struc(n)%dssat48%SW(3), &
                                  varname="DSSATSMD3", &
                                  wformat=wformat)

         ! Soil moisture D4 (m3/m3)
           call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                                  dssat48_struc(n)%dssat48%SW(4), &
                                  varname="DSSATSMD4", &
                                  wformat=wformat)

         ! Close the restart file
         if(wformat .eq. "binary") then
            call LIS_releaseUnitNumber(ftn)
         elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
            status = nf90_close(ftn)
            call LIS_verify(status, "Error in nf90_close in dssat48_readrst")
#endif
         endif
!
!         ! Assign the read-in states to the SnowModel local states:
!         do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!            col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
!            row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
!
!            snow_d(col,row)        = snowmodel_struc(n)%sm(t)%snow_d 
!            snow_depth(col,row)    = snowmodel_struc(n)%sm(t)%snow_depth
!            canopy_int(col,row)    = snowmodel_struc(n)%sm(t)%canopy_int 
!            soft_snow_d(col,row)   = snowmodel_struc(n)%sm(t)%soft_snow_d 
!            ro_snow_grid(col,row)  = snowmodel_struc(n)%sm(t)%ro_snow_grid 
!            swe_depth(col,row)     = snowmodel_struc(n)%sm(t)%swe_depth
!            ro_soft_snow_old(col,row) = snowmodel_struc(n)%sm(t)%ro_soft_snow_old 
!            snow_d_init(col,row)   = snowmodel_struc(n)%sm(t)%snow_d_init 
!            swe_depth_old(col,row) = snowmodel_struc(n)%sm(t)%swe_depth_old 
!            canopy_int_old(col,row)= snowmodel_struc(n)%sm(t)%canopy_int_old 
!            topo(col,row)          = snowmodel_struc(n)%sm(t)%topo 
!            sum_sprec(col,row)     = snowmodel_struc(n)%sm(t)%sum_sprec 
!         enddo

      endif
   enddo

end subroutine dssat48_readrst

