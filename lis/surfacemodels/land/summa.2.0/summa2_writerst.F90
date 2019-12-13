!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: summa2_writerst
! \label{summa2_writerst}
! 
! !REVISION HISTORY:
!   June 2018: James Geiger
! 20 Nov 2019: Zhuo Wang, added support for summa2
! 
! !INTERFACE:
subroutine summa2_writerst(n)
! !USES:
     use LIS_coreMod,    only : LIS_rc, LIS_masterproc
     use LIS_timeMgrMod, only : LIS_isAlarmRinging
     use LIS_logMod,     only : LIS_logunit, LIS_getNextUnitNumber, &
                                LIS_releaseUnitNumber , LIS_verify
     use LIS_fileIOMod,  only : LIS_create_output_directory, &
                                LIS_create_restart_filename
     use summa2_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!  This program writes restart files for summa2 model.
!  This includes all relevant water/energy storage and tile information.
!
!  Forcing-only (summa2) option for calling the write restart routines. 
! 
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory})\\
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename})\\
!  generates a timestamped restart filename
! \item[summa2\_dump\_restart](\ref{summa2_dump_restart})\\
!   writes the summa2 variables into the restart file
! \end{description}
!EOP

    character*100 :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status

    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "summa2 restart alarm")

    ! set restart file format (read from LIS configration file)
    wformat = trim(summa1_struc(n)%rformat)

    if(alarmCheck .or. (LIS_rc%endtime == 1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", &
                                            "SUMMA2",wformat=wformat)
            if(wformat .eq. "binary") then
                 ftn = LIS_getNextUnitNumber()
                 open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                 status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                 call LIS_verify(status, &
                      "Error in nf90_open in summa2_writerst")
#endif
#if (defined USE_NETCDF3)
                 status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                 call LIS_verify(status, &
                      "Error in nf90_open in summa2_writerst")
#endif
             endif
        endif

        call summa2_dump_restart(n, ftn, wformat)

        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, &
                     "Error in nf90_close in summa2_writerst")
#endif
             endif
             write(*,*)&
                  "[INFO] SUMMA.2.0 archive restart written: ",trim(filen)
         endif
    endif

end subroutine summa2_writerst

!BOP
!
! !ROUTINE: summa2_dump_restart
! \label{summa2_dump_restart}
!
! !REVISION HISTORY:
! 20 Nov 2019: Zhuo Wang, added support for summa2
! !INTERFACE:
subroutine summa2_dump_restart(n, ftn, wformat)

! !USES:
   use LIS_coreMod,    only : LIS_rc, LIS_masterproc
   use LIS_logMod,     only : LIS_logunit
   use LIS_historyMod
   use summa2_lsmMod
   USE nrtype

   use globalData,     only : gru_struc  ! gru-hru mapping structures
   use var_lookup,     only : iLookPROG  ! named variables for local state variables
   use var_lookup,     only : iLookINDEX

   implicit none
  
   integer, intent(in)                 :: ftn
   integer, intent(in)                 :: n
   character(len=*), intent(in)        :: wformat

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
!  The following is the list of variables written in the summa2 restart file:
!  \begin{verbatim}
!   nc, nr, ntiles     - grid and tile space dimensions
!   dt_init            - length of initial time step at start of next data interval
!   scalarCanopyIce    - mass of ice on the vegetation canopy
!   scalarCanopyLiq    - mass of liquid water on the vegetation canopy
!   scalarCanopyWat    - mass of total water on the vegetation canopy
!   scalarCanairTemp   - temperature of the canopy air space
!   scalarCanopyTemp   - temperature of the vegetation canopy
!   spectralSnowAlbedoDiffuse - diffuse snow albedo for individual spectral bands
!   scalarSnowAlbedo   - snow albedo for the entire spectral band
!   scalarSnowDepth    - total snow depth
!   scalarSWE          - snow water equivalent
!   scalarSfcMeltPond  - ponded water caused by melt of the snow without a layer
!   mLayerTemp         - temperature of each layer
!   mLayerVolFracIce   - volumetric fraction of ice in each layer
!   mLayerVolFracLiq   - volumetric fraction of liquid water in each layer
!   mLayerVolFracWat   - volumetric fraction of total water in each layer
!   mLayerMatricHead   - matric head of water in the soil
!   scalarAquiferStorage - water required to bring aquifer to the bottom of the soil profile
!   scalarSurfaceTemp  - surface temperature (just a copy of the upper-layer temperature)
!   mLayerDepth        - depth of each layer
!   mLayerHeight       - height of the layer mid-point (top of soil = 0)
!   iLayerHeight       - height of the layer interface (top of soil = 0)
!   nSnow              - number of snow layers
!   nSoil              - number of soil layers
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
    integer :: dimID(10)
    real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index))
    real    :: tmptilen_1d(LIS_rc%npatch(n, LIS_rc%lsm_index))
    integer :: nSnowData(LIS_rc%npatch(n, LIS_rc%lsm_index)) 
    integer :: nSoilData(LIS_rc%npatch(n, LIS_rc%lsm_index))

    ! state variables for vegetation
    integer :: dt_init_ID             ! length of initial time step at start of next data interval (s)
    integer :: CanopyIce_ID           ! mass of ice on the vegetation canopy (kg m-2)
    integer :: CanopyLiq_ID           ! mass of liquid water on the vegetation canopy (kg m-2)
    integer :: CanopyWat_ID           ! mass of total water on the vegetation canopy (kg m-2)
    integer :: CanairTemp_ID          ! temperature of the canopy air space (Pa)
    integer :: CanopyTemp_ID          ! temperature of the vegetation canopy (K)
    ! state variables for snow
    integer :: spectralSnowAlbedoDiffuse_ID ! diffuse snow albedo for individual spectral bands (-)
    integer :: SnowAlbedo_ID          ! snow albedo for the entire spectral band (-)
    integer :: SnowDepth_ID           ! total snow depth (m)
    integer :: SWE_ID                 ! snow water equivalent (kg m-2)
    integer :: SfcMeltPond_ID         ! ponded water caused by melt of the "snow without a layer" (kg m-2)
    ! state variables for the snow+soil domain
    integer :: mLayerTemp_ID                ! temperature of each layer (K)
    integer :: mLayerVolFracIce_ID          ! volumetric fraction of ice in each layer (-)
    integer :: mLayerVolFracLiq_ID          ! volumetric fraction of liquid water in each layer (-)
    integer :: mLayerVolFracWat_ID          ! volumetric fraction of total water in each layer (-)
    integer :: mLayerMatricHead_ID          ! matric head of water in the soil (m)
    ! other state variables
    integer :: AquiferStorage_ID      ! relative aquifer storage -- above bottom of the soil profile (m)
    integer :: SurfaceTemp_ID         ! surface temperature (K)
    ! coordinate variables
    integer :: mLayerDepth_ID               ! depth of each layer (m)
    integer :: mLayerHeight_ID              ! height at the mid-point of each layer (m)
    integer :: iLayerHeight_ID              ! height of the layer interface; top of soil = 0 (m)

    integer :: nSnow_ID                     ! number of snow layers
    integer :: nSoil_ID                     ! number of soil layers

    integer :: nSoil_max               ! number of soil layers
    integer :: nSnow_max               ! number of snow layers

    integer :: nlayer_real             ! rean snow number + nsoil = real snow number + 8

    integer :: spectral
    integer :: midSoil
    integer :: midToto
    integer :: ifcSoil
    integer :: ifcToto
    integer :: midSnow
    integer :: ifcSnow
!--------------------------------------------------------------------------------------------------
    ! write the header of the restart file

    nSnowData = summa1_struc(n)%indxStruct%gru(:)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
    nSoilData = summa1_struc(n)%indxStruct%gru(:)%hru(1)%var(iLookINDEX%nSoil)%dat(1)

    nSnow_max = 0 
    nSoil_max = 0

    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
      if(nSnowData(t) .gt. nSnow_max) then
         nSnow_max = nSnowData(t)
      endif

      if(nSoilData(t) .gt. nSoil_max) then
         nSoil_max = nSoilData(t)
      endif 
    enddo

    if(nSnow_max .gt. 0) then
      nSnow_max=5
    endif

    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, "SUMMA2",              &
                           dim1=2, dim2=8, dim3=13, dim4=9, dim5=14, dim6=5, dim7=6, dimID=dimID, output_format="netcdf")

    call LIS_writeHeader_restart(ftn, n, dimID, dt_init_ID, "dt_init", &
                                 "length of initial time step at start of next data interval", &
                                 "s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    call LIS_writeHeader_restart(ftn, n, dimID, CanopyIce_ID, "scalarCanopyIce", &
                                 "mass of ice on the vegetation canopy (instant)", &
                                 "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    call LIS_writeHeader_restart(ftn, n, dimID, CanopyLiq_ID, "scalarCanopyLiq", &
                                 "mass of liquid water on the vegetation canopy (instant)", &
                                 "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    call LIS_writeHeader_restart(ftn, n, dimID, CanopyWat_ID, "scalarCanopyWat", &
                                 "mass of total water on the vegetation canopy (instant)", &
                                 "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    call LIS_writeHeader_restart(ftn, n, dimID, CanairTemp_ID, "scalarCanairTemp", &
                                  "temperature of the canopy air space", &
                                  "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    call LIS_writeHeader_restart(ftn, n, dimID, CanopyTemp_ID, "scalarCanopyTemp", &
                                  "temperature of the vegetation canopy", &
                                  "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    call LIS_writeHeader_restart(ftn, n, dimID, spectralSnowAlbedoDiffuse_ID, "spectralSnowAlbedoDiffuse", &
                                  "diffuse snow albedo for individual spectral bands", &
                                  "-", vlevels=2, valid_min=0.0, valid_max=1.0, var_flag = "dim1")

    call LIS_writeHeader_restart(ftn, n, dimID, SnowAlbedo_ID, "scalarSnowAlbedo", &
                                  "snow albedo for the entire spectral band", &
                                  "-", vlevels=1, valid_min=0.0, valid_max=1.0)

    call LIS_writeHeader_restart(ftn, n, dimID, SnowDepth_ID, "scalarSnowDepth", &
                                  "total snow depth", &
                                  "m", vlevels=1, valid_min=0.0, valid_max=1.0)

    call LIS_writeHeader_restart(ftn, n, dimID, SWE_ID, "scalarSWE", &
                                  "snow water equivalent", &
                                  "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    call LIS_writeHeader_restart(ftn, n, dimID, SfcMeltPond_ID, "scalarSfcMeltPond", &
                                  "ponded water caused by melt of the snow without a layer", &
                                  "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    call LIS_writeHeader_restart(ftn, n, dimID, mLayerMatricHead_ID, "mLayerMatricHead", &
                                 "matric head of water in the soil", &
                                 "m", vlevels=8, &
                                 valid_min=-99999.0, valid_max=99999.0, var_flag = "dim2")
     
    call LIS_writeHeader_restart(ftn, n, dimID, AquiferStorage_ID, "scalarAquiferStorage", &
                                 "relative aquifer storage -- above bottom of the soil profile", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
     
    call LIS_writeHeader_restart(ftn, n, dimID, SurfaceTemp_ID, "scalarSurfaceTemp", &
                                 "surface temperature", &
                                 "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable with multiple layers
      call LIS_writeHeader_restart(ftn, n, dimID, mLayerTemp_ID, "mLayerTemp", &
                                   "temperature of each layer", &
                                   "K", vlevels=13, &
                                   valid_min=-99999.0, valid_max=99999.0, var_flag = "dim3")
     
       call LIS_writeHeader_restart(ftn, n, dimID, mLayerVolFracIce_ID, "mLayerVolFracIce", &
                                    "volumetric fraction of ice in each layer", &
                                    "-", vlevels=13, &
                                    valid_min=-99999.0, valid_max=99999.0, var_flag = "dim3")
     
       call LIS_writeHeader_restart(ftn, n, dimID, mLayerVolFracLiq_ID, "mLayerVolFracLiq", &
                                    "volumetric fraction of liquid water in each layer", &
                                    "-", vlevels=13, &
                                    valid_min=-99999.0, valid_max=99999.0, var_flag = "dim3")
     
       call LIS_writeHeader_restart(ftn, n, dimID, mLayerVolFracWat_ID, "mLayerVolFracWat", &
                                    "volumetric fraction of total water in each layer", &
                                    "-", vlevels=13, &
                                    valid_min=-99999.0, valid_max=99999.0, var_flag = "dim3")
     
       call LIS_writeHeader_restart(ftn, n, dimID, mLayerDepth_ID, "mLayerDepth", &
                                    "depth of each layer", &
                                    "m", vlevels=13, &
                                    valid_min=-99999.0, valid_max=99999.0, var_flag = "dim3")
     
       call LIS_writeHeader_restart(ftn, n, dimID, mLayerHeight_ID, "mLayerHeight", &
                                    "height of the layer mid-point (top of soil = 0)", &
                                    "m", vlevels=13, &
                                    valid_min=-99999.0, valid_max=99999.0, var_flag = "dim3")
     
       call LIS_writeHeader_restart(ftn, n, dimID, iLayerHeight_ID, "iLayerHeight", &
                                    "height of the layer interface (top of soil = 0)", &
                                    "m", vlevels=14, &
                                    valid_min=-99999.0, valid_max=99999.0, var_flag = "dim5")

   call LIS_writeHeader_restart(ftn, n, dimID, nSNow_ID, "nSnow", &
                                  "number of snow layers", &
                                  "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0) 

   call LIS_writeHeader_restart(ftn, n, dimID, nSoil_ID, "nSoil", &
                                  "number of soil layers", &
                                  "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0) 

    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, summa1_struc(n)%rstInterval)


    ! write state variables into restart file
    tmptilen_1d = 0

    tmptilen_1d =  summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%dt_init)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=dt_init_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarCanopyIce)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=CanopyIce_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarCanopyLiq)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=CanopyLiq_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarCanopyWat)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=CanopyWat_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarCanairTemp)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=CanairTemp_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarCanopyTemp)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=CanopyTemp_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarSnowAlbedo)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=SnowAlbedo_ID, dim=1, wformat=wformat)
        
    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarSnowDepth)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=SnowDepth_ID, dim=1, wformat=wformat)
        
    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarSWE)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=SWE_ID, dim=1, wformat=wformat)
        
    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarSfcMeltPond)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=SfcMeltPond_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarAquiferStorage)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=AquiferStorage_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = summa1_struc(n)%progStruct%gru(:)%hru(1)%var(iLookPROG%scalarSurfaceTemp)%dat(1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=SurfaceTemp_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = nSnowData
!   write(*,*)'nSnow=',tmptilen_1d(1:5)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=nSnow_ID, dim=1, wformat=wformat)

    tmptilen_1d = 0
    tmptilen_1d = nSoilData
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                              varid=nSoil_ID, dim=1, wformat=wformat)

    ! diffuse snow albedo for individual spectral bands
    do l=1, 2 ! spectral
      tmptilen = 0
      do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
       tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(l)
      enddo
      call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                 varid=spectralSnowAlbedoDiffuse_ID, dim=l, wformat=wformat)
    enddo
!----------------------------------------------
    ! temperature of each layer
      do l=1, 13 
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         nlayer_real = nSnowData(t) + nSoilData(t)
         if(l .le. nlayer_real) then
          tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%mLayerTemp)%dat(l)
         endif
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varid=mLayerTemp_ID, dim=l, wformat=wformat)
      enddo
!----------------------------------------------
    ! volumetric fraction of ice in each layer
       do l=1, 13 
         tmptilen = 0
         do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          nlayer_real = nSnowData(t) + nSoilData(t)
          if(l .le. nlayer_real) then
           tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat(l)
          endif
         enddo
         call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varid=mLayerVolFracIce_ID, dim=l, wformat=wformat)
        enddo
!----------------------------------------------
    ! volumetric fraction of liquid water in each layer
       do l=1, 13 ! midToto  (In original coldState.nc)
         tmptilen = 0
         do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
           nlayer_real = nSnowData(t) + nSoilData(t)
           if(l .le. nlayer_real) then
             tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat(l)
           endif
         enddo
         call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                   varid=mLayerVolFracLiq_ID, dim=l, wformat=wformat)
       enddo
!----------------------------------------------
    ! volumetric fraction of total water in each layer
      do l=1, 13 ! midToto  (In original coldState.nc)
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          nlayer_real = nSnowData(t) + nSoilData(t)
          if(l .le. nlayer_real) then
           tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat(l)
          endif
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=mLayerVolFracWat_ID, dim=l, wformat=wformat)
       enddo
!----------------------------------------------
    ! matric head of water in the soil
    do l=1, 8 ! midSoil
      tmptilen = 0
      do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
       tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat(l)
      enddo
      call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                 varid=mLayerMatricHead_ID, dim=l, wformat=wformat)
    enddo
!----------------------------------------------
    ! depth of each layer
       do l=1, 13
         tmptilen = 0
         do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          nlayer_real = nSnowData(t) + nSoilData(t)
          if(l .le. nlayer_real) then
           tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%mLayerDepth)%dat(l)
          endif
         enddo
         call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                   varid=mLayerDepth_ID, dim=l, wformat=wformat)
       enddo
!----------------------------------------------
    ! height at the mid-point of each layer
       do l=1, 13
         tmptilen = 0
         do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          nlayer_real = nSnowData(t) + nSoilData(t)
          if(l .le. nlayer_real) then
            tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%mLayerHeight)%dat(l)
          endif
         enddo
         call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                   varid=mLayerHeight_ID, dim=l, wformat=wformat)
        enddo
!     endif
!----------------------------------------------
    ! height of the layer interface; top of soil = 0 
      do l=1, 14  ! ifcToto (In original coldState.nc)
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         nlayer_real = nSnowData(t) + nSoilData(t)
         if(l .le. nlayer_real) then
           tmptilen(t) = summa1_struc(n)%progStruct%gru(t)%hru(1)%var(iLookPROG%iLayerHeight)%dat(l-1)
         endif
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=iLayerHeight_ID, dim=l, wformat=wformat)
      enddo
!----------------------------------------------

end subroutine summa2_dump_restart 
