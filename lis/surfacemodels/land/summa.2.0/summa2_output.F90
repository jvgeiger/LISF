!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: summa2_output
! \label{summa2_output}

! !REVISION HISTORY:
!   Jun 2019: James Geiger
!   Aug 2019: Zhuo Wang, added support for summa2
!   Aug 2020: David Mocko, cleaned up routine and added output variables
!
! !INTERFACE:
 subroutine summa2_output(n)
! !USES:
   use LIS_coreMod,    only : LIS_rc, LIS_surface
   use lis_historymod
   use LIS_histDataMod

   use summa2_lsmMod
   use nrtype,         only : i4b
   use globalData,     only : gru_struc  ! gru-hru mapping structures

   ! named variables
   use var_lookup,     only : iLookFLUX  ! named variables for model flux variables
   use var_lookup,     only : iLookFORCE ! named variables for forcing variables
   use var_lookup,     only : iLookPROG  ! named variables for prognostic state variables
   use var_lookup,     only : iLookDIAG  ! named variables for diagnostic variables
   use var_lookup,     only : iLookPARAM ! named variables for parameter variables
   use var_lookup,     only : iLookATTR  ! named variables for local attribute variables
   use var_lookup,     only : iLookTYPE  ! named variables for type classification variables
   use var_lookup,     only : iLookINDEX ! named variables for local column index variables
   use var_lookup,     only : iLookBVAR  ! named variables for basin-average model variables
   use var_lookup,     only : iLookID    ! named variables for local column model parameters

   implicit none
! !ARGUMENTS:
   integer, intent(in)              :: n
!
! !DESCRIPTION:
!
!  This is the entry point for processing SUMMA 2.0 LSM output.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

   real                             :: dat, dat2
   integer                          :: l, t
   integer                          :: nlayer_real      ! rean snow number + nsoil
   integer(i4b)                     :: iGRU, iHRU, nGRU
   integer(i4b)                     :: nHRUrun
   integer(i4b),allocatable         :: nSnowData(:)
   integer(i4b),allocatable         :: nSoilData(:)
! ---------------------------------------------------------------------------------------

   nGRU = LIS_rc%ntiles(n)
   nHRUrun = sum(gru_struc%hruCount)

   allocate(nSnowData(nHRUrun))
   allocate(nSoilData(nHRUrun))

   do iGRU = 1,nGRU
      do iHRU = 1,gru_struc(iGRU)%hruCount
         nSnowData(iGRU) = summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
         nSoilData(iGRU) = summa1_struc(n)%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)
      enddo
   enddo

   do iGRU = 1,nGRU
      do iHRU = 1,gru_struc(iGRU)%hruCount

!        if (gru_struc(iGRU)%hruCount /= 1) then
!           print *,'GREP: iGRU, hruCount', iGRU, gru_struc(iGRU)%hruCount
!        endif

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1) +       &
               summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SWNET,       &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarLWNetGround)%dat(1) +               &
               summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarLWNetCanopy)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_LWNET,       &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarCanopyIce)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOPICE,    &
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarCanopyLiq)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOPLIQ,    &
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarCanopyWat)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOPWAT,    &
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarCanopyWat)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOPINT,    &
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarSfcMeltPond)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SFCMELTPOND, &
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarSWE)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SWE,         &
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarSnowDepth)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SNOWDEPTH,   &
               value=dat,                                              &
               vlevel=1, unit="m ", direction="-",                     &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookDIAG%scalarGroundSnowFraction)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SNOWCOVER,   &
               value=dat,                                              &
               vlevel=1, unit="-", direction="-",                      &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarAquiferStorage)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_AQUIFERSTORAGE, &
               value=dat * 1000.,                                      &
               vlevel=1, unit="mm", direction="-",                     &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookDIAG%scalarSoilCompress)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SOILCOMPRESS,&
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookDIAG%scalarTotalSoilLiq)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SOILLIQ,     &
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookDIAG%scalarTotalSoilIce)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SOILICE,     &
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

! Note: LIS_MOC_TOTALSOILWAT is the sum of multi-layer soil moisture
!       LIS_MOC_SOILMOIST is the multi-layer soil moisture
         dat = summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookDIAG%scalarTotalSoilWat)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_TOTALSOILWAT,&
               value=dat,                                              &
               vlevel=1, unit="kg m-2", direction="-",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOPAIRNETFLUX, &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOPNETFLUX,&
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_GROUNDNETFLUX, &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SAV,         &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="IN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SAG,         &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="IN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarLWNetCanopy)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_IRC,         &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarLWNetGround)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_IRB,         &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSenHeatTotal)%dat(1) * -1.0
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_QH,          &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="UP",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSenHeatCanopy)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SHC,         &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="UP",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSenHeatGround)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SHB,         &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="UP",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarLatHeatTotal)%dat(1) * -1.0
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_QLE,         &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="UP",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_FCEV,        &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="UP",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_FCTR,        &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="UP",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarLatHeatGround)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_EVG,         &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="UP",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOPADVECT, &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_GROUNDADVECT,&
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopySublimation)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SUBCANOP,    &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSnowSublimation)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SUBSNOW,     &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopyTranspiration)%dat(1) * -1.0
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_TVEG,        &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="UP",            &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopyEvaporation)%dat(1) * -1.0
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_ECANOP,      &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="UP",            &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarGroundEvaporation)%dat(1) * -1.0
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_ESOIL,       &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="UP",            &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarThroughfallSnow)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_QST,         &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarThroughfallRain)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_THROUGHRAIN, &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOPYSNOWUNLOADING, &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOP_LIQDRAINAGE, &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarCanopyMeltFreeze)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_CANOP_MELT_FREEZE, &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarRainPlusMelt)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_RAIN_MELT,   &
               value=dat,                                              &
               vlevel=1, unit="m s-1", direction="-",                  &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSurfaceRunoff)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_QS,          &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarTotalRunoff)%dat(1) -               &
               summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSurfaceRunoff)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_QSB,         &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSoilBaseflow)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SOIL_BASEFLOW, &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSoilDrainage)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SOIL_DRAINAGE, &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarAquiferRecharge)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_AQUIFER_RECHARGE, &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarAquiferTranspire)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_AQUIFER_TRANSPIRE, &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarAquiferBaseflow)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_AQUIFER_BASEFLOW, &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="-",             &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarTotalET)%dat(1) * -1.0
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_EVAP,        &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="UP",            &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarTotalRunoff)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_QSTOT,       &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%bvarStruct%gru(iGRU)%                   &
               var(iLookBVAR%averageRoutedRunoff)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_ROUTEDRUNOFF,&
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSnowDrainage)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SNOWDRAINAGE,&
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarInfiltration)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_INFILTRATION,&
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="IN",            &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarNetRadiation)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_RNET,        &
               value=dat,                                              &
               vlevel=1, unit="W m-2", direction="DN",                 &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarRainfall)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_RAINF,       &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="DN",            &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%fluxStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookFLUX%scalarSnowfall)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SNOWF,       &
               value=dat,                                              &
               vlevel=1, unit="kg m-2 s-1", direction="DN",            &
               surface_type=LIS_rc%lsm_index)

!         dat = summa1_struc(n)%typeStruct%gru(iGRU)%hru(iHRU)%         &
!               var(iLookTYPE%hruId)
         dat = summa1_struc(n)%idStruct%gru(iGRU)%hru(iHRU)%           &
               var(iLookID%hruId)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_HRUID,       &
               value=dat,                                              &
               vlevel=1, unit="-", direction="-",                      &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%attrStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookATTR%elevation)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_ELEVATION,   &
               value=dat,                                              &
               vlevel=1, unit="m", direction="-",                      &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%mparStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPARAM%albedoMax)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_MXSNALBEDO,  &
               value=dat,                                              &
               vlevel=1, unit="-", direction="-",                      &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%diagStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookDIAG%scalarGreenVegFraction)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_GREENNESS,   &
               value=dat,                                              &
               vlevel=1, unit="-", direction="-",                      &
               surface_type=LIS_rc%lsm_index)

         dat = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%         &
               var(iLookPROG%scalarSurfaceTemp)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_AVGSURFT,    &
               value=dat,                                              &
               vlevel=1, unit="K", direction="-",                      &
               surface_type=LIS_rc%lsm_index)

! Output multi-layer variables
         nlayer_real = nSoilData(iGRU) + nSnowData(iGRU)

! Temperature of each layer (nsoil+nsnow) (K)
!         do l = 1,13
         do l = 1,8
            dat2 = 0
            if (l.le.nlayer_real) then
               dat2 = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%  &
                      var(iLookPROG%mLayerTemp)%dat(l)
            endif   ! l < nlayer_real
            call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_MLAYERTEMP,&
                  value = dat2,                                        &
                  vlevel=l, unit="K", direction="-",                   &
                  surface_type = LIS_rc%lsm_index)
         enddo ! l, layers

! Matric head of water in the soil
!         do l = 1,8
         do l = 1,3
            dat2 = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%     &
                   var(iLookPROG%mLayerMatricHead)%dat(l)
            call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_MLAYERMOIST,&
                  value = dat2,                                        &
                  vlevel=l, unit="m", direction="-",                   &
                  surface_type = LIS_rc%lsm_index)

            dat2 = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%     &
                   var(iLookPROG%mLayerVolFracWat)%dat(l)
            call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_SOILMOIST,&
                  value = dat2,                                        &
                  vlevel=l, unit="m^3 m-3", direction="-",             &
                  surface_type = LIS_rc%lsm_index)
         enddo ! l, layers

! Volumetric fraction of total water in each layer of nsoil+nsnow
!         do l = 1,13
         do l = 1,8
            dat2 = 0
            if (l.le.nlayer_real) then
               dat2 = summa1_struc(n)%progStruct%gru(iGRU)%hru(iHRU)%  &
                      var(iLookPROG%mLayerVolFracWat)%dat(l)
            endif   ! l < nlayer_real
            call LIS_diagnoseSurfaceOutputVar(n,iGRU,LIS_MOC_MLAYERVOLFRACWAT,&
                  value = dat2,                                        &
                  vlevel=l, unit="-", direction="-",                   &
                  surface_type = LIS_rc%lsm_index)
         enddo ! l, layers

      enddo ! iHRU
   enddo ! iGRU

   deallocate(nSnowData)
   deallocate(nSoilData)

 end subroutine summa2_output

