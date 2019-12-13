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
!  June 2019: James Geiger 
! 6 Aug 2019: Zhuo Wang, added support for summa2
!
! !INTERFACE:
subroutine summa2_output(n, summa1_struc)
! !USES:
   use LIS_coreMod,    only : LIS_rc, LIS_surface
   use lis_historymod
   use LIS_histDataMod
   use nrtype,         only : i4b
   use summa_type,     only : summa1_type_dec ! master summa data type
   use globalData,     only : gru_struc  ! gru-hru mapping structures

   ! named variables
   use var_lookup,     only : iLookFLUX  ! named variables for flux data structure
   use var_lookup,     only : iLookFORCE ! named variables for forcing data structure
   use var_lookup,     only : iLookPROG  ! named variables for local state variables
   use var_lookup,     only : iLookDIAG  ! named variables for local diagnostic variables
   use var_lookup,     only : iLookPARAM  ! named variables for local diagnostic variables
   use var_lookup,     only : iLookATTR  ! named variables for local diagnostic variables
   use var_lookup,     only : iLookTYPE  ! named variables for local diagnostic variables

   implicit none
! !ARGUMENTS: 
   integer, intent(in)                 :: n
   type(summa1_type_dec),intent(inout) :: summa1_struc  ! master summa data structure
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

   real                             :: dat
   integer(i4b)                     :: iGRU, iHRU, nGRU
! ---------------------------------------------------------------------------------------
! associate to elements in the data structure
  summaVars: associate(&
 
 ! primary data structures
 timeStruct           => summa1_struc%timeStruct  , & ! x%var(:)                   -- model time data
 forcStruct           => summa1_struc%forcStruct  , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
 indxStruct           => summa1_struc%indxStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
 progStruct           => summa1_struc%progStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
 diagStruct           => summa1_struc%diagStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
 fluxStruct           => summa1_struc%fluxStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
 bvarStruct           => summa1_struc%bvarStruct  , & ! x%gru(:)%var(:)%dat        -- basin-average variables

 ! miscellaneous variables
 nGRU                 => summa1_struc%nGRU        , & ! number of grouped response units
 nHRU                 => summa1_struc%nHRU          & ! number of global hydrologic response units
 
) ! assignment to variables in the data structures
! ---------------------------------------------------------------------------------------

   do iGRU=1,nGRU
      do iHRU=1,gru_struc(iGRU)%hruCount
         if ( gru_struc(iGRU)%hruCount /= 1 ) then
            print*,'GREP: iGRU, hruCount', iGRU, gru_struc(iGRU)%hruCount
         endif

         dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SWNET, &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="DN",           &
            surface_type=LIS_rc%lsm_index)

         dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFLUX%scalarLWNetGround)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n, iGRU ,LIS_MOC_LWNET, &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="DN",           &
            surface_type=LIS_rc%lsm_index)

! Zhuo Wang added new SUMMA output variables on 08/06/2019
          dat = summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookPROG%scalarCanopyIce)%dat(1)
   
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_CANOPICE,    &
            value=dat,                                        &
            vlevel=1, unit="kg m-2", direction="-",           &
            surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookPROG%scalarCanopyLiq)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_CANOPLIQ,    &
            value=dat,                                        &
            vlevel=1, unit="kg m-2", direction="-",           &
            surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%     &
                var(iLookPROG%scalarCanopyWat)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_CANOPWAT,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2", direction="-",           &
          surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookPROG%scalarSWE)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SWE,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%     &
                var(iLookPROG%scalarAquiferStorage)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_AQUIFERSTORAGE,    &
               value=dat,                                        &
               vlevel=1, unit="m", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%diagStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookDIAG%scalarSoilCompress)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SOILCOMPRESS,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%diagStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookDIAG%scalarTotalSoilLiq)%dat(1)

          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SOILLIQ,    &
               value=dat,                                      &
               vlevel=1, unit="kg m-2", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%diagStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookDIAG%scalarTotalSoilIce)%dat(1)
           
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SOILICE,    &
               value=dat,                                      &
               vlevel=1, unit="kg m-2", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%diagStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookDIAG%scalarTotalSoilWat)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SOILMOIST,    &
               value=dat,                                      &
               vlevel=1, unit="kg m-2", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_CANOPAIRNETFLUX,    &
               value=dat,                                      & 
               vlevel=1, unit="W m-2", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_CANOPNETFLUX,    &
               value=dat,                                      &
               vlevel=1, unit="W m-2", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_GROUNDNETFLUX,    &
               value=dat,                                      &
               vlevel=1, unit="W m-2", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SAV,    &
               value=dat,                                      &
               vlevel=1, unit="W m-2", direction="IN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SAG,    &
               value=dat,                                      & 
               vlevel=1, unit="W m-2", direction="IN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarLWNetCanopy)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_IRC,    &
               value=dat,                                      &              
               vlevel=1, unit="W m-2", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarLWNetGround)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_IRB,    &
               value=dat,                                      &              
               vlevel=1, unit="W m-2", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarSenHeatTotal)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QH,    &
               value=dat,                                      &      
               vlevel=1, unit="W m-2", direction="UP",           &    
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarSenHeatCanopy)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SHC,    &
               value=dat,                                      &
               vlevel=1, unit="W m-2", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarSenHeatGround)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SHB,    &
               value=dat,                                      &
               vlevel=1, unit="W m-2", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarLatHeatTotal)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QLE,    &
               value=dat,                                      & 
               vlevel=1, unit="W m-2", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_FCEV,    &
               value=dat,                                      &              
               vlevel=1, unit="W m-2", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_FCTR,    &
               value=dat,                                      &              
               vlevel=1, unit="W m-2", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarLatHeatGround)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_EVG,    &
               value=dat,                                      &              
               vlevel=1, unit="W m-2", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_CANOPADVECT,    &
               value=dat,                                      &              
               vlevel=1, unit="W m-2", direction="DN",           &
               surface_type=LIS_rc%lsm_index)
          
          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_GROUNDADVECT,    &
               value=dat,                                      &     
               vlevel=1, unit="W m-2", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopySublimation)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SUBCANOP,    &
               value=dat,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarSnowSublimation)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SUBSNOW,    &
               value=dat,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopyTranspiration)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_TVEG,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopyEvaporation)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_ECANOP,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="UP",           &
               surface_type=LIS_rc%lsm_index)
           
          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarGroundEvaporation)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_ESOIL,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
                var(iLookFLUX%scalarThroughfallSnow)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QST,    &
                value=dat,                                      &
                vlevel=1, unit="kg m-2 s-1", direction="-",           &
                surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarThroughfallRain)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_THROUGHRAIN,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SNOWUPLOADING,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_CANOP_LIQDRAINAGE,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarCanopyMeltFreeze)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_CANOP_MELT_FREEZE,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarRainPlusMelt)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_RAIN_MELT,    &
               value=dat,                                        &
               vlevel=1, unit="m s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
                var(iLookFLUX%scalarSurfaceRunoff)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QS,    &
                value=dat * 1000.,                                   &
                vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
                surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarSoilBaseflow)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SOIL_BASEFLOW,    &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)
 
           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarSoilDrainage)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SOIL_DRAINAGE,    &
               value=dat * 1000.,                                      &
               vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
               surface_type=LIS_rc%lsm_index)

           dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarAquiferRecharge)%dat(1)
           call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_AQUIFER_RECHARGE,    &
               value=dat,                                        &
               vlevel=1, unit="m s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarAquiferTranspire)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_AQUIFER_TRANSPIRE,    &
               value=dat,                                        &             
               vlevel=1, unit="m s-1", direction="-",           &
               surface_type=LIS_rc%lsm_index)
          
          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
                var(iLookFLUX%scalarAquiferBaseflow)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_AQUIFER_BASEFLOW,    &
                value=dat,                                        &    
                vlevel=1, unit="m s-1", direction="-",           &
                surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarTotalET)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_EVAP,    &    
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="UP",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
                var(iLookFLUX%scalarTotalRunoff)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QSTOT,    &
                value=dat * 1000.,                                     & 
                vlevel=1, unit="kg m-2 s-1", direction="OUT",           &
                surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarNetRadiation)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_RNET,    &
               value=dat,                                        &
               vlevel=1, unit="W m-2", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarRainfall)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_RAINF,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookFLUX%scalarSnowfall)%dat(1)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SNOWF,    &
               value=dat,                                        &
               vlevel=1, unit="kg m-2 s-1", direction="DN",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%typeStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookTYPE%hruId)

          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_HRUID, &
               value=dat,                                        &
               vlevel=1, unit="-", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookATTR%elevation)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_ELEVATION, &
               value=dat,                                        &
               vlevel=1, unit="m", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%mparStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookPARAM%albedoMax)%dat(1)

          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_MXSNALBEDO, &
               value=dat,                                        &
               vlevel=1, unit="-", direction="-",           &
               surface_type=LIS_rc%lsm_index)

          dat = summa1_struc%typeStruct%gru(iGRU)%hru(iHRU)%     &
               var(iLookTYPE%vegTypeIndex)
          call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_GREENNESS, &
               value=dat,                                        &
               vlevel=1, unit="-", direction="-",           &
               surface_type=LIS_rc%lsm_index)
   
      enddo    ! iHRU
   enddo       ! iGRU

end associate summaVars

end subroutine summa2_output
