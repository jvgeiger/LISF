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
!
! !INTERFACE:
subroutine summa2_output(n, summa1_struc)
! !USES:
   use LIS_coreMod,    only : LIS_rc, LIS_surface
   use lis_historymod
   use LIS_histDataMod
   use nrtype,         only : i4b
   use summa_type,     only : summa1_type_dec ! master summa data type
   use globalData,     only : gru_struc ! gru-hru mapping structures
   use var_lookup,     only : iLookFLUX! named variables for flux data structure
   use var_lookup,     only : iLookFORCE!named variables for flux data structure

   implicit none
! !ARGUMENTS: 
   integer, intent(in)                 :: n
   type(summa1_type_dec),intent(inout) :: summa1_struc
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

!<debug -- jim testing>
   real                             :: dat
   integer(i4b)                     :: iGRU, iHRU, nGRU
   integer                          :: jim_count=0
   character(len=3)                 :: sjim_count
   real,allocatable,dimension(:,:)  :: datarray
   logical                          :: first_time = .true.
!</debug -- jim testing>

#if 1
   nGRU = summa1_struc%nGRU
   print*,'GREP: nGRU',nGRU
!<debug -- jim testing>
   !write(unit=667,fmt=*) 'GREP: ----------------------------------'
   allocate(datarray(LIS_rc%gnc(n),LIS_rc%gnr(n)))
   datarray = LIS_rc%udef
   do iGRU=1,nGRU
      do iHRU=1,gru_struc(iGRU)%hruCount
         datarray(LIS_surface(n,LIS_rc%lsm_index)%tile(iGRU)%col,      &
            LIS_surface(n,LIS_rc%lsm_index)%tile(iGRU)%row) =    &
            summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%            &
            var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1)
      enddo
   enddo
   jim_count=jim_count+1
   write(sjim_count, '(i3.3)') jim_count
   open(unit=666,file='jim'//sjim_count//'.bin',access='direct',&
      recl=LIS_rc%gnc(1)*LIS_rc%gnr(1)*4)
   call LIS_writevar_bin(666, n, datarray, 1)
!</debug -- jim testing>
   do iGRU=1,nGRU
      do iHRU=1,gru_struc(iGRU)%hruCount
         if ( gru_struc(iGRU)%hruCount /= 1 ) then
            print*,'GREP: iGRU, hruCount', iGRU, gru_struc(iGRU)%hruCount
         endif

         dat = summa1_struc%forcStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFORCE%SWRadAtm)
!<debug -- jim testing>
         !write(unit=667,fmt=*) 'GREP: swradatm dat ', iGRU, dat
!</debug -- jim testing>
         call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QF,    &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="S2L",          &
            surface_type=LIS_rc%lsm_index)

         dat = summa1_struc%forcStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFORCE%LWRadAtm)
         call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QV,    &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="S2V",          &
            surface_type=LIS_rc%lsm_index)

         dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_SWNET, &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="DN",           &
            surface_type=LIS_rc%lsm_index)

!<debug -- jim testing>
         !write(unit=667,fmt=*) 'GREP: lwnetgrd dat ', iGRU, dat
!</debug -- jim testing>
         dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFLUX%scalarLWNetGround)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n, iGRU ,LIS_MOC_LWNET, &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="DN",           &
            surface_type=LIS_rc%lsm_index)

         dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFLUX%scalarLatHeatGround)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QLE,   &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="UP",           &
            surface_type=LIS_rc%lsm_index)

         dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFLUX%scalarSenHeatGround)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QH,    &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="UP",           &
            surface_type=LIS_rc%lsm_index)

         dat = summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%     &
            var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1)
         call LIS_diagnoseSurfaceOutputVar(n, iGRU, LIS_MOC_QG,    &
            value=dat,                                        &
            vlevel=1, unit="W m-2", direction="DN",           &
            surface_type=LIS_rc%lsm_index)
      enddo
   enddo
#endif
end subroutine summa2_output
