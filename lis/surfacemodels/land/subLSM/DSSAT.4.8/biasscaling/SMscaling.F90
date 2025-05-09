!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
! !ROUTINE: SMscaling
! \label{SMscaling}
!
! !REVISION HISTORY:
!  07 Oct 2021: Pang-Wei Liu; Subroutine for intra model SM bias correction

! !INTERFACE:
subroutine SMscaling(n, grid_index, lsm_sm)

! !USES:
   use LIS_coreMod
   use SMscaling_Mod, only: SMscaling_struc

   implicit none
! !ARGUMENTS:
   integer  :: n
   integer                :: grid_index
   real                   :: lsm_sm(4) !This is the SM from LSM
   integer                :: k, l
!-------------------------------------------------------------------------
!  Transform LSM soil moisture to DSSAT climatology using a scaling algorithms
!-------------------------------------------------------------------------

k=1 !Only have 1 model
if (SMscaling_struc(n)%bcmethod .eq. "CDF matching") then
      !if (SMscaling_struc(n)%ntimes .gt. 1 .and. SMscaling_struc(n)%cdf_read_opt .eq. 1) then
         !call LIS_rescale_with_CDF_matching( &
         !   n, k, &
         !   SMscaling_struc(n)%nbins, &
         !   1, &
         !   MAX_SM_VALUE, &
         !   MIN_SM_VALUE, &
         !   SMscaling_struc(n)%dssat_xrange, &
         !   SMscaling_struc(n)%lsm_xrange, &
         !   SMscaling_struc(n)%dssat_cdf1, &
         !   SMscaling_struc(n)%lsm_cdf1, &
         !   sm_current)
      !else
      !    call LIS_rescale_with_CDF_matching( &
      !      n, k, &
      !      SMscaling_struc(n)%nbins, &
      !      SMscaling_struc(n)%ntimes, &
      !      MAX_SM_VALUE, &
      !      MIN_SM_VALUE, &
      !      SMscaling_struc(n)%dssat_xrange, &
      !      SMscaling_struc(n)%lsm_xrange, &
      !      SMscaling_struc(n)%dssat_cdf, &
      !      SMscaling_struc(n)%lsm_cdf, &
      !      sm_current)
      !endif
   elseif(SMscaling_struc(n)%bcmethod.eq."Linear scaling") then
      !  call LIS_rescale_with_linear_scaling(    &
      !       n,                                   &
      !       k,                                   &
      !       SMscaling_struc(n)%nbins,         &
      !       SMscaling_struc(n)%ntimes,        &
      !       SMscaling_struc(n)%lsm_xrange,    &
      !       SMscaling_struc(n)%lsm_cdf,       &
      !       sm_current)
   elseif(SMscaling_struc(n)%bcmethod.eq."Anomaly scaling") then
        DO l=1,SMscaling_struc(n)%bclayer

            call LIS_DSSAT_rescale_with_anomaly(    &
                 n,                                 &
                 grid_index,                        &
                 SMscaling_struc(n)%ntimes,         &
                 SMscaling_struc(n)%lsm_mu(l,:,:),  &
                 SMscaling_struc(n)%dssat_mu(l,:,:),&
                 lsm_sm(l))
        ENDDO
   endif

end subroutine SMscaling

!BOP
! 
! !ROUTINE: LIS_DSSAT_rescale_with_anomaly
! \label{LIS_DSSAT_rescale_with_anomaly}
!
  subroutine LIS_DSSAT_rescale_with_anomaly(&
       n,             &
       ind,           &
       ntimes,        &
       obs_mu,        & !In the LIS-DSSAT case obs is LSM
       model_mu,      & !In the LIS-DSSAT case model is DSSAT
       obs_value)

    use LIS_coreMod

    implicit none

! !ARGUMENTS: 
    integer             :: n
    integer             :: ind
    integer             :: ntimes
    real                :: obs_mu(LIS_rc%ngrid(1),ntimes)
    real                :: model_mu(LIS_rc%ngrid(1),ntimes)
    real                :: obs_value
! !DESCRIPTION: 
! 
!   This routine rescales the input observation data 
!EOP
    integer             :: kk
    real                :: obs_tmp
    real                :: obs_anomaly

    if(ntimes.gt.1) then
       kk = LIS_rc%mo
    else
       kk = 1
    endif
       if(obs_value.ne.-9999.0) then

         obs_anomaly = obs_value - obs_mu(ind, kk)
         obs_tmp = model_mu(ind,kk) + obs_anomaly
          if (obs_tmp.lt.0.01) then
              obs_tmp = 0.01
          endif
          obs_value = obs_tmp
       else
          obs_value = LIS_rc%udef
       endif
end subroutine LIS_DSSAT_rescale_with_anomaly
!------------------------------------------------------------------
