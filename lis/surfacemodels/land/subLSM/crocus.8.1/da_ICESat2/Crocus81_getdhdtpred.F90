!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: Crocus81_getdhdtpred
! \label{Crocus81_getdhdtpred}
!
! !REVISION HISTORY:
! 8 Jan  2024 : Mahdi Navari ;Initial Specification for ICESat2 dhdt da
!
! !INTERFACE:
subroutine Crocus81_getdhdtpred(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use noahmp401_lsmMod
  use LIS_DAobservationsMod
  use Crocus81_dhdt_DAlogMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  real                   :: dh(LIS_rc%npatch(n,LIS_rc%lsm_index))
!EOP

! !DESCRIPTION:
!  This routine computes the obspred ('Hx') term for assimilation
!  instances.

  integer                :: t

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     dh(t) = Crocus81pred_struc(n)%model_dh(t) !Crocus81pred_struc(n)%model_dh(2,t) - Crocus81pred_struc(n)%model_dh(1,t) 
  enddo

  !Crocus81pred_struc(n)%model_dh(1,t) = 0.0
  !Crocus81pred_struc(n)%model_dh(1,t) = Crocus81pred_struc(n)%model_dh(2,t)

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       dh,&
       obs_pred)
  
end subroutine Crocus81_getdhdtpred

