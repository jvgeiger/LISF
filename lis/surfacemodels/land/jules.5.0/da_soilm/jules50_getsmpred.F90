!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: jules50_getsmpred
! \label{jules50_getsmpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 23Apr2018: Mahdi Navari: Modified for JULES 5.0 
!
! !INTERFACE:
subroutine jules50_getsmpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use jules50_lsmMod
  use jules50_dasoilm_Mod
  !use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use jules_soil_mod, only:  dzsoil 
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k 
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))

!
! !DESCRIPTION:
!
!  Returns the Soil moisture obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  real                   :: obs_tmp
  integer                :: i,t,m,gid,kk
  real                   :: inputs_tp(6), sm_out
  character*50           :: units_tp(6)
  real                   :: smc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                   :: tmp1(LIS_rc%nensem(n)),tmp2(LIS_rc%nensem(n)),tmp3(LIS_rc%nensem(n))

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
! apply the unit conversion. 
! Hx must be the same. we just need to change the unit of Hx. The Unit of Hx 
! from JULES is [kg/m2] therefore we need to change the unit to
! be consistent with that of y [m3/m3]. 
!     smc1(t) = jules50_struc(n)%jules50(t)%smcl_soilt(1) * 1/LIS_sfmodel_struc(n)%lyrthk(1)*1/1000  ! [kg/m2]*1/m*1/(kg/m3) --> [m3/m3] 

!MN: Eric changed the unit of LIS_sfmodel_struc%lyrthk so LIS_sfmodel_struc%lyrthk is in cm, we
!use dzsoil instead.

 smc1(t) = jules50_struc(n)%jules50(t)%smcl_soilt(1) * 1/dzsoil(1)*1/1000  ! [kg/m2]*1/m*1/(kg/m3) --> [m3/m3] 

  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc1,&
       obs_pred)




#if 0 
  tmp1 = LIS_rc%udef
  tmp2 = LIS_rc%udef
  tmp3 = LIS_rc%udef
   do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
    if (i.eq.1376)then
        do m=1,LIS_rc%nensem(n)
           t = (i-1)*LIS_rc%nensem(n)+m
                      tmp1(m) = jules50_struc(n)%jules50(t)%smcl_soilt(1)
	              tmp2(m) = smc1(t)
                      tmp3(m) = obs_pred(i,m)

	enddo
	             ! tmp3 = obs_pred(i)
     endif
   enddo
   print*, 'jules_predicted obs, obs_pred obs sapce, obs_pred meas. space', &
			sum(tmp1)/LIS_rc%nensem(n), &
			sum(tmp2)/LIS_rc%nensem(n), &
			sum(tmp3)/LIS_rc%nensem(n)!, &
                        ! LIS_sfmodel_struc(n)%lyrthk(1)
   print*, 'obs_pred meas. space', &
			tmp3
#endif



end subroutine jules50_getsmpred
