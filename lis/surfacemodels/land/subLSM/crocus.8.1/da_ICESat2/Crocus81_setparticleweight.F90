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
! !ROUTINE: Crocus81_setparticleweight
! \label{Crocus81_setparticleweight}
!
! !REVISION HISTORY:
! 19 Jan: Mahdi Navari; Initial Specification
!
!
! !INTERFACE:
subroutine Crocus81_setparticleweight(n, LIS_LSM_particle_weight)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  !use LIS_snowMod, only : LIS_snow_struc
  use LIS_logMod, only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_timeMgrMod
  use Crocus81_lsmMod
  use Crocus81_dhdt_DAlogMod
  use LIS_numerRecipesMod, only : LIS_rand_func

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LIS_LSM_particle_weight
! 
! !DESCRIPTION:
! 
!  This routine assigns the snow progognostic variables to noah's
!  model space. The state vector consists of total SWE and snow depth. 
!  This routine also updates other model prognostics (snice, snliq,
!  snow thickness, snow temperature) based on the update. 
! 
!EOP
  type(ESMF_Field)       :: ParticleWeightField
  !type(ESMF_Field)       :: snodField
  real, pointer          :: Pw(:)
  !real, pointer          :: snod(:)
  !real                   :: dsneqv,dsnowh
  integer                :: t
  integer                :: status
  integer                :: N_state, N_ens ,state_size, n_l
  integer                :: n_e, i, ii, jj, jk, n_id
  real                   :: ran_face
  real, allocatable      :: Pw_cumsum(:)
  real, allocatable      :: new_particle_idx(:), ens_id_SIR(:)
  real, allocatable      :: P_w_curr_ts(:)
  real                   :: current_p
  real                   :: iii
  type(ESMF_Time)        :: ATL15_StartTime
  type(ESMF_time)        :: currTime
  type(ESMF_TimeInterval):: obs_interval
  type(ESMF_TimeInterval):: tw
  integer                :: yr,mo,da,hr,mn,ss,ms
  REAL*8, allocatable    :: tmp_SNOWSWE(:,:)         
  REAL*8, allocatable    :: tmp_SNOWRHO(:,:)         
  REAL*8, allocatable    :: tmp_SNOWHEAT(:,:)        
  REAL*8, allocatable    :: tmp_SNOWALB(:) !,:)            
  REAL*8, allocatable    :: tmp_SNOWGRAN1(:,:)       
  REAL*8, allocatable    :: tmp_SNOWGRAN2(:,:)       
  REAL*8, allocatable    :: tmp_SNOWHIST(:,:)        
  REAL*8, allocatable    :: tmp_SNOWAGE(:,:)
  REAL*8, allocatable    :: tmp_SNOWLIQ(:,:)         
  REAL*8, allocatable    :: tmp_SNOWTEMP(:,:)        
  REAL*8, allocatable    :: tmp_SNOWDZ(:,:)          
  !REAL*8, allocatable    :: tmp_ALB(:,:)
  !REAL*8                 :: tmp_THRUFAL           
  REAL*8 , allocatable   :: tmp_GRNDFLUX(:)           
  REAL   , allocatable   :: tmp_SNDRIFT(:)            
  REAL   , allocatable   :: tmp_RI_n(:)               
  !REAL*8                 :: tmp_EMISNOW           
  REAL   , allocatable   :: tmp_CDSNOW(:)             
  REAL   , allocatable   :: tmp_USTARSNOW(:)          
  REAL   , allocatable   :: tmp_CHSNOW(:)
  REAL*8 , allocatable   :: tmp_SNOWMAK_dz(:)
  REAL*8 , allocatable   :: tmp_TG(:)
  REAL*8 , allocatable   :: tmp_XWGI(:)               
  REAL*8 , allocatable   :: tmp_XWG(:)    
  INTEGER                :: N_obs_in_tw
  !integer                :: yr1, mo1, da1, hr1, mn1
  integer                :: yr2, mo2, da2, hr2, mn2
  integer                :: yr3, mo3, da3, hr3, mn3
  integer                :: yr4, mo4, da4, hr4, mn4
  integer                :: yr5, mo5, da5, hr5, mn5
  type(ESMF_Time)        :: twObsTime1,twObsTime2,twObsTime3,twObsTime4,twObsTime5
  real, allocatable      :: Pw_combined(:)
  real, allocatable      :: Pw_norm(:)  


  call ESMF_StateGet(LIS_LSM_particle_weight,"ParticleWeight",ParticleWeightField,rc=status) 
  ! NOTE: "ParticleWeight" is hard coded in LIS_lsmMod.F90. LIS does not read that from attribute file.
  call LIS_verify(status)

  call ESMF_FieldGet(ParticleWeightField,localDE=0,farrayPtr=Pw,rc=status)
  call LIS_verify(status)

  !2018-10-01 22:30:00.00   2019-01-01 06:00:00.00   2019-04-02 13:30:00.00   2019-07-02 21:00:00.00
  !2019-10-02 04:30:00.00   2020-01-01 12:00:00.00   2020-04-01 19:30:00.00   2020-07-02 03:00:00.00
  !2020-10-01 10:30:00.00   2020-12-31 18:00:00.00   2021-04-02 01:30:00.00   2021-07-02 09:00:00.00

  !CROCUS81_struc(n)%NumOfObsPerAssimWindow

  call ESMF_ClockGet(LIS_clock, currTime = currTime, rc=status)

  call ESMF_TimeSet(ATL15_StartTime, yy = 2018, &
            mm = 10, &
            dd = 1, &
            h  = 12, &
            m  = mn,&
            s  = ss, &
            calendar = LIS_calendar, &
            rc = status)

  call ESMF_TimeGet(LIS_twStartTime, yy = yr, &
            mm = mo, &
            dd = da, &
            h  = hr, &
            m  = mn,&
            s  = ss, &
            calendar = LIS_calendar, &
            rc = status)
  if (mod(LIS_rc%twInterval,7889400.0) .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] The DA time window interval must be a multiple of the '
      write(LIS_logunit,*) '[ERR] obs interval of 7889400.0 [s]. '
      write(LIS_logunit,*) '[ERR] Program stopping.... '
      call LIS_endrun()
  endif 

  call ESMF_TimeintervalSet(tw,s=nint(LIS_rc%twInterval),rc=status)
  call ESMF_TimeintervalSet(obs_interval,s=7889400,rc=status)
  
  if (currTime .eq. ATL15_StartTime ) then        
     !reset timewondow
     call LIS_resetClockForPBSTimeWindow(LIS_rc)
  endif 

  !if ( (LIS_twStartTime .gt. (ATL15_StartTime - obs_interval/2)) .and. &
  !     (LIS_twStopTime  .gt.  ATL15_StartTime) ) then
  if (currTime .ge. ATL15_StartTime ) then 

    N_obs_in_tw = floor((LIS_twStopTime - LIS_twStartTime) / obs_interval) -1 ! Num. of dhdt obs  

    if (N_obs_in_tw .ge. 2) then  
      if (LIS_twStartTime .le. ATL15_StartTime) then ! assume DA simulation starts around the ATL15_StartTime 
          twObsTime1 = ATL15_StartTime    ! ATL15_StartTime means dh start time not dhdt start time 
      else
          twObsTime1 = ATL15_StartTime + tw   !- obs_interval
      endif
     
     SELECT CASE (N_obs_in_tw)
     CASE (2)
        twObsTime2  = twObsTime1 + obs_interval
        call ESMF_TimeGet(twObsTime2,yy=yr2,mm=mo2,dd=da2, h=hr2, m=mn2 ,calendar=LIS_calendar,&
               rc=status)

     CASE (3)
        twObsTime2  = twObsTime1 + obs_interval
        twObsTime3  = twObsTime1 + 2*obs_interval 
        call ESMF_TimeGet(twObsTime2,yy=yr2, mm=mo2, dd=da2, h=hr2, m=mn2 ,calendar=LIS_calendar,&
               rc=status)
        call ESMF_TimeGet(twObsTime3,yy=yr3, mm=mo3, dd=da3, h=hr3, m=mn3 ,calendar=LIS_calendar,&
               rc=status)
     CASE (4)
        twObsTime2  = twObsTime1 + obs_interval
        twObsTime3  = twObsTime1 + 2*obs_interval
        twObsTime4  = twObsTime1 + 3*obs_interval
        call ESMF_TimeGet(twObsTime2,yy=yr2, mm=mo2, dd=da2, h=hr2, m=mn2 ,calendar=LIS_calendar,&
               rc=status)
        call ESMF_TimeGet(twObsTime3,yy=yr3, mm=mo3, dd=da3, h=hr3, m=mn3 ,calendar=LIS_calendar,&
               rc=status)
        call ESMF_TimeGet(twObsTime4,yy=yr4, mm=mo4, dd=da4, h=hr4, m=mn4 ,calendar=LIS_calendar,&
               rc=status)

     CASE (5)
        twObsTime2  = twObsTime1 + obs_interval
        twObsTime3  = twObsTime1 + 2*obs_interval
        twObsTime4  = twObsTime1 + 3*obs_interval
        twObsTime5  = twObsTime1 + 4*obs_interval
        call ESMF_TimeGet(twObsTime2,yy=yr2, mm=mo2, dd=da2, h=hr2, m=mn2 ,calendar=LIS_calendar,&
               rc=status)
        call ESMF_TimeGet(twObsTime3,yy=yr3, mm=mo3, dd=da3, h=hr3, m=mn3 ,calendar=LIS_calendar,&
               rc=status)
        call ESMF_TimeGet(twObsTime4,yy=yr4, mm=mo4, dd=da4, h=hr4, m=mn4 ,calendar=LIS_calendar,&
               rc=status)
        call ESMF_TimeGet(twObsTime5,yy=yr5, mm=mo5, dd=da5, h=hr5, m=mn5 ,calendar=LIS_calendar,&
               rc=status)

     CASE (6)
        write(LIS_logunit,*)'[ERR] Please reduce the assimilation window to one year or less'
        write(LIS_logunit,*)'[ERR] program stopping ...'
        call LIS_endrun 
     END SELECT

     N_ens = LIS_rc%nensem(n)
     if(.not.allocated(Crocus81pred_struc(n)%Pw_combined)) then 
        allocate(Crocus81pred_struc(n)%Pw_combined(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     endif 
     if(.not.allocated(Crocus81pred_struc(n)%ens_id_SIR)) then
        allocate(Crocus81pred_struc(n)%ens_id_SIR(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     endif   
 
     allocate(Pw_combined(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     allocate(Pw_norm(LIS_rc%npatch(n,LIS_rc%lsm_index)))
 
     !account for all observations in this time window, 
     !by taking the product of the weights for
     !each ensemble member derived from all assimilated observations 
    
     !if (currTime .le. LIS_twStartTime) then
     if (currTime .eq. twObsTime2) then ! if (N_obs_in_tw .eq. 2) then  
        Crocus81pred_struc(n)%Pw_combined(:) = 0.0
        Crocus81pred_struc(n)%Pw_combined(:) = Pw
        Crocus81pred_struc(n)%ens_id_SIR(:) = 0   

     else if ( (currTime .eq. twObsTime3 ) .or.  &
               (currTime .eq. twObsTime4 ) .or.  &
               (currTime .eq. twObsTime5 )  ) then  
        Pw_combined = Crocus81pred_struc(n)%Pw_combined*Pw
        !normalize into probability "weights" that account for all observations
        if(sum(Pw_combined).ne.0) then
           Pw_norm=Pw_combined/sum(Pw_combined)
           Crocus81pred_struc(n)%Pw_combined = Pw_norm 
        endif
     
     endif 

     !if (LIS_twStopTime .ge. ATL15_StartTime ) then  
     !if (currTime .le. LIS_twStopTime ) then      
     if ((currTime .eq. twObsTime2 ) .or.  &  
         (currTime .eq. twObsTime3 ) .or.  &
         (currTime .eq. twObsTime4 ) .or.  &
         (currTime .eq. twObsTime5 )  ) then  
 
     if(.not.allocated(P_w_curr_ts)) then
        allocate(P_w_curr_ts(N_ens))
     endif

     if(.not.allocated(new_particle_idx)) then
        allocate(new_particle_idx(N_ens))
     endif

     if(.not.allocated(ens_id_SIR)) then
        allocate(ens_id_SIR(N_ens))
     endif

        do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
           P_w_curr_ts = Crocus81pred_struc(n)%Pw_combined(((i-1)*N_ens+1):((i-1)*N_ens+N_ens))
           
           !sequential importance resampling (SIR)  
           do n_e=1,N_ens
              Pw_cumsum(n_e)=sum(P_w_curr_ts(1:n_e))
           end do

           do n_e=1,N_ens
              !generate random number between 0-1 (random dice roll)
              call LIS_rand_func(1,ran_face)
              !select first insance of the random number 
              !being exceeded by the posterior cumulative
              !sumation of the posterior PDF

              !vector of 0's and 1's. 1's begin in the first 
              !instance of ran_face exceeding the Pw_cumsum vector
              do jk=1,N_ens
                 if (Pw_cumsum(jk).le.ran_face) then
                    new_particle_idx(jk)=0
                 else if (Pw_cumsum(jk).gt.ran_face) then
                    new_particle_idx(jk)=1
                 end if
              end do

              !select the first instance of the logical statement 
              !occuring (or the first non-zero instance) and resample
              !from this instance
              iii=0

              do ii=1,N_ens
                 current_p=new_particle_idx(ii)
                 if (current_p.gt.0 .AND. iii.eq.0) then
                    !do jj=1,N_state
                    ! Note: We only assimilate elevation change. Since we update the ensemble weights all 
                    ! states get the same weight.
                    ! Unlike other DA methods, We can not simultaneously assimilate different observations. 
                    ! Because the magnitude of (obs-state) affects the results of the likelihood function. 
                    ! Margulis has done this using CDF matching. He matched the magnitude of other states 
                    ! and observations to one of them using CDF matching. 
                       !updated_state(jj,n_e)=State_incr(jj,i)
                       !ens_id_SIR(jj,n_e) = i
                       ens_id_SIR(n_e) = ii
                    !end do
                    iii=iii+1
                 end if
              end do

           end do ! N_ens
           Crocus81pred_struc(n)%ens_id_SIR(((i-1)*N_ens+1):((i-1)*N_ens+N_ens)) = ens_id_SIR
        enddo ! i
         
     else if (currTime .eq. LIS_twStopTime) then
        !reset timewondow
        call LIS_resetClockForPBSTimeWindow(LIS_rc)
        !apply the update 

        allocate( tmp_SNOWSWE  (LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWRHO  (LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWHEAT (LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWGRAN1(LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWGRAN2(LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWHIST (LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWAGE  (LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWALB  (LIS_rc%nensem(n) ) )
        allocate( tmp_SNOWLIQ  (LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWTEMP (LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNOWDZ   (LIS_rc%nensem(n), CROCUS81_struc(n)%nsnow ) )
        allocate( tmp_SNDRIFT  (LIS_rc%nensem(n)))  
        allocate( tmp_RI_n     (LIS_rc%nensem(n)))
        allocate( tmp_CDSNOW   (LIS_rc%nensem(n)))
        allocate( tmp_USTARSNOW(LIS_rc%nensem(n)))
        allocate( tmp_CHSNOW   (LIS_rc%nensem(n)))
        allocate( tmp_SNOWMAK_dz(LIS_rc%nensem(n)))
        allocate( tmp_TG       (LIS_rc%nensem(n)))
        allocate( tmp_XWGI     (LIS_rc%nensem(n)))
        allocate( tmp_XWG      (LIS_rc%nensem(n)))


        do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)

           ! save state variables in temporary variables 
           do n_e=1,N_ens
              do n_l=1,CROCUS81_struc(n)%nsnow
              tmp_SNOWSWE(n_e,n_l)  = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWSWE(n_l)
              tmp_SNOWRHO(n_e,n_l)  = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWRHO(n_l)
              tmp_SNOWHEAT(n_e,n_l) = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWHEAT(n_l)
              tmp_SNOWGRAN1(n_e,n_l)= CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWGRAN1(n_l)
              tmp_SNOWGRAN2(n_e,n_l)= CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWGRAN2(n_l)
              tmp_SNOWHIST(n_e,n_l) = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWHIST(n_l)
              tmp_SNOWAGE(n_e,n_l)  = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWAGE(n_l)
              tmp_SNOWLIQ(n_e,n_l)  = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWLIQ(n_l)
              tmp_SNOWTEMP(n_e,n_l) = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWTEMP(n_l)
              tmp_SNOWDZ(n_e,n_l)   = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWDZ(n_l)
              enddo
              tmp_SNOWALB(n_e)      = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWALB 
              tmp_GRNDFLUX(n_e)     = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%GRNDFLUX
              tmp_SNDRIFT(n_e)      = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNDRIFT
              tmp_RI_n(n_e)         = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%RI_n
              tmp_CDSNOW(n_e)       = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%CDSNOW
              tmp_USTARSNOW(n_e)    = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%USTARSNOW
              tmp_CHSNOW(n_e)       = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%CHSNOW
              tmp_SNOWMAK_dz(n_e)   = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%SNOWMAK_dz
              tmp_TG(n_e)           = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%TG
              tmp_XWGI(n_e)         = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%XWGI
              tmp_XWG(n_e)          = CROCUS81_struc(n)%crocus81((i-1)*N_ens+n_e)%XWG

           enddo

           ! Now, discard the un-survived ensembles and replace them with survived ensembles 
           ! using sequential importance resampling (SIR)

           ens_id_SIR = Crocus81pred_struc(n)%ens_id_SIR(((i-1)*N_ens+1):((i-1)*N_ens+N_ens)) 
           do n_e=1,N_ens
              n_id = ens_id_SIR(n_e)
              do n_l=1,CROCUS81_struc(n)%nsnow
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWSWE(n_l)   = tmp_SNOWSWE(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWRHO(n_l)   = tmp_SNOWRHO(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWHEAT(n_l)  = tmp_SNOWHEAT(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWGRAN1(n_l) = tmp_SNOWGRAN1(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWGRAN2(n_l) = tmp_SNOWGRAN2(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWHIST(n_l)  = tmp_SNOWHIST(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWAGE(n_l)   = tmp_SNOWAGE(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWLIQ(n_l)   = tmp_SNOWLIQ(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWTEMP(n_l)  = tmp_SNOWTEMP(n_id,n_l)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWDZ(n_l)    = tmp_SNOWDZ(n_id,n_l)
              enddo
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWALB      = tmp_SNOWALB(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%GRNDFLUX     = tmp_GRNDFLUX(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNDRIFT      = tmp_SNDRIFT(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%RI_n         = tmp_RI_n(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%CDSNOW       = tmp_CDSNOW(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%USTARSNOW    = tmp_USTARSNOW(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%CHSNOW       = tmp_CHSNOW(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%SNOWMAK_dz   = tmp_SNOWMAK_dz(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%TG           = tmp_TG(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%XWGI         = tmp_XWGI(n_id)
              CROCUS81_struc(n)%crocus81(((i-1)*N_ens+1):((i-1)*N_ens+n_e))%XWG          = tmp_XWG(n_id)
           enddo
        enddo
        deallocate( tmp_SNOWSWE )
        deallocate( tmp_SNOWRHO )
        deallocate( tmp_SNOWHEAT )
        deallocate( tmp_SNOWGRAN1 )
        deallocate( tmp_SNOWGRAN2 )
        deallocate( tmp_SNOWHIST )
        deallocate( tmp_SNOWAGE )
        deallocate( tmp_SNOWALB )
        deallocate( tmp_SNOWLIQ )
        deallocate( tmp_SNOWTEMP )
        deallocate( tmp_SNOWDZ )  
        deallocate( tmp_SNDRIFT ) 
        deallocate( tmp_RI_n )
        deallocate( tmp_CDSNOW )
        deallocate( tmp_USTARSNOW )
        deallocate( tmp_CHSNOW )
        deallocate( tmp_SNOWMAK_dz )
        deallocate( tmp_TG )
        deallocate( tmp_XWGI )
        deallocate( tmp_XWG )

        deallocate(Crocus81pred_struc(n)%Pw_combined)
        deallocate(Crocus81pred_struc(n)%ens_id_SIR)
        deallocate(Pw_combined)
        deallocate(Pw_norm)
        deallocate(ens_id_SIR)
        deallocate(new_particle_idx)
        deallocate(P_w_curr_ts) 
     endif 
    endif  
  endif 


end subroutine Crocus81_setparticleweight


