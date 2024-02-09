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
! !ROUTINE: read_ATL15_GrISobs
! \label{read_ATL15_GrISobs}
!
! !REVISION HISTORY: 
!  13 Dec 2023    Mahdi Navari;   Initial Specification
!
! !INTERFACE: 
subroutine read_ATL15_GrISobs(n,k, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use ATL15_GrISobs_module
  use LIS_fileIOMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use ATL15_GrISobs_module
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the ATL15_GrIS height-change observations 
!  The processed data is packaged into an ESMF State
!  for later use within the DA algorithm
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]                index of the nest
!  \item[k]                index of the data assimilation instance
!  \item[OBS\_State]       observations state
!  \item[OBS\_Pert\_State] observations perturbation state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: dhdtfield
  type(ESMF_Time)     :: startdate, currentTime
!  integer             :: iret

  real,    pointer    :: obsl(:)
!  real                :: gracegrid(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
!  real                :: gracegrid_glb(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k))
!  real                :: graceerr(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
!  real                :: graceerr_glb(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k))
!  real                :: gridDesc(6)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: ATL15_GrISobsdir ! BMc, change 200 to LIS_CONST_PATH_LEN
  character(len=LIS_CONST_PATH_LEN) :: name ! BMc, change 200 to LIS_CONST_PATH_LEN
  logical             :: file_exists
!  integer             :: col,row
  logical             :: data_upd
  logical             :: data_upd_flag_local
  logical             :: data_upd_flag(LIS_npes)
  integer             :: fnd
  logical             :: readflag
  logical             :: alarmCheck
  integer             :: status !,ftn
!  real, allocatable   :: ssdev(:)
  integer             :: t,p,r,c
!  integer             :: nc, nr, num_obs
  integer             :: dimid, ncid, varid
  integer             :: leng
  real, allocatable   ::obs_time(:)
  real, allocatable   ::dhdt(:,:,:) ! "time y x" but Fortran reads in reverse order (x,y,time)
  real, allocatable   ::dhdt_local_grid(:,:,:)
  integer             :: yy,mm,dd,h,m,s
  integer             :: doy,ts
  real                :: gmt
  real*8              :: timenow
  real*8              :: start_date, start_date_tmp
  integer             :: offset

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       ATL15_GrISobsdir, rc=status)
  call LIS_verify(status)
  
  call ESMF_AttributeSet(OBS_State,"Data Update Status",&
       .false., rc=status)
  call LIS_verify(status)
  
  call ATL15_GrIS_filename(name,ATL15_GrISobsdir)
  
  inquire(file=name,exist=file_exists)
  if(file_exists) then 
     call ESMF_AttributeSet(OBS_State,"File Status",&
          .true., rc=status)
     call LIS_verify(status)
  else
     call ESMF_AttributeSet(OBS_State,"File Status",&
          .false., rc=status)
     call LIS_verify(status)
     write(LIS_logunit,*)  '[ERR] ATL15 hight change data',trim(name)
     write(LIS_logunit,*)  'is missing. Please check the path and ATL15 file name'
     call LIS_endrun()
  endif


! TODO add a condition to check the CROCUS81_struc(n)%ts == 15mn otherwise the model bypass the alarmCheck


! ATL15 time "days since 2018-01-01". 
! First and second observations for fo dh are @ 2018-10-01 22:30 and 2019-01-01 06:00
! First observation for dhdt is @ 2018-11-16 14:15 ( in the middle of the 1st and 2st dh obervation times)   
! Time interval is 91.3125 days or 7889400 seconds
! TODO: For PBS we need to compute the state VAR between the 1st and 2st dh obervation times

  yy = LIS_rc%yr
  mm = LIS_rc%mo
  dd = LIS_rc%da
  h  = LIS_rc%hr
  m  = LIS_rc%mn
  s  = 0 ! LIS_rc%sss
  ts=0
  call LIS_tick(timenow,doy,gmt,yy,mm,dd,h,m,s,real(ts))
  !! ICESat2 dhdt start date
  !call LIS_tick(start_date_tmp,doy,gmt,2018,11,16,12,0,0,0.0)
  !call LIS_tick(start_date,doy,gmt,2018,11,16,14,15,0,0.0)   ! dhdt start date 
  call LIS_tick(start_date,doy,gmt,2018,10,01,12,0,0,0.0)     ! dh   start date

  !call ESMF_TimeSet(startdate, yy=2018,&
  !     mm = 11, dd=16, h =14, &
  !     m = 15, s = 0, calendar = LIS_calendar,&
  !     rc=status)
  !call LIS_verify(status, 'ESMF_TimeSet failed in read_ATL15_GrISobs')

  !call ESMF_TimeSet(currentTime, yy=LIS_rc%yr,&
  !     mm = LIS_rc%mo, dd=LIS_rc%da, h =LIS_rc%hr, &
  !     m = LIS_rc%mn, s = 0, calendar = LIS_calendar,&
  !     rc=status)
  !call LIS_verify(status, 'ESMF_TimeSet failed in read_ATL15_GrISobs')

  !alarmcheck = (mod(currentTime-startdate, 7889400.0).eq.0)
  alarmcheck = (mod(timenow-start_date, 7889400.0).eq.0)

  !if(alarmCheck .and. currentTime.ge.start_date_tmp ) then
  if(alarmCheck .and. (timenow .ge. start_date)) then ! _tmp ) then
      !ATL15_GrIS_struc(n)%startMode = .false.

      data_upd = .false. 

      !if(LIS_rc%DAincrMode(n).eq.1) then  ! TODO what is this? hard codded to 1 in read_config
         
      call ESMF_AttributeSet(OBS_State,"File Status",&
            .true., rc=status)
      call LIS_verify(status)
           !call ESMF_AttributeSet(OBS_State,&
           !     name="Data averaging factor",&
           !     value=float(days(LIS_rc%mo)),rc=status)
           !call LIS_verify(status)

      write(LIS_logunit,*)  '[INFO] Reading ATL15_GrIS data ',trim(name)
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
      status = nf90_open(trim(name), nf90_NoWrite, ncid)
      if(status/=0) then
         if(LIS_masterproc) then
            write(LIS_logunit,*)'[ERR] Problem opening file: ',trim(filename),status
            call LIS_endrun
         endif
      else
         if(LIS_masterproc) then
            write(LIS_logunit,*)'[INFO] Opened file: ',trim(filename)
         endif
      endif

      call LIS_verify(nf90_inq_dimid(ncid, "dhdt_lag1/x",dimid),&
           'nf90_inq_dimid failed in dhdt_lag1/x, read_ATL15_GrIS')
      call LIS_verify(nf90_inquire_dimension(ncid, dimid, len=leng),&
           'nf90_inquire_dimension failed in dhdt_lag1/x, read_ATL15_GrIS')
      ATL15_GrIS_struc(n)%nc = leng
 
      call LIS_verify(nf90_inq_dimid(ncid, "dhdt_lag1/y",dimid),&
           'nf90_inq_dimid failed in dhdt_lag1/y, read_ATL15_GrIS')
      call LIS_verify(nf90_inquire_dimension(ncid, dimid, len=leng),&
           'nf90_inquire_dimension failed in dhdt_lag1/y, read_ATL15_GrIS')
      ATL15_GrIS_struc(n)%nr = leng

      call LIS_verify(nf90_inq_dimid(ncid, "dhdt_lag1/time",dimid),&
           'nf90_inq_dimid failed in dhdt_lag1/y, read_ATL15_GrIS')
      call LIS_verify(nf90_inquire_dimension(ncid, dimid, len=leng),&
           'nf90_inquire_dimension failed in dhdt_lag1/time, read_ATL15_GrIS')
      ATL15_GrIS_struc(n)%num_obs = leng

      !ATL15_GrIS_struc(n)%nc = 80
      !ATL15_GrIS_struc(n)%nr = 140
      !ATL15_GrIS_struc(n)%time = 20
      !allocate(X(nc))
      !allocate(Y(nr))
      allocate(obs_time(ATL15_GrIS_struc(n)%num_obs))
      allocate(dhdt(ATL15_GrIS_struc(n)%nc,ATL15_GrIS_struc(n)%nr,ATL15_GrIS_struc(n)%num_obs))
      allocate(dhdt_local_grid(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k),dhdt_local_grid))
 
      call LIS_verify(nf90_inq_varid(ncid,'dhdt_lag1/time',varid), &
           'nf90_inq_varid failed in dhdt_lag1/time, read_ATL15_GrIS')
      call LIS_verify(nf90_get_var(ncid,varid,obs_time, &
           'nf90_get_var failed in dhdt_lag1/time, read_ATL15_GrIS')

      call LIS_verify(nf90_inq_varid(ncid,'dhdt_lag1/dhdt',varid), &
          'nf90_inq_varid failed in dhdt_lag1/dhdt, read_ATL15_GrIS')
      call LIS_verify(nf90_get_var(ncid,varid,dhdt, &
          'nf90_get_var failed in dhdt_lag1/dhdt, read_ATL15_GrIS')
#endif

! TODO for PBS: How to store all obs in assimilation window? 
! For now let's implement PF 
! How to read in local domain?

      dhdt_local_grid = dhdt(&
             LIS_ews_obs_halo_ind(n,LIS_localPet+1):&
             LIS_ewe_obs_halo_ind(n,LIS_localPet+1), &
             LIS_nss_obs_halo_ind(n,LIS_localPet+1): &
             LIS_nse_obs_halo_ind(n,LIS_localPet+1),:)

!-------------------------------------------------------------------------
!  Extract data for the current time
!-------------------------------------------------------------------------     
      call ESMF_StateGet(OBS_State,"Observation01",dhdtfield,&
           rc=status)
      call LIS_verify(status, 'ESMF_StateGet failed in read_ATL15_GrISobs')
     
      call ESMF_FieldGet(dhdtfield,localDE=0,farrayPtr=obsl,rc=status)
      call LIS_verify(status,'ESMF_FieldGet failed in read_ATL15_GrISobs')
  
      obsl(:) = -9999.0

      !offset = int((currentTime - startTime)/ATL15_GrIS_struc(n)%ts) + 1 ! nint
      offset = floor((timenow - start_date)/7889400.0) + 1 ! nint
      ! we might be able to use obs_time insted of offset.   

      ! fnd: flag to indicate if there are valid observations in the local 
      !      processor's domain

      if(offset.gt.0) then
          fnd = 0 
          !data_upd_flag_local = .false. 
          do r=1,LIS_rc%obs_lnr(k)
             do c=1,LIS_rc%obs_lnc(k)
                if(dhdt_local_grid(c,r,offset).ne.LIS_rc%udef) then 
                   fnd = 1
                endif
             enddo
          enddo

          if(fnd.eq.0) then 
             obsl = LIS_rc%udef
          else
             do r=1,LIS_rc%obs_lnr(k)
                do c=1,LIS_rc%obs_lnc(k)
                   if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                      if(dhdt_local_grid(c,r,offset).gt.0.0) then
                         obsl(LIS_obs_domain(n,k)%gindex(c,r)) = dhdt_local_grid(c,r,offset)
                      end if
                   endif
                end do
             end do
          endif

          deallocate(obs_time)
          deallocate(dhdt)
          deallocate(dhdt_local_grid)

! TODO: do we need LSM based quality control and screening of observations
     !call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
     !     //trim(LIS_ATL15_GrISobsId)//char(0),n, k, OBS_state)
     !call LIS_checkForValidObs(n,k,obsl,fnd,dhdt_current)

          if(fnd.eq.0) then 
             data_upd_flag_local = .false.
          else
             data_upd_flag_local = .true. 
          endif
#if(defined SPMD)
          call MPI_ALLGATHER(data_upd_flag_local,1,&
             MPI_LOGICAL, data_upd_flag(:),&
             1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
          data_upd = .false.
          do p=1,LIS_npes
             data_upd = data_upd.or.data_upd_flag(p)
          enddo

          if(data_upd) then              
             do t=1,LIS_rc%obs_ngrid(k)
                gid(t) = t
                if(obsl(t).ne.-9999.0) then 
                   assimflag(t) = 1
                else
                   assimflag(t) = 0
                endif
             enddo

             call ESMF_AttributeSet(OBS_State,"Data Update Status",&
                  .true., rc=status)
             call LIS_verify(status)

             if(LIS_rc%obs_ngrid(k).gt.0) then 
                call ESMF_AttributeSet(dhdtfield,"Grid Number",&
                     gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(dhdtfield,"Assimilation Flag",&
                     assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)
             endif

          else
             call ESMF_AttributeSet(OBS_State,"Data Update Status",&
                 .false., rc=status)
             call LIS_verify(status)
          endif ! data_upd
      else
         call ESMF_AttributeSet(OBS_State,"Data Update Status",&
              .false., rc=status)
         call LIS_verify(status)
         return
      endif ! offset
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  endif !  alarmCheck
end subroutine read_ATL15_GrISobs

subroutine ATL15_GrIS_filename(name, ndir)
  
  implicit none
  character(len=*)  :: name
  character (len=*) :: ndir
! First and second observations for fo dh are @ 2018-10-01 22:30 and 2019-01-01 06:00
! First observation for dhdt is @ 2018-11-16 14:15 ( in the middle of the 1st and 2st dh obervation times)   
! Time interval is 91.3125 days or 7889400 seconds
! in the meta data it mentioned: height-change maps at 3 month intervals since 29 March 2019 are in a single file (BUT IT IS NOT CORRECT) 
  name = trim(ndir)//'/ATL15_GL_0311_20km_001_01.nc'
end subroutine ATL15_GrIS_filename



