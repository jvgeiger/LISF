!------------------------------------------------------------------------------
!BOP

    module LIS_perform_regridMod


    ! ESMF Framework module
    use ESMF
    use ESMF_TestMod
    use LIS_coreMod,    only : LIS_rc, LIS_vm, LIS_masterproc
    use LIS_logMod,     only : LIS_verify, LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod
    use nldas2_forcingMod, only : nldas2_struc

    use LIS_timeMgrMod,     only : LIS_tick
    use LIS_metforcing_componentMod, only :   set_nldas2_file_name
    use LIS_metforcing_componentMod, only :   met_forcing_register
    use LIS_metforcing_componentMod, only : setFor_nest_index, setFor_forcast_member_index
    use LIS_regrid_couplerMod,       only : regridCoupler_register
    use LIS_model_componentMod,      only :         model_register
    use LIS_model_componentMod,      only : setMod_nest_index, setMod_forcast_member_index

    implicit none

    private
    public  ::        run_esmfRegrid
    public  ::   finalize_esmfRegrid
    public  :: initialize_esmfRegrid

    character(len=ESMF_MAXSTR) :: forcing_name
    character(len=ESMF_MAXSTR) :: model_name
    character(len=ESMF_MAXSTR) :: coupler_name
    type(ESMF_State)           :: forcing_export
    type(ESMF_State)           :: model_import
    type(ESMF_GridComp)        :: forcing_comp
    type(ESMF_GridComp)        :: model_comp
    type(ESMF_CplComp)         :: cpl
    type(ESMF_Clock)           :: clock
    integer                    :: rc, localrc, userrc
    integer                    :: new_day

    integer, save :: findtime1, findtime2
    integer, save :: rstflag(3) = 1
    real,    save :: nldas2time1, nldas2time2
!EOP
!------------------------------------------------------------------------------
    CONTAINS
!------------------------------------------------------------------------------
!BOP

      SUBROUTINE initialize_esmfRegrid(n)
!
      integer, intent(in) :: n
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (LIS_masterproc) print *, "Start of regrid_test:"

      !------------------------------------------
      ! Create the 2 model components and coupler
      !------------------------------------------
      forcing_name = "Forcing Data"
      forcing_comp = ESMF_GridCompCreate(name=forcing_name, rc=localrc)
      call LIS_verify(localrc, 'ESMF_GridCompCreate failed for forcing data')
      if (LIS_masterproc) print *, "---> Created component ", trim(forcing_name), "rc =", rc

      model_name = "Model Data"
      model_comp = ESMF_GridCompCreate(name=model_name, rc=localrc)
      call LIS_verify(localrc, 'ESMF_GridCompCreate failed for model data')
      if (LIS_masterproc) print *, "---> Created component ", trim(model_name), "rc =", rc

      coupler_name = "user one-way coupler"
      cpl = ESMF_CplCompCreate(name=coupler_name, rc=localrc)
      call LIS_verify(localrc, 'ESMF_CplCompCreate failed')
      if (LIS_masterproc) print *, "---> Created component ", trim(coupler_name), ", rc =", rc

      if (LIS_masterproc) print *, "---> Comp Creates finished"

      !-------------------
      !  Register section
      !-------------------
      call ESMF_GridCompSetServices(forcing_comp, met_forcing_register, userRc=userrc, rc=localrc)
      call LIS_verify(localrc, 'ESMF_GridCompSetServices failed for forcing data')
      if (LIS_masterproc) print *, "---> Forcing Comp SetServices finished, rc= ", rc, userrc

      call ESMF_GridCompSetServices(model_comp, model_register, userRc=userrc, rc=localrc)
      call LIS_verify(localrc, 'ESMF_GridCompSetServices failed for model data')
      if (LIS_masterproc) print *, "---> Model Comp SetServices finished, rc= ", rc, userrc

      call ESMF_CplCompSetServices(cpl, regridCoupler_register, userRc=userrc, rc=localrc)
      call LIS_verify(localrc, 'ESMF_CplCompSetServices failed')
      if (LIS_masterproc) print *, "---> Coupler Comp SetServices finished, rc= ", rc, userrc

      !---------------
      !  Init section
      !---------------

      forcing_export = ESMF_StateCreate(name="forcing_comp export",  &
                               stateintent=ESMF_STATEINTENT_EXPORT, rc=localrc)
      call LIS_verify(localrc, ' ESMF_StateCreate failed')

      call ESMF_GridCompInitialize(forcing_comp, exportState=forcing_export, clock=clock, &
                                                 userRc=userrc, rc=localrc)
      call LIS_verify(localrc, 'ESMF_GridCompInitialize failed')
      if (LIS_masterproc) print *, "Forcing Comp Initialize finished, rc =", rc

      model_import = ESMF_StateCreate(name="model_comp import",  &
                               stateintent=ESMF_STATEINTENT_IMPORT, rc=localrc)
      call LIS_verify(localrc, 'ESMF_StateCreat failed')

      call ESMF_GridCompInitialize(model_comp, importState=model_import, clock=clock, &
                                               userRc=userrc, rc=localrc)
      call LIS_verify(localrc, 'ESMF_GridCompInitialize failed')
      if (LIS_masterproc) print *, "Model Comp Initialize finished, rc =", rc

      ! note that the coupler's import is forcing_comp's export
      call ESMF_CplCompInitialize(cpl, importState=forcing_export, &
                                       exportState=model_import, clock=clock, &
                                       userRc=userrc, rc=localrc)
      call LIS_verify(localrc, 'ESMF_CplCompInitialize failed')
      if (LIS_masterproc) print *, "Coupler Initialize finished, rc =", rc

      nldas2time1 = nldas2_struc(n)%nldas2time1
      nldas2time2 = nldas2_struc(n)%nldas2time2

      new_day = -999

      END SUBROUTINE initialize_esmfRegrid
!EOC
!------------------------------------------------------------------------------
!BOP

      SUBROUTINE run_esmfRegrid(n, findex)
!EOP
         integer, intent(in) :: n
         integer, intent(in) :: findex

         integer :: c,f,ferrora,ferrorb,ferror,try
         integer :: order
         integer :: readbfile
         real*8  :: time1,time2,timenow
         real*8  :: dtime1, dtime2
         integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
         integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
         character*100 :: name_a,name_b
         real    :: gmt1,gmt2,ts1,ts2
         integer :: movetime     ! 1=move time 2 data into time 1  
         integer :: kk           ! Forecast member index

!------------------------------------------------------------------------------
!BOC

!====Assumption will be not to find or move any data
  nldas2_struc(n)%findtime1=0
  nldas2_struc(n)%findtime2=0
  movetime=0

!=== Determine Required NLDAS-2 Data Times (The previous hour and the future hour)
  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  if(LIS_rc%ts.gt.3600) then
     write(LIS_logunit,*) '[ERR] The model timestep is > forcing data timestep'
     write(LIS_logunit,*) '[ERR] LIS does not support this mode currently'
     write(LIS_logunit,*) '[ERR] Program stopping ...'
     call LIS_endrun()
  endif

  if(mod(nint(LIS_rc%ts),3600).eq.0) then
     !if(timenow.ge.nldas2time2) then
     if(timenow.ge.nldas2_struc(n)%nldas2time2) then
        yr1 = LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        hr1=LIS_rc%hr
        mn1=0
        ss1=0
        ts1=-60*60
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=LIS_rc%hr
        mn2=0
        ss2=0
        ts2=0
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        nldas2_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.nldas2_struc(n)%nldas2time2) then
        yr1 = LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        hr1=LIS_rc%hr
        mn1=0
        ss1=0
        ts1=0
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=LIS_rc%hr
        mn2=0
        ss2=0
        ts2=60*60
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

        movetime = 1
        nldas2_struc(n)%findtime2 = 1
     endif
  endif

  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then    !beginning of the run 
     nldas2_struc(n)%findtime1=1
     nldas2_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif

  if(movetime.eq.1) then
     nldas2_struc(n)%nldas2time1=nldas2_struc(n)%nldas2time2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           nldas2_struc(n)%metdata1(:,f,c)=nldas2_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif    !end of movetime=1

  if(nldas2_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferrora=0
     ferrorb=0
     ferror=0
     try=0
     ts1=-60*60*24
     do
        if ( ferror /= 0 ) exit
        try=try+1
        !- Obtaining NLDAS-2 File-A:
        do kk= nldas2_struc(n)%st_iterid, nldas2_struc(n)%en_iterid
           CALL run_esmfRegrid_oneStep('A', n, kk, yr1, mo1, da1, &
                                      doy1, hr1)
        enddo

        if( readbfile .gt. 0) then
          do kk= nldas2_struc(n)%st_iterid, nldas2_struc(n)%en_iterid
             CALL run_esmfRegrid_oneStep('B', n, kk, yr1, mo1, da1, &
                                      doy1, hr1)
          enddo
        else
           ferrorb = 1
        endif
        ferror = ferrora + ferrorb
        if(ferror.ge.1) nldas2_struc(n)%nldas2time1=time1
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(*,*)'error: NLDAS-2 data gap exceeds 10 days on file 1'
           stop
        endif
     enddo
!=== end of data search
  endif   !end of LIS_rc%findtime=1 

  if(nldas2_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferrora=0
     ferrorb=0
     ferror=0
     try=0
     ts2=-60*60*24
     do
        if ( ferror /= 0 ) exit
        try=try+1

     !- Obtaining NLDAS-2 File-A:
        do kk= nldas2_struc(n)%st_iterid, nldas2_struc(n)%en_iterid
           CALL run_esmfRegrid_oneStep('A', n, kk, yr2, mo2, da2, &
                                      doy2, hr2)
        end do

        if( readbfile .gt. 0) then
          do kk= nldas2_struc(n)%st_iterid, nldas2_struc(n)%en_iterid
             CALL run_esmfRegrid_oneStep('B', n, kk, yr2, mo2, da2, &
                                      doy2, hr2)
          enddo
        else
           ferrorb = 1
        endif
        ferror = ferrora + ferrorb
        if(ferror.ge.1) then
           nldas2_struc(n)%nldas2time2=time2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(*,*)'error: NLDAS-2 data gap exceeds 10 days on file 2'
           stop
        endif
     enddo
     !=== end of data search
  endif   ! end of findtime2=1


      CONTAINS
         subroutine run_esmfRegrid_oneStep(file_type, n, kk, curYear, curMonth, curDay, &
                                      curDayOfYear, curHour)

            character(len=1) :: file_type
            integer :: n, kk, curYear, curMonth, curDay, curDayOfYear, curHour

            CALL set_nldas2_file_name(file_type, n, curYear, curMonth, curDay, &
                                      curDayOfYear, curHour)

            CALL setFor_forcast_member_index(kk)
            CALL setFor_nest_index(n)
  
            call ESMF_GridCompRun(forcing_comp, exportState=forcing_export, clock=clock, &
                                                 userRc=userrc, rc=localrc)
            call LIS_verify(localrc, 'ESMF_GridCompRun failed')

            call ESMF_CplCompRun(cpl, importState=forcing_export, &
                                      exportState=model_import, clock=clock, &
                                      userRc=userrc, rc=localrc)
            call LIS_verify(localrc, 'ESMF_CplCompRun failed')

            CALL setMod_forcast_member_index(kk)
            CALL setMod_nest_index(n)

            call ESMF_GridCompRun(model_comp, importState=model_import, clock=clock, &
                                           userRc=userrc, rc=localrc)
            call LIS_verify(localrc, 'ESMF_GridCompRun failed')
         end subroutine run_esmfRegrid_oneStep

      END SUBROUTINE run_esmfRegrid
!EOC
!------------------------------------------------------------------------------
!BOP

      SUBROUTINE finalize_esmfRegrid()

    call ESMF_GridCompFinalize(forcing_comp, exportState=forcing_export, clock=clock, &
                                             userRc=userrc, rc=localrc)
    call LIS_verify(localrc, 'ESMF_GridCompFinalize failed')
    if (LIS_masterproc) print *, "---> Forcing Comp Finalize finished, rc =", rc

    call ESMF_GridCompFinalize(model_comp, importState=model_import, clock=clock, &
                                           userRc=userrc, rc=localrc)
    call LIS_verify(localrc, 'ESMF_GridCompFinalize failed')
    if (LIS_masterproc) print *, "---> Model Comp Finalize finished, rc =", rc

    call ESMF_CplCompFinalize(cpl, importState=forcing_export, &
                                   exportState=model_import, clock=clock, &
                                   userRc=userrc, rc=localrc)
    call LIS_verify(localrc, 'ESMF_CplCompFinalize failed')
    if (LIS_masterproc) print *, "---> Coupler Finalize finished, rc =", rc

    if (LIS_masterproc) print *, "Comp Finalize returned"

!
    call ESMF_StateDestroy(forcing_export, rc=rc)
    call LIS_verify(localrc, 'ESMF_StateDestroy failed')

    call ESMF_StateDestroy(model_import, rc=rc)
    call LIS_verify(localrc, 'ESMF_StateDestroy failed')

    call ESMF_ClockDestroy(clock, rc=rc)
    call LIS_verify(localrc, 'ESMF_ClockDestroy failed')

    call ESMF_GridCompDestroy(forcing_comp, rc=rc)
    call LIS_verify(localrc, 'ESMF_GridCompDestroy failed')

    call ESMF_GridCompDestroy(model_comp, rc=rc)
    call LIS_verify(localrc, 'ESMF_GridCompDestroy failed')

    call ESMF_CplCompDestroy(cpl, rc=rc)
    call LIS_verify(localrc, 'ESMF_CplCompDestroy failed')

    if (LIS_masterproc) print *, "All Destroy routines done"

      END SUBROUTINE finalize_esmfRegrid
!EOC
!-------------------------------------------------------------------------
    end module LIS_perform_regridMod

