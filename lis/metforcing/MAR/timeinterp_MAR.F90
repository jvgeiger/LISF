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
! !ROUTINE: timeinterp_MAR
! \label{timeinterp_MAR}
!
! !REVISION HISTORY:
!  22 Aug 2023: Mahdi Navari; Initial Specification 
!
! !INTERFACE:
subroutine timeinterp_MAR(n, findex)
! !USES:
  use ESMF
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc, LIS_FORC_Base_State
  use LIS_coreMod,       only : LIS_rc,LIS_domain
  use LIS_constantsMod,  only : LIS_CONST_SOLAR
  use LIS_timeMgrMod,    only : LIS_tick, LIS_time2date
  use LIS_logMod,        only : LIS_logunit, LIS_verify
  use MAR_forcingMod, only : MAR_struc
 
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation and longwave radiation
!  are not temporally interpolated, and the previous 1 hourly value
!  is used. All other variables are linearly interpolated between 
!  the 1 hourly blocks. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing source
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[LIS\_tick](\ref{LIS_tick}) \newline
!    advances or retracts time by the specified amount
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
  integer          :: zdoy
  real             :: zw1, zw2
  real             :: czm, cze, czb
  real             :: wt1, wt2,swt1,swt2
  real             :: gmt1, gmt2
  integer          :: t,index1
  integer          :: bdoy,byr,bmo,bda,bhr,bmn
  real*8           :: btime,newtime1,newtime2
  real             :: tempgmt1,tempgmt2
  integer          :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer          :: tempbss
  real             :: tempbts
  integer          :: status
  type(ESMF_Field) :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field) :: psurfField,rainField,snowField
  real,pointer     :: tair(:),qair(:),uwind(:),vwind(:)
  real,pointer     :: swd(:),lwd(:),psurf(:),rainf(:),snowf(:)
  !real,pointer     :: fheight(:),acond(:)

  !logical          :: forcing_z, forcing_ch
  
  btime=MAR_struc(n)%MARtime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  
  tempbdoy=bdoy
  tempgmt1=gmt1
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime1,tempbdoy,tempgmt1,& 
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn, & 
       tempbss,tempbts)
  
  btime=MAR_struc(n)%MARtime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  tempbdoy=bdoy
  tempgmt2=gmt2
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime2,tempbdoy,tempgmt2,&
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn,&
       tempbss,tempbts)
  
!=== Interpolate Data in time      
  wt1=(MAR_struc(n)%MARtime2-LIS_rc%time)/ & 
      (MAR_struc(n)%MARtime2-MAR_struc(n)%MARtime1)
  wt2=1.0-wt1
 
  ! Solar zenith angle weights for SWdown:
 swt1=(newtime2-LIS_rc%time)/(newtime2-newtime1)
  swt2=1.0-swt1


  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),q2Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),psurfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),rainField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Snowf%varname(1),snowField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Snowf in the forcing variables list')

#if 0
  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Forc_Hgt%varname(1),&
          fhgtField, rc=status)
     forcing_z = .true.
  else
     forcing_z = .false.
  endif

  if(LIS_FORC_Ch%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Ch%varname(1),acondField,&
          rc=status)
     forcing_ch = .true.
  else
     forcing_ch = .false.
  endif
#endif
    !     MAR       LIS 
    ! TT  (C)       K       tmp     <--> LIS_forc%metdata1(1,:)
    ! QQ  (g/kg)    kg/kg   q2      <--> LIS_forc%metdata1(2,:)
    ! SWD (W/m2)    W/m2    swd     <--> LIS_forc%metdata1(3,:)
    ! LWD (W/m2)    W/m2    lwd     <--> LIS_forc%metdata1(4,:)
    ! UU  (m/s)     m/s     uwind   <--> LIS_forc%metdata1(5,:)
    ! VV  (m/s)     m/s     vwind   <--> LIS_forc%metdata1(6,:)
    ! SP  (hpa)     pa      psurf   <--> LIS_forc%metdata1(7,:)
    ! RF  (mmWE/day)kg/m2s  (railfall)  <--> LIS_forc%metdata1(8,:)    
    ! SF  (mmWE/day)kg/m2s  (snowfall)  <--> LIS_forc%metdata1(9,:)


  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     zdoy=LIS_rc%doy
     !  compute and apply zenith angle weights
     call zterp(1,LIS_domain(n)%grid(index1)%lat, &
                LIS_domain(n)%grid(index1)%lon,   &
                gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)

     if(MAR_struc(n)%metdata1(3,index1).ne.LIS_rc%udef .and. &
        MAR_struc(n)%metdata2(3,index1).ne.LIS_rc%udef) then 
        swd(t) = MAR_struc(n)%metdata1(3,index1)*zw1+&
                 MAR_struc(n)%metdata2(3,index1)*zw2
!       In cases of small cos(zenith) angles, use linear weighting
!       to avoid overly large weights

        if((swd(t).gt.MAR_struc(n)%metdata1(3,index1)  .and. & 
            swd(t).gt.MAR_struc(n)%metdata2(3,index1)) .and. & 
           (czb.lt.0.1.or.cze.lt.0.1))then
           swd(t) = MAR_struc(n)%metdata1(3,index1)*swt1 + & 
                    MAR_struc(n)%metdata2(3,index1)*swt2
        endif
     endif
     
     if (swd(t).gt.LIS_CONST_SOLAR) then
        write(unit=LIS_logunit,fmt=*)'[WARN] SW radiation too high!!'
        write(unit=LIS_logunit,fmt=*)'[WARN] it is', swd(t)
        write(unit=LIS_logunit,fmt=*)'[WARN] data1=',MAR_struc(n)%metdata1(3,index1)
        write(unit=LIS_logunit,fmt=*)'[WARN] data2=',MAR_struc(n)%metdata2(3,index1)
        write(unit=LIS_logunit,fmt=*)'[WARN] zw1=',zw1,'zw2=',zw2
        write(unit=LIS_logunit,fmt=*)'[WARN] swt1=',swt1,'swt2=',swt2
 
!SWD is not set.. It is assumed to be filled from a valid
!baseforcing value.        
!        swd(t) = MAR_struc(n)%metdata1(3,index1)*swt1+ & 
!                 MAR_struc(n)%metdata2(3,index1)*swt2
!        swd(t) = LIS_rc%udef
!        write(unit=LIS_logunit,fmt=*)'forcing set to ',swd(t) 
     endif
  enddo

  !do block precipitation interpolation
  !WRF writes precipitation accumlated from the beginning of the
  !simulation.  Therefore compute the difference between the bookends

  !to find the accumulated precipitation for the hour.
#if 0
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if( (MAR_struc(n)%metdata1(8,index1).ne.LIS_rc%udef) .and. &
         (MAR_struc(n)%metdata2(8,index1).ne.LIS_rc%udef) ) then 
        pcp(t) = MAR_struc(n)%metdata2(8,index1) - &
                 MAR_struc(n)%metdata1(8,index1)
        pcp(t) = pcp(t)/(60.0*60.0)
        pcp(t) = max(pcp(t),0.0) ! added to eliminate round-off errors      
     endif
  enddo
#endif

  call ESMF_FieldGet(rainField,localDE=0,farrayPtr=rainf,rc=status)
  call LIS_verify(status)
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
        !rainf(t) = MAR_struc(n)%metdata2(8,index1) - MAR_struc(n)%metdata1(8,index1)
        rainf(t) = MAR_struc(n)%metdata1(8,index1)
        rainf(t) = rainf(t)/(60*60*24) ! (mmWE/day) --*1/(60*60*24)--> kg/m2s   ! mmWE/day  !/(60.0*60.0)
        rainf(t) = max(rainf(t),0.0) ! added to eliminate round-off errors      
  enddo          
        
  call ESMF_FieldGet(snowField,localDE=0,farrayPtr=snowf,rc=status)
  call LIS_verify(status)
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
        !snowf(t) = MAR_struc(n)%metdata2(9,index1) - MAR_struc(n)%metdata1(9,index1)
        snowf(t) = MAR_struc(n)%metdata1(9,index1)
        snowf(t) = snowf(t)/(60*60*24)! (mmWE/day) --*1/(60*60*24)--> kg/m2s   ! mmWE/day  !/(60.0*60.0)
        snowf(t) = max(snowf(t),0.0) ! added to eliminate round-off errors      
  enddo

  !linearly interpolate everything else

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tair,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((MAR_struc(n)%metdata1(1,index1).ne.LIS_rc%udef) .and. &
        (MAR_struc(n)%metdata2(1,index1).ne.LIS_rc%udef)) then 
        tair(t) = (MAR_struc(n)%metdata1(1,index1)*wt1 + & 
                   MAR_struc(n)%metdata2(1,index1)*wt2)+ 273.15  ! C --> K
     endif
  enddo

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=qair,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((MAR_struc(n)%metdata1(2,index1).ne.LIS_rc%udef) .and. &
        (MAR_struc(n)%metdata2(2,index1).ne.LIS_rc%udef)) then 
        qair(t) = (MAR_struc(n)%metdata1(2,index1)*wt1 + & 
                MAR_struc(n)%metdata2(2,index1)*wt2)/1000 !   g/kg --/1000-->   kg/kg
     endif
  enddo
     
  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((MAR_struc(n)%metdata1(4,index1).ne.LIS_rc%udef) .and. &
        (MAR_struc(n)%metdata2(4,index1).ne.LIS_rc%udef)) then 
        lwd(t) = MAR_struc(n)%metdata1(4,index1)*wt1 + & 
                 MAR_struc(n)%metdata2(4,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((MAR_struc(n)%metdata1(5,index1).ne.LIS_rc%udef) .and. &
        (MAR_struc(n)%metdata2(5,index1).ne.LIS_rc%udef)) then 
        uwind(t) = MAR_struc(n)%metdata1(5,index1)*wt1 + & 
                   MAR_struc(n)%metdata2(5,index1)*wt2
     endif
  enddo
  
  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((MAR_struc(n)%metdata1(6,index1).ne.LIS_rc%udef) .and. &
        (MAR_struc(n)%metdata2(6,index1).ne.LIS_rc%udef)) then 
        vwind(t) = MAR_struc(n)%metdata1(6,index1)*wt1 + & 
                   MAR_struc(n)%metdata2(6,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((MAR_struc(n)%metdata1(7,index1).ne.LIS_rc%udef) .and. &
        (MAR_struc(n)%metdata2(7,index1).ne.LIS_rc%udef)) then 
        psurf(t) = (MAR_struc(n)%metdata1(7,index1)*wt1 + & 
                   MAR_struc(n)%metdata2(7,index1)*wt2)*100 !  hPa/millibar --*100--> Pa
     endif
  enddo

#if 0
  if ( forcing_z ) then
     call ESMF_FieldGet(fhgtField,localDE=0,farrayPtr=fheight,rc=status)
     call LIS_verify(status)
  
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index 
        if((MAR_struc(n)%metdata1(10,index1).ne.LIS_rc%udef) .and. &
           (MAR_struc(n)%metdata2(10,index1).ne.LIS_rc%udef)) then 
           fheight(t) = MAR_struc(n)%metdata1(10,index1)*wt1 + & 
                        MAR_struc(n)%metdata2(10,index1)*wt2
        endif
     enddo
  endif

  if ( forcing_ch ) then
     call ESMF_FieldGet(acondField,localDE=0,farrayPtr=acond,rc=status)
     call LIS_verify(status)

     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index 
        if((MAR_struc(n)%metdata1(11,index1).ne.LIS_rc%udef) .and. &
           (MAR_struc(n)%metdata2(11,index1).ne.LIS_rc%udef)) then 
           acond(t) = MAR_struc(n)%metdata1(11,index1)*wt1 + & 
                      MAR_struc(n)%metdata2(11,index1)*wt2
        endif
     enddo
  endif
#endif

end subroutine timeinterp_MAR
 
