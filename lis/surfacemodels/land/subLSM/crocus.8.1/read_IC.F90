!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

!BOP
!
! !ROUTINE: read_IC
! \label{read_IC}
!
! !DESCRIPTION:
!  The code in this file read the Ice sheet initial profile from the MAR forcing data
!
! !REVISION HISTORY:
!  8 Nov 2023: Mahdi Navari; Initial version
!
subroutine read_IC(n)

! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_timeMgrMod, only : LIS_tick
  use LIS_logMod
  use LIS_coreMod, only      : LIS_rc, LIS_masterproc, LIS_localPet, LIS_config
  use Crocus81_lsmMod

  implicit none
! !INPUT PARAMETERS:
!
  integer, intent(in) :: n
  character(len=LIS_CONST_PATH_LEN) :: filename, MARdir
  integer           :: yy,mm,dd,h,m,s
  integer           :: doy, ts
  integer           :: rc
  real              :: gmt
  real*8            :: timenow
  integer           :: timestep
  character*4       :: cyr
  character*2       :: cmo
  logical           :: file_exists
  external          :: read_IC_from_MAR

  if((LIS_rc%startcode .eq. "coldstart") .and.  CROCUS81_struc(n)%Init_Profile_from_MAR_BOOL) then

  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour need to be read.
  !----------------------------------------------------------------
     yy = LIS_rc%syr
     mm = LIS_rc%smo
     dd = LIS_rc%sda
     h  = LIS_rc%shr
     m  = 0 ! LIS_rc%smn
     s  = 0 ! LIS_rc%sss

     ts=0
     call LIS_tick(timenow,doy,gmt,yy,mm,dd,h,m,s,real(ts))

    ! One file per month-- get to the time record
    ! start hr in monthly and yearly file is 1 (not 0)
    ! hr 0-23;  
    timestep = 24*(dd - 1) + h  ! for forcing we use timestep = 24*(da - 1) + hr
                                ! e.g., start hr = 5 LIS forcing reader reads timestep 5  and 6
                                ! read_IC reads inital condition from timestep 5 
    if(LIS_masterproc) then
       write(LIS_logunit,*)'[INFO] timestep, LIS day, LIS hr ::', timestep, dd, h
    endif

    call ESMF_ConfigFindLabel(LIS_config,"MAR forcing directory:",rc=rc)
    call LIS_verify(rc, 'MAR forcing directory: not defined')
    call ESMF_ConfigGetAttribute(LIS_config,MARdir,rc=rc)

    write(unit=cyr,fmt='(i4.4)') yy
    write(unit=cmo,fmt='(i2.2)') mm
    filename = trim(MARdir)//'/'//&
               'ICE.'//trim(cyr)//trim(cmo)//'.b85.nc'

!     'ICE.'//trim(cyr)//'.01-12.b85.nc'
!     'MARv3.12-ERA5-25km-hourly-'//trim(cyr)//'.nc'
!      ICE.202207.b85.nc 

    inquire(file=filename, exist=file_exists)
    if(.not. file_exists)then
       if(LIS_masterproc) then
          write(LIS_logunit,*) "Restart file (MAR forcing) ", filename," does not exist "
          write(LIS_logunit,*) "Program stopping ..."
          call LIS_endrun
       endif 
    else
       if(LIS_masterproc) then
          write(LIS_logunit,*)'[INFO] Reading inital profile from'
          write(LIS_logunit,*)'[INFO] MAR forcing data .. ',trim(filename)
       endif   
       call read_IC_from_MAR(n,filename, timestep)
    endif
  endif ! coldstart and Init_Profile_from_MAR_BOOL
end subroutine read_IC

!BOP
!
! !ROUTINE: read_IC_from_MAR
! \label{read_IC_from_MAR}

subroutine read_IC_from_MAR(n,filename, timestep) ! TODO 

! !USES:
    use ESMF
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use LIS_coreMod, only    : LIS_rc, LIS_domain, LIS_masterproc,&
                               LIS_localPet, LIS_config, &
                               LIS_ews_halo_ind, LIS_ewe_halo_ind,&
                               LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify
    use LIS_timeMgrMod, only : LIS_tick
    use Crocus81_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

! !DESCRIPTION:
!  This program reads restart files for Crocus81 from MAR forcing. This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the Crocus81
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    SNOWSWE                    - Crocus81 Snow layer(s) liquid Water Equivalent (SWE:kg m-2) [kg/m2]
!    SNOWRHO                    - Crocus81 Snow layer(s) averaged density (kg/m3) [kg/m3]
!    SNOWHEAT                   - Crocus81 Snow layer(s) heat content (J/m2) [J/m2]
!    SNOWALB                    - Crocus81 snow surface albedo [-]
!    SNOWGRAN1                  - Crocus81 Snow layers grain feature 1 [-]
!    SNOWGRAN2                  - Crocus81 Snow layer grain feature 2 [-]
!    SNOWHIST                   - Crocus81 Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5} [-]
!    SNOWAGE                    - Crocus81 Age since snowfall (day) [day]
!    SNOWLIQ                    - Crocus81 Snow layer(s) liquid water content (m) [m]
!    SNOWTEMP                   - Crocus81 Snow layer(s) temperature (K) [K]
!    SNOWDZ                     - Crocus81 Snow layer(s) thickness (m) [m]
!    GRNDFLUX                   - Crocus81 Soil/snow interface heat flux (W/m2) [W/m2]
!    SNDRIFT                    - Crocus81 Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592) [kg/m2/s]
!    RI_n                       - Crocus81 Richardson number (-)  NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90) [-]
!    CDSNOW                     - Crocus81 Drag coefficient for momentum over snow (-) [-]
!    USTARSNOW                  - Crocus81 Friction velocity over snow (m/s); [m/s]
!    CHSNOW                     - Crocus81 Drag coefficient for heat over snow  (-) [-]
!    SNOWMAK_dz                 - Crocus81 Snowmaking thickness (m) [m]
!  \end{verbatim}
!
!EOP

  implicit none

  integer           :: v, ly, t, l
  integer           :: nc, nr, SNOLAY
  integer           :: n,rc
  integer           :: ftn
  integer           :: status
  integer           :: ios              ! set to non-zero if there's an error
  character (len=LIS_CONST_PATH_LEN)        :: filename !, MARdir
  integer           :: yy,mm,dd,h,m,s
  integer           :: doy, ts
  real              :: gmt
  real*8            :: timenow
  integer           :: timestep
  character*4       :: cyr
  character*2       :: cmo
  integer           :: dimid, AlbId, ncid, varid, N_active_layer
  integer           :: leng
  ! integer           :: c, j, i, findex, iret
  real, allocatable :: tmpAlb(:,:), Alb2Cro(:,:), tmpZSCAP(:,:,:)
  integer           :: ncol,nrow,nlayer   ! MAR bottom layer is layer=1 and top layer is layer=21
  REAL              :: XRHOLW             ! Volumic mass of liquid water
  REAL              :: XCL,XCI            ! Cl (liquid), Ci (ice)
  REAL              :: XTT                ! Triple point temperature
  REAL              :: XLVTT              ! Vaporization heat constant
  REAL              :: XLSTT              ! Sublimation heat constant
  REAL              :: XLMTT              ! Melting heat constant
  integer, parameter:: N_IC = 11  ! # of IC variables
  character(8), dimension(N_IC), parameter :: MAR_ICv = (/  &
         'DZSN1   ',&
         'ROSN1   ',&
         'TISN1   ',&
         'WASN1   ',&
         'G1SN1   ',&
         'G2SN1   ',&
         'AGSN1   ',&
         'NHSN1   ',&
         'LWC     ',&
         'SWE     ',&
         'SnowHeat'/)
  !       SD ,Rho ,SnowT ,SnowHum, G1, G2, Age ,History, LWC, SWE, SnowHeat

  real,allocatable  :: datain(:,:,:)      ! input data 
  real,allocatable  :: tmp2Cro(:,:,:,:)
  real,allocatable  :: data_in_local_grid(:,:,:)
  real,parameter    :: missing_value = -1.e+34  ! LIS_rc%udef
!print*,'read_IC.F90 LIS_localPet', LIS_localPet

  XRHOLW = 1000.
  XCL    = 4.218E+3
  XCI    = 2.106E+3
  XTT    = 273.16
  XLVTT  = 2.5008E+6
  XLSTT  = 2.8345E+6
  XLMTT  = XLSTT - XLVTT

!-----------------------------------------------------------------------    
! Extract dim form MAR forcing files 
!----------------------------------------------------------------------- 
! Open netCDF file.
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  status = nf90_open(trim(filename), nf90_NoWrite, ncid)
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

  call LIS_verify(nf90_inq_dimid(ncid, "X10_85",dimid),&
       'nf90_inq_dimid failed in X10_85, read_IC_from_MAR')
  call LIS_verify(nf90_inquire_dimension(ncid, dimid, len=leng),&
       'nf90_inquire_dimension failed in X10_85, Crocus81_IC_Mod')
  nc = leng

  call LIS_verify(nf90_inq_dimid(ncid, "Y20_155",dimid),&
       'nf90_inq_dimid failed in Y20_155,Crocus81_IC_Mod')
  call LIS_verify(nf90_inquire_dimension(ncid, dimid, len=leng),&
       'nf90_inquire_dimension failed in Y20_155, read_IC_from_MAR')
  nr = leng

  call LIS_verify(nf90_inq_dimid(ncid, 'SNOLAY',dimid),&
       'nf90_inq_dimid failed in SNOLAY, Crocus81_IC_Mod')
  call LIS_verify(nf90_inquire_dimension(ncid, dimid, len=leng),&
       'nf90_inquire_dimension failed in SNOLAY, Crocus81_IC_Mod')
  SNOLAY = leng
! Close netCDF file.
  call LIS_verify(nf90_close(ncid))
#endif

!-----------------------------------------------------------------------    
! Allocate memory
!-----------------------------------------------------------------------
  allocate(tmpAlb(nc,nr))
  allocate(Alb2Cro(LIS_rc%lnc(n),LIS_rc%lnr(n)))
  allocate(datain(nc,nr,SNOLAY))
  allocate(data_in_local_grid(LIS_rc%lnc(n),LIS_rc%lnr(n),SNOLAY))
  allocate(tmp2Cro(LIS_rc%lnc(n),LIS_rc%lnr(n),SNOLAY,N_IC), stat=ios)
  if(ios.ne.0) then
    write(LIS_logunit,*) '[ERR] Error allocating tmp2Cro,',LIS_localPet
    call LIS_endrun
  endif

!=== Open MAR forcing files ===
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  status = nf90_open(trim(filename), nf90_NoWrite, ncid)
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

!-----------------------------------------------------------------------    
! Read albedo
!-----------------------------------------------------------------------
  call LIS_verify(nf90_inq_varid(ncid,'AL2',AlbId), &
       'nf90_inq_varid failed for Albedo in MAR')
  call LIS_verify(nf90_get_var(ncid,AlbId, tmpAlb, &
       start=(/1,1,1,timestep/), &
       count=(/nc,nr,1,1/)),&
       'nf90_get_var failed for Albedo in MAR')
  Alb2Cro(:,:) = tmpAlb(LIS_ews_halo_ind(n,LIS_localPet+1):&
                LIS_ewe_halo_ind(n,LIS_localPet+1), &
                LIS_nss_halo_ind(n,LIS_localPet+1): &
                LIS_nse_halo_ind(n,LIS_localPet+1))

!-----------------------------------------------------------------------    
! Read SD,Rho,SnowT,SnowHum,G1,G2,Age,and Hist from MAR forcing
!-----------------------------------------------------------------------
  do v = 1,8 
     call LIS_verify(nf90_inq_varid(ncid,trim(MAR_ICv(v)),varId), &
          'nf90_inq_varid failed for '//trim(MAR_ICv(v))//' in MAR')
     call LIS_verify(nf90_get_var(ncid,varId, datain, &
          start=(/1,1,1,timestep/), &
          count=(/nc,nr,SNOLAY,1/)),&
          'nf90_get_var failed for getting IC from MAR')

      data_in_local_grid(:,:,:) = datain(LIS_ews_halo_ind(n,LIS_localPet+1):&
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1),:)
      tmp2Cro(:,:,:,v) = data_in_local_grid(:,:,:)
 end do !v loop
 ! Close netCDF file.
 call LIS_verify(nf90_close(ncid))
#endif

!-----------------------------------------------------------------------    
! Calculate  LWC, SWE, SnowHeat content
!-----------------------------------------------------------------------
! Crocus needs LWC to be in [m] 
! 1) liquid water content can be deduced from WASN1  LWC (in %) = 0.1*ROSN1*WASN1 
! Xavier --> The factor 0.1 is in fact 100 (for being in %) /1000 (density of water). 
! [100 * kg_w/kg_ice * kg_ice/m3 / 1000 (m3/kg_w) = 0.1*ROSN1*WASN1  [%]  IGNOR  THIS WE NEED LWC in meter
! 2) LWC (in m) = ROSN1*WASN1*DZSN1  
!   [kg_snow/m3] * [kg_w/kg_snow] * [m] = [kg_w/m2] which is mm --> /1000 --> m
    
! '1.DZSN1(SD)', '2.ROSN1(Rho)', '3.TISN1(SnowT)', '4.WASN1(SnowHum)', '5.G1SN1', '6.G2SN1', '7.AGSN1(Age)', '8.NHSN1(His)'
  do nrow=1,LIS_rc%lnr(n)
     do ncol=1,LIS_rc%lnc(n)
        do nlayer=1,SNOLAY
           if (tmp2Cro(ncol,nrow,nlayer,1) .gt. 0) then
              tmp2Cro(ncol,nrow,nlayer,9) = &
              tmp2Cro(ncol,nrow,nlayer,2) * &
              tmp2Cro(ncol,nrow,nlayer,4) * &
              tmp2Cro(ncol,nrow,nlayer,1) / 1000 ! [m]
           endif
        enddo
     enddo
  enddo

! SWE can be computed using ROSN1*DZSN1
  do nrow=1,LIS_rc%lnr(n)
     do ncol=1,LIS_rc%lnc(n)
        do nlayer=1,SNOLAY
           if (tmp2Cro(ncol,nrow,nlayer,1) .gt. 0) then
              tmp2Cro(ncol,nrow,nlayer,10) = &
              tmp2Cro(ncol,nrow,nlayer,2)* &
              tmp2Cro(ncol,nrow,nlayer,1)
           endif
        enddo
     enddo
  enddo

! snow heat content can be computed using Crocus equations
! 1) Update snow heat content (J/m2) using dry density instead of total density:
!     ZSCAP(JJ,JST) = ( PSNOWRHO(JJ,JST) - &
!                      PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)) * XCI
!    PSNOWHEAT(JJ,JST) = PSNOWDZ(JJ,JST) * &
!                        ( ZSCAP(JJ,JST)*(ZSNOWTEMP(JJ,JST)-XTT) - XLMTT*PSNOWRHO(JJ,JST) ) + &
!                        XLMTT * XRHOLW * PSNOWLIQ(JJ,JST)
! 2) Heat content using total density
!    ZSCAP     (JJ,JST) = PSNOWRHO(JJ,JST) * XCI
!    PSNOWHEAT (JJ,JST) = PSNOWDZ(JJ,JST) * &
!                        ( ZSCAP(JJ,JST)*(PSNOWTEMP(JJ,JST)-XTT) - XLMTT*PSNOWRHO(JJ,JST) ) + &
!                        XLMTT * XRHOLW * PSNOWLIQ(JJ,JST)

  allocate(tmpZSCAP(LIS_rc%lnc(n),LIS_rc%lnr(n),SNOLAY))
  tmpZSCAP = LIS_rc%udef
  do nrow=1,LIS_rc%lnr(n)
     do ncol=1,LIS_rc%lnc(n)
        do nlayer=1,SNOLAY
           if (tmp2Cro(ncol,nrow,nlayer,1) .gt. 0) then
              ! convert snow temp to K 
              tmp2Cro(ncol,nrow,nlayer,3) = tmp2Cro(ncol,nrow,nlayer,3) + 273.16
              ! TODO check the value. It should be nagative.
              tmpZSCAP(ncol,nrow,nlayer)   = (tmp2Cro(ncol,nrow,nlayer,2) - &
                                             tmp2Cro(ncol,nrow,nlayer,9) * & 
                                             XRHOLW / tmp2Cro(ncol,nrow,nlayer,1)) * XCI
              tmp2Cro(ncol,nrow,nlayer,11) = tmp2Cro(ncol,nrow,nlayer,1) * &
                                             (tmpZSCAP(ncol,nrow,nlayer) * & 
                                             (tmp2Cro(ncol,nrow,nlayer,3) - XTT) - &
                                             XLMTT * tmp2Cro(ncol,nrow,nlayer,2)) + &
                                             XLMTT * XRHOLW * tmp2Cro(ncol,nrow,nlayer,9)
           endif
        enddo
     enddo
  enddo


!-----------------------------------------------------------------------    
! Set IC to Crocus IC structure
!----------------------------------------------------------------------- 
! Note: 
! In MAR bottom layer is layer=1 and top layer is layer=21
! Crocus expect top layer to be layer=1 

! '1.DZSN1(SD)', '2.ROSN1(Rho)', '3.TISN1(SnowT)', '4.WASN1(SnowHum)', '5.G1SN1', '6.G2SN2', '7.AGSN1(Age)', '8.NHSN1(His)'
!  9.LWC          10. SWE         11. SnowHeat

  do nrow=1,LIS_rc%lnr(n)
     do ncol=1,LIS_rc%lnc(n)
       !if (ncol .eq. 35 .and. nrow .eq. 9)then
       !print*, 'ncol, nrow', ncol, nrow
       !endif 
        t = LIS_domain(n)%gindex(ncol,nrow)! there is one IC for each grid 
        !do l = 1,SNOLAY  
        !   if (tmp2Cro(ncol,nrow,l,1) .LE. 0.000001) then
        !      tmp2Cro(ncol,nrow,l,2) = 10000000000
        !   endif
        !enddo 
        if (LIS_domain(n)%gindex(ncol,nrow).ne.-1) then
           N_active_layer = 0 !SNOLAY 
           do l = 1,SNOLAY  ! from bottom to top  
              if (tmp2Cro(ncol,nrow,l,1) .GE. 0.000001) then ! XSNOWDMIN = 0.000001  ! (m) to prevent numerical problems as snow becomes vanishingly thin.
                 !if (tmp2Cro(ncol,nrow,l,2) .eq. 0) then
                 !tmp2Cro(ncol,nrow,l,2) = 10000000000
                 !endif
                 N_active_layer = l
              endif
           enddo
              CROCUS81_struc(n)%crocus81(t)%SNOWRHO(:) = 10000000000
              if (N_active_layer .gt. 0) then 
              do ly=1, N_active_layer !Crocus81_IC_struc(n)%SNOLAY 
                 CROCUS81_struc(n)%crocus81(t)%SNOWDZ(ly)   = tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,1)!Snow layer(s) thickness (m)
                 CROCUS81_struc(n)%crocus81(t)%SNOWRHO(ly)  = tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,2)!Snow layer(s) averaged density (kg/m3)
!if (LIS_localPet .eq. 1) then
!  print*,'t, LPet, SNOWRHO', t, LIS_localPet, tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,2)            
!endif
                 CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(ly) = tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,3)!Snow layer(s) temperature (K)
                 CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(ly)= tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,5)!Snow layers grain feature 1 
                 CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(ly)= tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,6)!Snow layers grain feature 2 
                 CROCUS81_struc(n)%crocus81(t)%SNOWAGE(ly)  = tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,7)!Age since snowfall (day) 
                 CROCUS81_struc(n)%crocus81(t)%SNOWHIST(ly) = tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,8)!Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5}
                 CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(ly)  = tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,9)!Snow layer(s) liquid water content (m)
                 CROCUS81_struc(n)%crocus81(t)%SNOWSWE(ly)  = tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,10)!Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
                 CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(ly) = tmp2Cro(ncol,nrow,(N_active_layer+1)-ly,11)!Snow layer(s) heat content (J/m2)
              enddo ! ly
              else
              CROCUS81_struc(n)%crocus81(t)%SNOWDZ(:)   = 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWRHO(:)  = 10000000000
              CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(:) = 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(:)= 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(:)= 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(:)= 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWAGE(:)  = 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWAGE(:)  = 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWHIST(:) = 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(:)  = 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWSWE(:)  = 0.
              CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(:) = 0.
              endif 
                 CROCUS81_struc(n)%crocus81(t)%SNOWALB = Alb2Cro(ncol,nrow)!snow surface albedo

                 !TODO : check the following values 
                 !read: Soil/snow interface heat flux (W/m2)
                 CROCUS81_struc(n)%crocus81(t)%GRNDFLUX = 0.0

                 !read: Blowing snow sublimation (kg/m2/s) 
                 !NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592)
                 CROCUS81_struc(n)%crocus81(t)%SNDRIFT = 0.0

                 !read: Richardson number (-)  
                 !NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90)
                 CROCUS81_struc(n)%crocus81(t)%RI_n = 0.01 ! @2m for u=3 to 5 m/s % https://doi.org/10.3189/S0022143000034882

                 !read: Drag coefficient for momentum over snow (-)
                 CROCUS81_struc(n)%crocus81(t)%CDSNOW = 0

                 !read: Friction velocity over snow (m/s);
                 CROCUS81_struc(n)%crocus81(t)%USTARSNOW = 0

                 !read: Drag coefficient for heat over snow  (-)
                 CROCUS81_struc(n)%crocus81(t)%CHSNOW = 0

                 !read: Snowmaking thickness (m)
                 CROCUS81_struc(n)%crocus81(t)%SNOWMAK_dz = 0
        endif
     enddo ! c
  enddo ! r

     deallocate(tmpAlb)
     deallocate(Alb2Cro)
     deallocate(tmpZSCAP)
     deallocate(datain)
     deallocate(data_in_local_grid)
     deallocate(tmp2Cro)

end subroutine read_IC_from_MAR






