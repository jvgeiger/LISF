!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_MAR
!  \label{read_MAR}
!
! !REVISION HISTORY:
! 22 Aug  2023 : Mahdi Navari ; Initial Specification
! 
! !INTERFACE:
subroutine read_MAR(n, findex, order, yr, mon, da, hr, ferror)
! !USES:
  use LIS_coreMod
  use LIS_metforcingMod, only : LIS_forc
  use LIS_timeMgrMod,       only : LIS_tick
  use LIS_logMod
  use LIS_constantsMod,     only : LIS_CONST_PATH_LEN
  use MAR_forcingMod, only : MAR_struc
  use LIS_mpiMod
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: findex
  integer, intent(in)          :: order
  integer, intent(in)    :: yr,mon,da,hr   
  integer, intent(out)         :: ferror
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  MAR files, transforms into 9 LIS forcing 
!  parameters.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing source
!  \item[order]
!    flag indicating which data to be read (order=1, read for the previous 
!    1hr bookend, order=2, read for the next 1hr bookend)
!  \item[fname]
!    name of the MAR file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
!EOP
  integer, parameter :: NF = 9   ! # of S forcing variables

  character(11), dimension(NF), parameter :: MAR_fv = (/  &
       'TTZ     ',  &    ! metdata(1) == tmp   @2m    (TIME, ZTQLEV, Y20_155, X10_85)   ZTQLEV = 2
       'QQZ     ',  &    ! metdata(2) == q2    @2m    (TIME, ZTQLEV, Y20_155, X10_85) 
       'SWD     ',  &    ! metdata(3) == swd          (TIME, Y20_155, X10_85)
       'LWD     ',  &    ! metdata(4) == lwd          (TIME, Y20_155, X10_85)
       'UUZ     ',  &    ! metdata(5) == uwind @10m   (TIME, ZUVLEV, Y20_155, X10_85)   ZUVLEV = 2
       'VVZ     ',  &    ! metdata(6) == vwind @10m   (TIME, ZUVLEV, Y20_155, X10_85)   
       'SP      ',  &    ! metdata(7) == psurf        (TIME, Y20_155, X10_85)
       'RF      ',  &    ! metdata(8) == rainf        (TIME, Y20_155, X10_85)
       'SF      '  /)    ! metdata(9) == snowf        (TIME, Y20_155, X10_85)

  character(len=LIS_CONST_PATH_LEN) :: infile

! netcdf variables
  integer :: ncid, varid, status
  integer :: latid, lonid

  real*8  :: timenow
  real    :: gmt
  integer :: doy,mn,ss,ts
  integer :: timestep
  character :: cyr*4
  character :: cmo*2
  character :: cda*2

  integer :: i, j, v, ii, r, k, eindex, x, y
  integer :: ios            ! set to non-zero if there's an error
  integer :: ierr
  integer :: iret, c
  real    :: gridDesco(50)         ! Input,output grid info arrays

  integer :: kk                    ! forecast index
  integer :: mo

  real,allocatable :: lat(:,:)     ! input data (longitude,latitude)
  real,allocatable :: lon(:,:)     ! input data (longitude,latitude)

  real,allocatable :: datain(:,:)  ! input data (longitude,latitude)
  real,allocatable :: temp2MAR(:,:,:)
  real,allocatable :: data_local_grid(:,:)
  !real,allocatable :: tempinterp(:,:,:)
  !real,allocatable :: f(:)         ! 1D in fields
  !real,allocatable :: go(:)        ! 1D out fields
  !real,allocatable :: tg(:,:)      ! Interpolated 2D data field
  !logical*1,allocatable :: lb(:)   ! input bitmap
  !logical*1,allocatable :: lo(:)   ! output bitmaps
  real,   parameter       ::  missing_value = -1.e+34  ! LIS_rc%udef

! __________________________________

#if 0
  integer                 :: c,r,gindex
  integer                 :: ftn, ios
  logical                 :: file_exists !, rainc_exists
  integer                 :: tmpId, qId, swdId, lwdId
  integer                 :: uId, vId, psfcId, rainId, snowId
  !integer                 :: raincID
  real                    :: gvar(LIS_rc%gnc(n),LIS_rc%gnr(n),1)
  real                    :: tair(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: qair(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: swd(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: lwd(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: u(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: v(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: psfc(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: rainf(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: snowf(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: tot_p(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,   parameter       ::  missing_value = -1.e+34  ! LIS_rc%udef
  real*8  :: timenow
  real    :: gmt
  integer :: doy,mn,ss,ts
  integer :: timestep
#endif

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   ! If a problem, ferror is set to zero
   ferror = 1
   mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

   ! Allocate memory
   !allocate(lat(MAR_struc(n)%nc,MAR_struc(n)%nr))
   !allocate(lon(MAR_struc(n)%nc,MAR_struc(n)%nr))
   allocate(datain(MAR_struc(n)%nc,MAR_struc(n)%nr))
   allocate(temp2MAR(LIS_rc%lnc(n),LIS_rc%lnr(n),NF)) ! MAR_struc(n)%nc,MAR_struc(n)%nr,NF))
   allocate(data_local_grid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
   !allocate(f(MAR_struc(n)%nc*MAR_struc(n)%nr))
   !allocate(go(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
   !allocate(lb(MAR_struc(n)%nc*MAR_struc(n)%nr))
   !allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
   !allocate(tg(LIS_rc%lnc(n),LIS_rc%lnr(n)))

   !allocate(tempinterp(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nf), stat=ios)
   !if(ios.ne.0) then
   !  write(LIS_logunit,*) '[ERR] Error allocating tempinterp,',LIS_localPet
   !  call LIS_endrun
   !endif

   temp2MAR = 0.0     ! initialize


  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour need to be read.
  !----------------------------------------------------------------
   mn=LIS_rc%mn ! Time of input file
   ss=0
   ts=0
   call LIS_tick(timenow,doy,gmt,yr,mon,da,hr,mn,ss,real(ts))

   ! One file per year--use day of year to get to the time record
   !timestep = 24*(doy - 1) + (1 + hr) ??????? why  (1 + hr)

   ! One file per month-- get to the time record
   ! start hr in monthly and yearly file is 1 (not 0)
   ! hr 0-23;  
   timestep = 24*(da - 1) + hr  ! - 1  ! + 1

   if(LIS_masterproc) then
      write(LIS_logunit,*)'[INFO] Order, timestep, doy, hr ::', order, timestep, doy, hr
   endif

   write(unit=cyr,fmt='(i4.4)') yr
   write(unit=cmo,fmt='(i2.2)') mon
     ! File name for data MARv3.12-ERA5-25km-hourly-year.nc
   infile = trim(MAR_struc(n)%MARdir)//'/'//&
          'ICE.'//trim(cyr)//trim(cmo)//'.b85.nc'
!          'ICE.'//trim(cyr)//'.01-12.b85.nc'
!         'MARv3.12-ERA5-25km-hourly-'//trim(cyr)//'.nc'

!=== Open MAR out forcing files ===
     ! Open netCDF file.
     status = nf90_open(trim(infile), nf90_NoWrite, ncid)
     if(status/=0) then
       if(LIS_masterproc) then
          write(LIS_logunit,*)'[ERR] Problem opening file: ',trim(infile),status
          call LIS_endrun
       endif
     else
       if(LIS_masterproc) then
         write(LIS_logunit,*)'[INFO] Opened file: ',trim(infile)
       endif
     end if

     ! Forcing variable loop:
     do v = 1, LIS_rc%met_nf(findex)  ! Number of met fields in MAR data

       datain = LIS_rc%udef
       status = nf90_inq_varid(ncid, trim(MAR_fv(v)), varid)
       if ((v .eq. 1) .or. (v .eq. 2)) then  
          status = nf90_get_var(ncid, varid, datain, &
                     start=(/1,1,1,timestep/), &
                     count=(/MAR_struc(n)%nc,MAR_struc(n)%nr,1,1/))

       elseif ((v .eq. 5) .or. (v .eq. 6)) then
          status = nf90_get_var(ncid, varid, datain, &
                     start=(/1,1,2,timestep/), &
                     count=(/MAR_struc(n)%nc,MAR_struc(n)%nr,1,1/))
       else 
          status = nf90_get_var(ncid, varid, datain, &
                     start=(/1,1,timestep/), &
                     count=(/MAR_struc(n)%nc,MAR_struc(n)%nr,1/))
      endif    
      data_local_grid(:,:) = datain(LIS_ews_halo_ind(n,LIS_localPet+1):&
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1))

      temp2MAR(:,:,v) = data_local_grid(:,:)
!print*,'MAR_fv(v)', v,  MAR_fv(v)
#if 0
       if( trim(WRFAK_fv(v)) == "T2" .and. LIS_rc%tscount(n) == 1 ) then
         ! Read in lat and lon info from the file:
         lat = LIS_rc%udef
         status = nf90_inq_varid(ncid, "XLAT", latid)
         status = nf90_get_var(ncid, latid, lat, &
                       start=(/1,1/), &
                       count=(/MAR_struc(n)%nc,MAR_struc(n)%nr/))

         lon = LIS_rc%udef
         status = nf90_inq_varid(ncid, "XLONG", lonid)
         status = nf90_get_var(ncid, lonid, lon, &
                       start=(/1,1/), &
                       count=(/MAR_struc(n)%nc,MAR_struc(n)%nr/))
       endif
      !-----------------------------------------------------------------
      ! Filter out any unrealistic forcing values.
      ! Transfer WRF out forcing fields to LIS format
      !-----------------------------------------------------------------
       do j=1,MAR_struc(n)%nr
         do i=1,MAR_struc(n)%nc

           select case (v)
            case (1) ! Tair
              temp2MAR(i,j,1) = datain(i,j)

            case (2) ! Qair
              temp2MAR(i,j,2) = datain(i,j)

            case (3) ! Shortwave
              if (datain(i,j) < 0.0001) then
                datain(i,j) = 0.0001
              endif
              temp2MAR(i,j,3) = datain(i,j)

            case (4) ! Longwave
              temp2MAR(i,j,4) = datain(i,j)

            case (5) ! U-vector wind
              if (datain(i,j) < 0.0001) then
                datain(i,j) = 0.0001
              endif
              temp2MAR(i,j,5) = datain(i,j)   

            case (6) ! V-vector wind
              if (datain(i,j) < 0.0001) then
                datain(i,j) = 0.0001
              endif
              temp2MAR(i,j,6) = datain(i,j)  

            case (7) ! Pressure
              temp2MAR(i,j,7) = datain(i,j)

            case (8) ! rainf
              if (datain(i,j) < 0.0) then
                datain(i,j) = 0.0
              endif
              temp2MAR(i,j,8) = datain(i,j)

            case (9) ! snowf
              if (datain(i,j) < 0.0) then
                datain(i,j) = 0.0
              endif
              temp2MAR(i,j,9) = datain(i,j)
           end select
         enddo
       enddo
#endif
     end do !v or variable loop

     ! Close netCDF file.
     call LIS_verify(nf90_close(ncid))


! The data is expected to be in the same map projection
! and resolution as the current LIS run.

! TODO: We can test and activate the following options in the future.
#if 0
! NOTE: Bilinear interpolation (bilinear_interp_input.F90), 
! conservative interpolation (conserv_interp_input.F90), 
! and neighbor interpolation (neighbor_interp_input.F90) issue calls to get_fieldpos.F90. However, 
! please be aware that this subroutine does not support polar stereographic projection.

     !--------------------------------------------------------------
     ! Interpolate each forcing variable to LIS domain
     !--------------------------------------------------------------

     !=== Initialize input & output grid arrays
     gridDesco = 0

     !=== Set input & output grid array values (WRF out to LIS)
     gridDesco = LIS_rc%gridDesc(n,:)

     !== valid value over land and ocean for WRFAK data
     lb = .false. ! .true.
     lo = .false.

     tempinterp(:,:,LIS_rc%nf) = 0.0

     do v=1,NF 

        !== Transferring current data to 1-D array for interpolation
        ! handling missing values  
        c=0
        do j=1,MAR_struc(n)%nr
           do i=1,MAR_struc(n)%nc
              c = c + 1
              if (temp2MAR(i,j,v) .ne. missing_value) then
                 f(c) = temp2MAR(i,j,v)
                 lb(c) = .true.
              else 
                 f(c) = LIS_rc%udef
              endif 
           enddo
        enddo

        !== Interpolate forcing to model grids
        select case( LIS_rc%met_interp(findex) )

          case( "bilinear" )
            call bilinear_interp(gridDesco,lb,f,lo,go,&
                 MAR_struc(n)%mi,mo, &
                 LIS_domain(n)%lat, LIS_domain(n)%lon,&
                 MAR_struc(n)%w111,MAR_struc(n)%w121,&
                 MAR_struc(n)%w211,MAR_struc(n)%w221,&
                 MAR_struc(n)%n111,MAR_struc(n)%n121,&
                 MAR_struc(n)%n211,MAR_struc(n)%n221,&
                 LIS_rc%udef, iret)

         case( "budget-bilinear" )
          if(v.eq.8 .or. v.eq.9) then
            call conserv_interp(gridDesco,lb,f,lo,go,&
                 MAR_struc(n)%mi,mo, &
                 LIS_domain(n)%lat, LIS_domain(n)%lon,&
                 MAR_struc(n)%w112,MAR_struc(n)%w122,&
                 MAR_struc(n)%w212,MAR_struc(n)%w222,&
                 MAR_struc(n)%n112,MAR_struc(n)%n122,&
                 MAR_struc(n)%n212,MAR_struc(n)%n222,&
                 LIS_rc%udef,iret)
          else
            call bilinear_interp(gridDesco,lb,f,lo,go,MAR_struc(n)%mi,mo, &
                 LIS_domain(n)%lat, LIS_domain(n)%lon,&
                 MAR_struc(n)%w111,MAR_struc(n)%w121,&
                 MAR_struc(n)%w211,MAR_struc(n)%w221,&
                 MAR_struc(n)%n111,MAR_struc(n)%n121,&
                 MAR_struc(n)%n211,MAR_struc(n)%n221,LIS_rc%udef, iret)
          endif

         case( "neighbor" )
           call neighbor_interp(gridDesco,lb,f,lo,go,&
                MAR_struc(n)%mi,mo,&
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                MAR_struc(n)%n113,LIS_rc%udef,iret)

        end select
        !== Convert data to original 3D array & a 2D array to 
        !== Fill in of missing points due to geography difference  
        tg = 0.0
        c = 0
        do j = 1, LIS_rc%lnr(n)
            do i = 1, LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(i,j).ne.-1) then
                tg(i,j) = go(i+c)
                if(tg(i,j) > 0. ) then
                endif
              endif
            enddo
            c = c + LIS_rc%lnc(n)
         enddo

         do j = 1, LIS_rc%lnr(n)
           do i = 1, LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(i,j).ne.-1) then
                if ((tg(i,j) .lt. 0) .and. (tg(i,j) .ne. LIS_rc%udef)) then
                  write(LIS_logunit,*)'[ERR] No nearest neighbors, v, i, j',v,i,j,tg(i,j)
                  call LIS_endrun()
                endif
                tempinterp(i,j,v) = tg(i,j)
              endif
           end do ! c
        enddo     ! r

    enddo  ! end v=variable loop
#endif

    do v= 1, NF
      ! Fill in undefined and ocean points
      do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r).ne.-1) then
              if(order.eq.1)then
                 MAR_struc(n)%metdata1(v,LIS_domain(n)%gindex(c,r))=temp2MAR(c,r,v) ! tempinterp(c,r,v) <- for activation of interpolation
              else
                 MAR_struc(n)%metdata2(v,LIS_domain(n)%gindex(c,r))=temp2MAR(c,r,v) !tempinterp(c,r,v) <- for activation of interpolation
              endif
           endif
        enddo !c
      enddo   !r
    enddo     !v

  !end do  ! End forecast member index loop

  ! Deallocate local interp-based variables:
  deallocate(datain)
  deallocate(temp2MAR)
  deallocate(data_local_grid) 
  !deallocate(f)
  !deallocate(go)
  !deallocate(lb)
  !deallocate(lo)
  !deallocate(tg)
  !deallocate(tempinterp)

! NetCDF
#endif  
end subroutine read_MAR
