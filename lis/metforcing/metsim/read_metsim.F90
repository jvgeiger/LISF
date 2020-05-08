!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_metsim 
!  \label{read_metsim}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Beaudoing; Initial Specification adopted from 
!                                 LIS/retberg.F90
!  16 Feb 2016: Hiroko Beaudoing; Fixed indexes for precip and pressure fields 
!                                 so budget-bilinear interpolation applies to 
!                                 correct individual fields.
!  15 May 2017: Bailing Li; Added changes for reading in version 2.2 data
!                           that is in 3D array (4D in version 1 & 2).
!  22 Oct 2018: Daniel Sarmiento; Added changes to support version 3 data
!  29 Apr 2020: Zhuo Wang: Added MetSim forcing for LIS-SUMMA 
!
! !INTERFACE:
subroutine read_metsim( order, n, findex, yr, mon, da, hr, ferror )

! !USES:
  use LIS_coreMod,          only : LIS_rc,LIS_domain, LIS_localPet, &
                                   LIS_masterproc
  use LIS_metforcingMod,    only : LIS_forc
  use LIS_timeMgrMod,       only : LIS_tick
  use LIS_logMod,           only : LIS_logunit, LIS_endrun, LIS_verify, LIS_warning
  use metsim_forcingMod,    only : metsim_struc
  use LIS_forecastMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: order     ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: n         ! nest
  integer, intent(in)    :: findex    ! forcing index
  integer, intent(in)    :: yr,mon,da,hr     ! data and hour (multiple of 3)
  integer, intent(inout) :: ferror           ! set to zero if there's an error
!
! !DESCRIPTION:
!  For the given time, reads the parameters from 0.125 degree
!  MetSim data, transforms into 7 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  MetSim variables used to force LIS are:
!  mean values starting at timestep, available every 3 hours \newline
!
!  NOTE-1: be aware that MetSim has only total precipitation.
!  NOTE 2: only one wind component, it is magnitude \newline
!
!  MetSim FORCING VARIABLES: \newline
!  1. airtemp   Air Temperature [K] \newline
!  2. spechum   Specific humidity [g g-1] \newline
!  3. SWRadAtm  Downward shortwave radiation [W m-2] \newline
!  4. LWRadAtm  Downward longwave radiation [W m-2] \newline
!  5. windspd   Wind speed [m s-1] \newline
!  6. airpres   Surface pressure [Pa] \newline
!  7. pptrate   Precipitation [kg m-2 s-1] \newline
!
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    3 hourly instance, order=2, read the next 3 hourly instance)
!  \item[n]
!    index of the nest
!  \item[yr]
!    current year
!  \item[mon]
!    current month
!  \item[da]
!    current day of the year
!  \item[hr]
!    current hour of day
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[metsimgrid\_2\_lisgrid](\ref{metsimgrid_2_lisgrid}) \newline
!    transform the MetSim data to the LIS grid 
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!
!EOP

! Specify MetSim forcing parameters & file paths
  integer,parameter ::  nvars =  7          ! number of variables: last
                                             ! variable is of interest
  integer,parameter ::  ndims =  3          ! number of dimensions

  integer,parameter ::  z =      1
  real,   parameter ::  no_data_value = 2.e+20

  integer, parameter :: N_PF = 7 ! # of MetSim forcing variables
  integer, parameter :: NF = 9   ! # of GLDAS forcing variables

  character(12), dimension(N_PF), parameter :: metsim_fv = (/  &
!       'tas         ',    &
!       'shum        ',    &
!       'dswrf       ',    &
!       'dlwrf       ',    &
!       'wind        ',    &
!       'pres        ',    &
!       'prcp        '     /)

!      'airtemp     ',    &
!      'spechum     ',    &
!      'SWRadAtm    ',    &
!      'LWRadAtm    ',    &
!      'windspd     ',    &
!      'airpres     ',    &
!      'pptrate     '     /)

       'temp        ',    &
       'spec_humid  ',    &
       'shortwave   ',    &
       'longwave    ',    &
       'wind        ',    &
       'air_pressure',    &
       'prec        '    /)

! integer, dimension(N_PF), parameter :: fnum = (/ &
!      31, 33, 34, 35, 36, 37, 38 /) 

  character*100 :: infile

! netcdf variables
  integer :: ncid, varid, status

  integer :: mo
  real*8  :: timenow
  real    :: gmt
  integer :: doy,mn,ss,ts
  character :: cyr*4, cmon*4
  integer :: timestep

  integer :: i, j, v, ii, r, k, eindex, x, y
  integer :: ios                   ! set to non-zero if there's an error
  integer :: gldas, nmetsim        ! Size of I/O 1D fields
  integer :: iret, c
  real    :: gridDesco(50)         ! Input,output grid info arrays

! integer :: kk                    ! forecast index
  integer :: num_met               ! number of fields

  real,allocatable :: temp(:,:)
  real,allocatable :: spechum(:,:)
  real,allocatable :: shortwave(:,:)
  real,allocatable :: longwave(:,:)
  real,allocatable :: wind(:,:)
  real,allocatable :: air_pressure(:,:)
  real,allocatable :: prec(:,:)
! real,allocatable :: datain(:,:)  ! input data (longitude,latitude)
  real,allocatable :: temp2metsim(:,:,:)
! real,allocatable :: templdas(:,:,:)
! real,allocatable :: f(:)         ! 1D in fields
! real,allocatable :: go(:)        ! 1D out fields
! real,allocatable :: tg(:,:)      ! Interpolated 2D data field
! logical*1,allocatable :: lb(:)   ! input bitmap
! logical*1,allocatable :: lo(:)   ! output bitmaps
! ______________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   ! If there is a problem, then ferror is set to zero
   ferror = 1
   mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
   num_met = 7

   ! Allocate memory
!  allocate(datain(metsim_struc(n)%ncold,metsim_struc(n)%nrold))
   allocate(temp2metsim(metsim_struc(n)%ncold,metsim_struc(n)%nrold,NF))
!  allocate(f(metsim_struc(n)%ncold*metsim_struc(n)%nrold))
!  allocate(go(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(lb(metsim_struc(n)%ncold*metsim_struc(n)%nrold))
!  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(tg(LIS_rc%lnc(n),LIS_rc%lnr(n)))  
!  allocate(templdas(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nf), stat=ios)
 
   temp2metsim = 0.0     ! initialize

  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour correspond to
  ! and which file to open.
  !----------------------------------------------------------------
   mn=LIS_rc%mn ! Time of input file
   ss=0
   ts=0

   ! Determine the time, time in GMT, and the day of year based on 
   ! the value of year, month, day of month, hour of the day,
   !  the day, minute and second.
   call LIS_tick(timenow,doy,gmt,yr,mon,da,hr,mn,ss,real(ts))

   ! One file per year--use day of year to get to the time record
   ! MetSim is one file per month. Need to revise to get the current timestep ???
   ! yr,mon,da,hr are input !!! 

! For annual file
!  timestep = 8*(doy - 1) + (1 + hr/3)

! For monthly file 
   timestep = 8*(da - 1) + (1 + hr/3)

   if(LIS_masterproc) then
!     write(LIS_logunit,*)'[INFO] Order and month-timestep', order, timestep, doy, hr
      write(LIS_logunit,*)'[INFO] Order and month-timestep', order, timestep, mon, da, hr
   endif

!=== Open MetSim forcing files ===

  ! Loop over forecast index:
  do kk= metsim_struc(n)%st_iterid, metsim_struc(n)%en_iterid

     ! Modified by KRA to implement forecast mode:
!    if(LIS_rc%forecastMode.eq.0) then ! hindcast run
        write(cyr, '(i4.4)') yr
        write(cmon, '(i4.4)') mon

!    else !forecast mode
!       !sample yr, mo, da
!       call LIS_sample_forecastDate(n, kk, findex, yr, mo, da) 
!       write(cyr, '(i4.4)') yr
!    endif

     ! Forcing variable loop:
     ! ??? If one monthly file includes all forcing variables, then we don't need to 
     ! do variable loop
     ! Check ../merra-land/read_merraland.F90 (lines 125-182) how to read all variables from one file

!    do v = 1, num_met  ! Number of met fields in MetSim data

       ! File name for data year/variable(v)_3hourly_year-year.nc
       ! ???? Forcing directory ????
!      infile=trim(metsim_struc(n)%metsimdir)//'/'//cyr//'/'//&
!             trim(metsim_fv(v))//'_3hourly_'//cyr//'-'//cmon//'.nc'
       infile=trim(metsim_struc(n)%metsimdir)//'/'&
              //'icar_canesm2_'//cyr//'-'//cmon//'.nc'

       ! Open netCDF file.
       status = nf90_open(infile, nf90_NoWrite, ncid)
!      status = nf90_inq_varid(ncid, trim(metsim_fv(v)), varid)
       status = nf90_inq_varid(ncid, 'temp', tempId)
       status = nf90_inq_varid(ncid, 'spec_humid', spechumId)
       status = nf90_inq_varid(ncid, 'shortwave', shortwaveId)
       status = nf90_inq_varid(ncid, 'longwave', longwaveId)
       status = nf90_inq_varid(ncid, 'wind', windId)
       status = nf90_inq_varid(ncid, 'air_pressure', airPresId)
       status = nf90_inq_varid(ncid, 'prec',precId)

       if(status/=0) then
         if(LIS_masterproc) then 
            write(LIS_logunit,*)'[ERR] Problem opening file: ',infile,status
            write(LIS_logunit,*)'[ERR]  Stopping...'
            call LIS_endrun
         endif
         call LIS_endrun
       else
         if(LIS_masterproc) write(LIS_logunit,*)'[INFO] Opened file: ',infile
       end if

!      status = nf90_get_var(ncid, varid, datain,start=(/1,1,timestep/), & 
!      count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/))

       status = nf90_get_var(ncid, tempId, temp,start=(/1,1,timestep/), &
       count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/))

       status = nf90_get_var(ncid, spechumId, spechum,start=(/1,1,timestep/), &
       count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/))

       status = nf90_get_var(ncid, shortwaveId, shortwave,start=(/1,1,timestep/), &
       count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/))

       status = nf90_get_var(ncid, longwaveId, longwave,start=(/1,1,timestep/), &
       count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/))

       status = nf90_get_var(ncid, windId, wind,start=(/1,1,timestep/), &
       count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/))

       status = nf90_get_var(ncid, airPresId,air_pressure,start=(/1,1,timestep/), &
       count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/))

       status = nf90_get_var(ncid, precId, prec,start=(/1,1,timestep/), &
       count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/))

       ! Close netCDF file.
       status=nf90_close(ncid)

      !----------------------------------------------------------------
      ! Change data from MetSim grid convention to LIS-domain
      ! Shift longitudes by 180deg. 
      !----------------------------------------------------------------
      ! Filter out any unrealistic forcing values.
      ! Transfer MetSim forcing fields to LIS format
      !-----------------------------------------------------------------
     do v = 1, num_met  ! Number of met fields in MetSim data
       do j=1,metsim_struc(n)%nrold
         do i=1,metsim_struc(n)%ncold

           select case (v)
            case (1)! temp
              temp2metsim(i,j,1) = temp(i,j)
            case (2)! spec_humid 
              temp2metsim(i,j,2) = spechum(i,j)
            case (3)! shortwave
              if (shortwave(i,j) < 0.0001) then
                shortwave(i,j) = 0.0001
              endif
              temp2metsim(i,j,3) = shortwave(i,j)
            case (4)! longwave 
              temp2metsim(i,j,4) = longwave(i,j)
            case (5)! wind
              if (wind(i,j) < 0.0001) then
                wind(i,j) = 0.0001
              endif
              temp2metsim(i,j,5) = wind(i,j)  !Since absolute wind speed 
                                                   ! let U=WIND and V=0.0
            case (6)! air_pressure 
              temp2metsim(i,j,7) = air_pressure(i,j)    
            case (7)! Total precipitation prec: convective precp=0
              if (prec(i,j) < 0.0) then
                prec(i,j) = 0.0
              endif
           end select
         enddo
       enddo
     end do !v or variable loop

     !--------------------------------------------------------------
     ! Interpolate each forcing variable to LIS domain
     !--------------------------------------------------------------

     !=== Initialize input & output grid arrays
!    gridDesco = 0
           
     !=== Set input & output grid array values (reanlMetSim to LIS)
!    gridDesco = LIS_rc%gridDesc(n,:)

     !=== Define input & output data bitmaps
!    nmetsim = metsim_struc(n)%ncold*metsim_struc(n)%nrold
!    gldas  = LIS_rc%lnc(n)*LIS_rc%lnr(n)

     !== valid value over land and ocean for MetSim data
!    lb = .true.
!    lo = .false.
   !  tmask = .false.

!    templdas(:,:,9) = 0.0

!    do v=1,NF-1      ! do not process convective precip field, which is set to 0
!       if (v .ne. 6) then ! not the v-wind component, which is set to zero.
      
!         !=== Transferring current data to 1-D array for interpolation
!         c=0
!         do i=1,metsim_struc(n)%nrold
!            do j=1,metsim_struc(n)%ncold
!               c = c + 1
!               f(c) = temp2metsim(j,i,v)
!                if(v.eq.1 .and. f(c).ne.-9999.0) write(LIS_logunit,*) c,f(c)
!            enddo
!         enddo

          !=== Interpolate if forcing and model grids are not both one deg.
          !=== ???? Is 1 degree grid for MerSim data?

!         if( LIS_rc%gridDesc(n,9) .ne. 1.0 ) then   ???
!         if( LIS_rc%gridDesc(n,9) .ne. 0.125) then

!          select case( LIS_rc%met_interp(findex) )

!            case( "bilinear" )
!              call bilinear_interp(gridDesco,lb,f,lo,go,&
!                 metsim_struc(n)%mi,mo, & 
!                 LIS_domain(n)%lat, LIS_domain(n)%lon,&
!                 metsim_struc(n)%w111,metsim_struc(n)%w121,&
!                 metsim_struc(n)%w211,metsim_struc(n)%w221,& 
!                 metsim_struc(n)%n111,metsim_struc(n)%n121,&
!                 metsim_struc(n)%n211,metsim_struc(n)%n221,&
!                 LIS_rc%udef, iret)

!            case( "budget-bilinear" )
!              if(v.eq.8) then 
!                call conserv_interp(gridDesco,lb,f,lo,go,&
!                    metsim_struc(n)%mi,mo, & 
!                    LIS_domain(n)%lat, LIS_domain(n)%lon,&
!                    metsim_struc(n)%w112,metsim_struc(n)%w122,&
!                    metsim_struc(n)%w212,metsim_struc(n)%w222,& 
!                    metsim_struc(n)%n112,metsim_struc(n)%n122,&
!                    metsim_struc(n)%n212,metsim_struc(n)%n222,&
!                    LIS_rc%udef,iret)
!              else
!                call bilinear_interp(gridDesco,lb,f,lo,go,metsim_struc(n)%mi,mo, & 
!                    LIS_domain(n)%lat, LIS_domain(n)%lon,&
!                    metsim_struc(n)%w111,metsim_struc(n)%w121,&
!                    metsim_struc(n)%w211,metsim_struc(n)%w221, & 
!                    metsim_struc(n)%n111,metsim_struc(n)%n121,&
!                    metsim_struc(n)%n211,metsim_struc(n)%n221,LIS_rc%udef, iret)
!              endif

!          end select

!        else ! forcing and model grids both 0.125 degree
!           k = 0
!           do r=1,LIS_rc%lnr(n)
              !y = r + (LIS_rc%gridDesc(n,4) + 89.50) / 1.000
!              y = r + (LIS_rc%gridDesc(n,4) + 52.9375) / 0.125    ! ????
!              do c=1,LIS_rc%lnc(n)
                 !x = c + (LIS_rc%gridDesc(n,5) + 179.50) / 1.000
!                 x = c + (LIS_rc%gridDesc(n,5) - 67.0625) / 0.125 ! ???
!                 k = k + 1
                 !eindex = ((y - 1) * 360) + x
!                 eindex = ((y - 1) * 464) + x     ! ???
!                 go(k) = f(eindex)
!                 lo(k) = lb(eindex)
!              end do
!           end do

!        end if ! LIS_rc%domain 
    
         !=== Convert data to original 3D array & a 2D array to 
         !=== fill in of missing points due to geography difference  
!         tg = -9999.0
!        tg = 0.0 
!        c = 0
!        do j = 1, LIS_rc%lnr(n)
!           do i = 1, LIS_rc%lnc(n)
!              write(LIS_logunit,*) i,j,gindex(i,j),lo(i+c)
!             if(LIS_domain(n)%gindex(i,j).ne.-1) then 
!                geogmask(i,j) = lo(i+c)
!               tg(i,j) = go(i+c)
!             endif
!           enddo
!           c = c + LIS_rc%lnc(n)
!        enddo

         !== no need to fill in for MetSim
         ! call geogfill2(n, LIS_rc%lnc(n),LIS_rc%lnr(n),geogmask,tg,v,tmask)
         !       write(LIS_logunit,*) gindex(21,3),tg(21,3),v

!        do j = 1, LIS_rc%lnr(n)
!          do i = 1, LIS_rc%lnc(n)
!             if(LIS_domain(n)%gindex(i,j).ne.-1) then 
!               if ((tg(i,j) .lt. 0) .and. (tg(i,j) .ne. LIS_rc%udef)) then
!                  !For version 3, replace the value with an undefined value instead of ending the run.
!                  !For all other versions, keep old call to end the program.
!                  if(metsim_struc(n)%version == "3") then
!                     write(LIS_logunit,*)'[WARN] No nearest neighbors, v, i, j',v,i,j,tg(i,j)
!                     write(LIS_logunit,*)'[WARN] Check output to make sure the missing neighbor'
!                     write(LIS_logunit,*)'is isolated and does not distort the output.'      
!                     tg(i,j) = LIS_rc%udef ! New code to change data to undefined
!                  else
!                     write(LIS_logunit,*)'[WARN] No nearest neighbors, v, i, j',v,i,j,tg(i,j)
!                     call LIS_endrun()  ! Old code to call a fatal error
!                  endif
!               endif
!               templdas(i,j,v) = tg(i,j)
!             endif
!          end do ! c
!       enddo    ! r
!
!     else ! v==6, v-wind component, always zero
!        templdas(:,:,6) = 0.0
!     endif
!   enddo  ! end v=variable loop

!   do v= 1,NF
      ! Fill in undefined and ocean points
!     do r = 1,LIS_rc%lnr(n)
!       do c = 1,LIS_rc%lnc(n)
!         if (LIS_domain(n)%gindex(c,r).ne.-1) then
!          if(order.eq.1)then
!             metsim_struc(n)%metdata1(kk,v,LIS_domain(n)%gindex(c,r))=templdas(c,r,v)
!          else
!             metsim_struc(n)%metdata2(kk,v,LIS_domain(n)%gindex(c,r))=templdas(c,r,v)
!          endif
!         endif
!       enddo !c
!     enddo   !r
!   enddo     !v

! end do  ! End forecast member index loop

  ! Deallocate local interp-based variables:
! deallocate(datain)
  deallocate(temp2metsim)
! deallocate(f)
! deallocate(go)
! deallocate(lb)
! deallocate(lo)
! deallocate(tg)
! deallocate(templdas)
#endif

end subroutine read_metsim
