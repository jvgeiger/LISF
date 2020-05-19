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
!  29 Apr 2020: Zhuo Wang: Added MetSim forcing for LIS-SUMMA
!
! !INTERFACE:
subroutine read_metsim(n, kk, findex, order, mo, da, hr, filename, ferror)

! !USES:
   use LIS_coreMod,          only : LIS_rc, LIS_domain
   use LIS_logMod,           only : LIS_logunit, LIS_endrun, LIS_verify, &
                                    LIS_warning
   use LIS_metforcingMod,    only : LIS_forc
   use LIS_timeMgrMod,       only : LIS_tick
   use metsim_forcingMod,    only : metsim_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

   implicit none
! !ARGUMENTS:
   integer, intent(in)    :: n        ! nest
   integer, intent(in)    :: kk       ! forecast member index
   integer, intent(in)    :: findex   ! forcing index
   integer, intent(in)    :: order    ! lower(1) or upper(2) time interval bdry
   integer, intent(in)    :: mo
   integer, intent(in)    :: da
   integer, intent(in)    :: hr
   character(len=100), intent(in) :: filename ! name of input MetSim file
   integer, intent(inout) :: ferror   ! set to zero if there's an error
!
! !DESCRIPTION:
!  For the given time, reads the parameters from 0.125 degree
!  MetSim data, transforms into 7 LIS forcing
!  parameters and interpolates to the LIS domain.
!
!  MetSim variables used to force LIS are:
!  mean values starting at timestep, available every 3 hours \newline
!
!  NOTE 1: be aware that MetSim has only total precipitation.
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
!  \item[n]
!    index of the nest
!  \item[kk]
!    index of the forecast ensemble member
!  \item[findex]
!    index of the forcing source
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    3 hourly instance, order=2, read the next 3 hourly instance)
!  \item[mo]
!    current month of the year
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
   integer,parameter ::  ndims =  3          ! number of dimensions
   integer, parameter :: N_MS = 7 ! # of MetSim forcing variables

   character(len=12), dimension(N_MS), parameter :: metsim_fv = (/  &
      'temp        ',    &
      'spec_humid  ',    &
      'shortwave   ',    &
      'longwave    ',    &
      'wind        ',    &
      'air_pressure',    &
      'prec        '    /)
   character(len=12) :: varname


   integer :: ncid, varid, status
   integer :: timestep, day

   integer :: v, c, r, t, iret
   integer :: input_size
   real    :: gridDesco(50)         ! Input,output grid info arrays

   real,allocatable, dimension(:,:) :: vartemp
   real,allocatable, dimension(:) :: vartemp1d
   real, allocatable, dimension(:,:) :: varfield
   logical*1,allocatable :: input_bitmask(:)   ! input bitmap
   ! ______________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   ! If there is a problem, then ferror is set to zero
   ferror = 1
   input_size = metsim_struc(n)%ncold*metsim_struc(n)%nrold

   allocate(vartemp(metsim_struc(n)%ncold,metsim_struc(n)%nrold))
   allocate(varfield(LIS_rc%lnc(n), LIS_rc%lnr(n)))
   allocate(input_bitmask(metsim_struc(n)%ncold*metsim_struc(n)%nrold))
   allocate(vartemp1d(metsim_struc(n)%ncold*metsim_struc(n)%nrold))

   !----------------------------------------------------------------
   ! Establish which file timestep the date & hour correspond to
   ! and which file to open.
   !----------------------------------------------------------------

   ! MetSim data do not contain leap days.  If it is 29 Feb, then reset the
   ! day to 28 and reuse that data.
   if ( mo == 2 .and. da == 29 ) then
      day = 28
   else
      day = da
   endif

   ! For monthly file
   timestep = 8*(day - 1) + (1 + hr/3)

   status = nf90_open(filename, nf90_NoWrite, ncid)

   if ( status /= 0 ) then
      write(LIS_logunit,*)'[ERR] Problem opening file: ',trim(filename),status
      write(LIS_logunit,*)'[ERR]  Stopping...'
      call LIS_endrun
   else
      write(LIS_logunit,*)'[INFO] Opened file: ',trim(filename)
   endif

   do v = 1, N_MS
      varname = metsim_fv(v)
      call LIS_verify(nf90_inq_varid(ncid, trim(varname), varId), &
         'nf90_inq_varid failed for '//trim(varname)//' in read_metsim')
      call LIS_verify(nf90_get_var(ncid,varId,vartemp,start=(/1,1,timestep/), &
         count=(/metsim_struc(n)%ncold,metsim_struc(n)%nrold,1/)), &
         'nf90_get_var failed for '//trim(varname)//' in read_metsim')

      ! Filter out any unrealistic forcing values.
      if ( varname == 'prec' ) then
         do r=1,metsim_struc(n)%nrold
            do c=1,metsim_struc(n)%ncold
               if (vartemp(c,r) < 0.0) then
                  vartemp(c,r) = 0.0
               endif
            enddo
         enddo
      endif

      vartemp1d = reshape(vartemp, (/input_size/))
      ! TODO: review missingValue
      input_bitmask = .false.
      where ( vartemp1d /= LIS_rc%udef )
         input_bitmask = .true.
      endwhere

      call interp_metsim(n, findex, input_size, varname, &
         vartemp1d, input_bitmask, LIS_rc%gridDesc(n,:), &
         LIS_rc%lnc(n),LIS_rc%lnr(n),varfield)

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            t = LIS_domain(n)%gindex(c,r)
            if ( t /= -1 ) then
               metsim_struc(n)%metdata2(kk,v,t) = varfield(c,r)
            endif
         enddo
      enddo
   enddo


   ! Close netCDF file.
   status=nf90_close(ncid)
   deallocate(input_bitmask)
   deallocate(vartemp)
   deallocate(varfield)
#endif

end subroutine read_metsim

!BOP
! !ROUTINE: interp_metsim
! \label{interp_metsim}
!
! !INTERFACE:
subroutine interp_metsim(n, findex, input_size, input_name, input_data, &
      input_bitmask, lis_gds, nc, nr, varfield)
! !USES:
   use LIS_coreMod, only : LIS_rc, LIS_domain
   use metsim_forcingMod, only : metsim_struc

   implicit none

! !ARGUMENTS:
   integer, intent(in) :: n, findex, input_size, nc, nr
   character(len=12), intent(in) :: input_name
   real, intent(in) :: input_data(input_size)
   logical*1, intent(in) :: input_bitmask(input_size)
   real, intent(in) :: lis_gds(50)
   real, intent(out) :: varfield(nc, nr)

! !DESCRIPTION:
!   This subroutine interpolates a given MetSim field
!   to the LIS grid.
!  The arguments are:
!  \begin{description}
! \item[n]
!  index of the nest
! \item[findex]
!  index of MetSim in the list of all enabled metforcing sources.
! \item[input\_size]
!  size of the input\_data (as a 1-d array)
! \item[input\_data]
!  1-d array containing the MetSim field to spatially interpolate
! \item[input\_bitmask]
!  logical bitmask of the locations of valid input data
! \item[lis\_gds]
!  gridDesc array of the LIS running grid
! \item[nc]
!  number of columns (in the west-east dimension) in the LIS running grid
! \item[nr]
!  number of rows (in the south-north dimension) in the LIS running grid
! \item[varfield]
!  output array containing the spatially interpolated field
!  \end{description}

   integer :: iret
   integer :: output_size
   integer :: count1,i,j
   real, allocatable, dimension(:) :: output_data
   logical*1, allocatable, dimension(:) :: output_bitmask(:)

   output_size = nc*nr
   allocate(output_bitmask(output_size))
   allocate(output_data(output_size))
   output_bitmask = .true.


   select case( LIS_rc%met_interp(findex) )

   case( "bilinear" )
      call bilinear_interp(lis_gds,input_bitmask,input_data, &
         output_bitmask,output_data,input_size,output_size, &
         LIS_domain(n)%lat, LIS_domain(n)%lon,&
         metsim_struc(n)%w111,metsim_struc(n)%w121,&
         metsim_struc(n)%w211,metsim_struc(n)%w221,&
         metsim_struc(n)%n111,metsim_struc(n)%n121,&
         metsim_struc(n)%n211,metsim_struc(n)%n221,&
         LIS_rc%udef, iret)

   case( "budget-bilinear" )
      call conserv_interp(lis_gds,input_bitmask,input_data, &
         output_bitmask,output_data,input_size,output_size, &
         LIS_domain(n)%lat, LIS_domain(n)%lon,&
         metsim_struc(n)%w112,metsim_struc(n)%w122,&
         metsim_struc(n)%w212,metsim_struc(n)%w222,&
         metsim_struc(n)%n112,metsim_struc(n)%n122,&
         metsim_struc(n)%n212,metsim_struc(n)%n222,&
         LIS_rc%udef, iret)
   case( "neighbor" )
      call neighbor_interp(lis_gds,input_bitmask,input_data, &
         output_bitmask,output_data,input_size,output_size, &
         LIS_domain(n)%lat, LIS_domain(n)%lon,&
         metsim_struc(n)%n113,LIS_rc%udef, iret)
   end select

   ! convert the interpolated data to 2d.
   count1 = 0
   do j = 1, nr
      do i = 1, nc
         varfield(i,j) = output_data(i+count1)
      enddo
      count1 = count1 + nc
   enddo

   deallocate(output_bitmask)
   deallocate(output_data)
end subroutine interp_metsim
