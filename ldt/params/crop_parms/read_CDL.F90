!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_CDL
!  \label{read_CDL}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!  28 Mar 2023: J. Erlingis: Modify UMD crop reader for USDA Cropland
!                            Data Layer (CDL)
! !INTERFACE:
subroutine read_CDL(n, num_types, fgrd)

! !USES:
  use ESMF
  use LDT_coreMod,     only : LDT_config, LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_LSMCropModifier_Mod

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMCrop_struc(n)%croptype%vlevels)

!
! !DESCRIPTION:
!  This subroutine reads the USGS Cropland Data Layer data and returns 
!  the distribution of vegetation in each grid cell, in a lat/lon
!  projection. The CDL varies year-to-year, so the input file generated
!  by LDT must have a time dimension. This reader handles data that has
!  already been resampled to 0.01 x 0.01 degrees instead of the native
!  30-meter CDL data.  
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[begin_year]
!     beginning year of CDL data
!   \item[end_year]
!     ending year of CDL data
!   \item[array1]
!     output array with the dominant crop type in each grid cell
!   \end{description}
!EOP      

   integer :: ftn, ierr, ios1
   logical :: file_exists
   integer :: begin_year, end_year                  ! CDL varies year-to-year
   integer :: jj, i, t, c, r, line, counter
   integer :: glpnc, glpnr                          ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr                        ! Parameter subsetted columns and rows
   integer :: rc
   integer :: mi                                    ! Total number of input param grid array points
   integer :: mo                                    ! Total number of output LIS grid array points
   integer*2 :: read_cropmap(6254,2866              ! Input grid)
   real    :: param_gridDesc(20)                    ! Input parameter grid desc fgrd
   real    :: subparam_gridDesc(20)                 ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)

   integer   :: file_size
   integer   :: file_dim
   integer   :: lonid, latid, cropid                ! CDL netCDF variable IDs
   integer   :: nrows, ncols                        ! spatial dimensions of CDL
   character(300) :: tempfile
   character(15)   :: myyear
!__________________________________________________________________

! Cropland Data Layer (postprocessed and aggregated to 0.01 degrees)
   integer, parameter :: IN_cols = 6254
   integer, parameter :: IN_rows = 2866
   real,    parameter :: IN_xres = 0.01 ! degrees
   real,    parameter :: IN_yres = 0.01 ! degrees
 
   begin_year = LDT_LSMCrop_struc(n)%begin_year
   end_year = LDT_LSMCrop_struc(n)%end_year

   !if( allocated(domcdl) ) then
   !   deallocate(domcdl)
   !endif

   !allocate(domcdl(LDT_rc%lnc(n),LDT_rc%lnr(n),((end_year-begin_year)+1)))
   !domcdl = LDT_rc%udef

   !water_class  = 33

   !vegcnt  = 0.
   !vegtype = float(water_class)
   !crop_array = 0.0
   !fgrd    = 0.0

!- Set parameter grid array inputs:
   LDT_LSMCrop_struc(n)%crop_proj = "latlon"
   param_gridDesc(1)  = 0.          ! Latlon
   param_gridDesc(2)  = 6254 !701 !5781 !6254
   param_gridDesc(3)  = 2866 !401 !2385 !2866
   param_gridDesc(4)  = 22.9451216 !39.9501215683218 !25.1101215683218 !22.9451216     ! LL lat 
   param_gridDesc(5)  = -127.8872122 !-96.8822121796901 !-124.88221217969 !-127.8872122   ! LL lon 
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = 51.6051216 !43.9501215683218 !48.9501215683218 !51.6051216     ! UR lat
   param_gridDesc(8)  = -65.3472122 !-89.8822121796901 !-67.0822121796901 !-65.3472122    ! UR lon
   param_gridDesc(9)  = 0.01
   param_gridDesc(10) = 0.01
   param_gridDesc(20) = 64

! Loop through years

!fgrd(:,:,:) = LDT_rc%udef
counter = 0

do jj=begin_year,end_year
   counter=counter+1
   print*,jj
   write(myyear,'(i4)') jj
   print*,myyear
   tempfile = trim(LDT_LSMCrop_struc(n)%croptfile)//trim(myyear)//"_1km_cdls_wgs84_int16.nc"
   !- Check if land cover file exists:
      inquire( file=trim(tempfile), exist=file_exists ) 
      if(.not. file_exists) then 
         write(LDT_logunit,*) "Crop type map: ",trim(tempfile)," does not exist. "
         write(LDT_logunit,*) "program stopping ..."
         call LDT_endrun
      endif
      write(LDT_logunit,*)"[INFO] Reading CDL: ",trim(tempfile)
      write(LDT_logunit,*)"** NOTE: Only 'mode' spatial transform currently works"
      write(LDT_logunit,*)"         with the 'CDL' crop type map source." 

      ftn = LDT_getNextUnitNumber()

      !if( allocated(read_cropmap) ) then
      !    deallocate(read_cropmap)
      !endif

      #if (defined USE_NETCDF3 || defined USE_NETCDF4)
         ierr = nf90_open(path=trim(tempfile),mode=NF90_NOWRITE,ncid=ftn)
         call LDT_verify(ierr,'error opening CDL file '//trim(tempfile))

         ierr = nf90_inq_dimid(ftn,'lon',lonid)
         call LDT_verify(ierr,'nf90_inq_dimid failed for longitude in read_CDL')
         ierr = nf90_inq_dimid(ftn,'lat',latid)
         call LDT_verify(ierr,'nf90_inq_dimid failed for latitude in read_CDL')
         ierr = nf90_inquire_dimension(ftn,lonid,len=ncols)
         call LDT_verify(ierr,'nf90_inquire_dimension for longitude')
         ierr = nf90_inquire_dimension(ftn,latid,len=nrows)
         call LDT_verify(ierr,'nf90_inquire_dimension for latitude')

         print*,nrows
         print*,ncols

         !allocate( read_cropmap(nrows, ncols) )
         !read_cropmap = LDT_rc%udef

         ierr = nf90_inq_varid(ftn,'Band1',cropid)
         call LDT_verify(ierr, 'nf90_inq_varid failed for crop type in read_CDL')
         ierr = nf90_get_var(ftn, cropid, read_cropmap)
         call LDT_verify(ierr, 'nf90_get_var failed for read_CDL')
         ierr = nf90_close(ftn)
         call LDT_verify(ierr, 'nf90_close failed in read_CDL')

         !print*,read_cropmap
         !
         ! Need to subset here
         !

         !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
         subparam_gridDesc = 0.
         call LDT_RunDomainPts( n, LDT_LSMCrop_struc(n)%crop_proj, param_gridDesc, &
                     glpnc, glpnr, subpnc, subpnr,  &
                     subparam_gridDesc, lat_line, lon_line )

!         print*,'glpnc ',glpnc
!         print*,'glpnr ',glpnr
!         print*,'subpnc ',subpnc
!         print*,'subprn ',subpnr
!         print*,'latline ',lat_line
!         print*,'lonline ',lon_line
         !array1(:,:,jj) = read_cropmap
       #endif

       do r = 1, subpnr
          do c = 1, subpnc
             fgrd(c,r,counter) = LDT_rc%udef
             !print*,r
             !print*,c
             !print*,lon_line(c,r)
             !print*,lat_line(c,r)
             !print*,read_cropmap(lon_line(c,r),lat_line(c,r))
             fgrd(c,r,counter) = float(read_cropmap(lon_line(c,r),lat_line(c,r)))
             !print*,fgrd(c,r,counter)
          end do ! 
       end do !
end do !Loop over years

!print*,fgrd

!JE Commenting out rescaling for now
!   ! -------------------------------------------------------------------
!   !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
!   ! -------------------------------------------------------------------
!   !- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
!      subparam_gridDesc = 0.
!      call LDT_RunDomainPts( n, LDT_LSMCrop_struc(n)%crop_proj, param_gridDesc, &
!                     glpnc, glpnr, subpnc, subpnr,  &
!                     subparam_gridDesc, lat_line, lon_line )
!
!   !- Determine if crop map file is 2D or 3D:
!      inquire( file=trim(LDT_LSMCrop_struc(n)%croptfile), size=file_size )
!      if( file_size > glpnc*glpnr*4 ) then
!         write(unit=LDT_logunit,fmt=*) '[INFO] Opening/Reading 3-D crop map ...'
!         file_dim = 3   ! 3-D fields
!      else 
!         file_dim = 2   ! 2-D fields
!      end if
!
!      if( file_dim == 3 .and. LDT_LSMCrop_struc(n)%crop_gridtransform == "mode" ) then  
!        if( subparam_gridDesc(9) .ne. (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)) ) then
!          write(*,*) "[WARN] IF 3-D TILED CROP MAP DOES NOT HAVE THE SAME RESOLUTION "
!          write(*,*) "  AS THE LIS RUN-DOMAIN, THEN 'MODE' OPTION CANNOT BE SELECTED,"
!          write(*,*) "  OR YOU NEED TO SELECT 'TILE' AS YOUR OPTION."
!          write(*,*) " Stopping ..."
!          call LDT_endrun
!        endif
!      endif
!enddo ! Year loop
!! -------------------------------------------------------------------
!!    READ IN LAND COVER PARAMETER FIELDS (NON-TILED/TILED OPTIONS)
!! -------------------------------------------------------------------
!
!!- Initialize parameter read-in array:
!  select case ( LDT_LSMCrop_struc(n)%crop_gridtransform ) 
!
!!    case ( "mode", "tile" )
!    case ( "mode" )
!      line = 0
!      do t = 1, num_types
!         do r = 1, subpnr
!            do c = 1, subpnc
!               line = (lat_line(c,r)-1)*glpnc + lon_line(c,r) + (t-1)*(glpnc*glpnr)
!               read(ftn,rec=line) vegcnt(c,r,t)
!            enddo
!         enddo
!      enddo
!
!    case default
!       write(LDT_logunit,*)" No other aggregation types are currently supported for UMD/CROPLAND." 
!       write(LDT_logunit,*)" -- Program stopping ..."
!       call LDT_endrun
!    end select  
!    deallocate( lat_line, lon_line )
!
!! -------------------------------------------------------------------
!!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
!! -------------------------------------------------------------------
!   mi = subpnc*subpnr
!   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
!   allocate( gi(mi), li(mi) )
!   !gi = float(water_class)
!   li = .false.
!   lo1 = .false.;  lo2 = .false.
!
!   select case( LDT_LSMCrop_struc(n)%crop_gridtransform ) 
!
!  !- (a) Estimate dominant landcover/crop types:
!     case( "mode" )
!
!    !- Aggregate:
!       if( file_dim == 2 ) then   ! 2-D fields
!  
!       !- Transform parameter from original grid to LIS output grid:
!          call LDT_transform_paramgrid(n, LDT_LSMCrop_struc(n)%crop_gridtransform, &
!                   subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )
!
!       !- Convert 1D dominant veg type to 2D grid arrays:
!          i = 0
!          do r = 1, LDT_rc%lnr(n)
!             do c = 1, LDT_rc%lnc(n)
!                i = i + 1
!                vegtype(c,r) = go1(i)
!             enddo
!          enddo
! 
!    !- No spatial aggregation needed; just dom. class selected:
!       elseif( file_dim == 3 ) then   ! 3-D fields
!
!          do r = 1, subpnr
!             do c = 1, subpnc
!                do t = 1, num_types
!                   if( t > 13 ) crop_array(t) = vegcnt(c,r,t)
!                end do
!                if( sum(crop_array(14:num_types)) == 0. ) then
!                  fgrd(c,r,1) = LDT_rc%udef
!                else
!                  fgrd(c,r,1) = maxloc(crop_array,1)
!                end if
!             enddo
!          enddo
!          return  ! Return dominant vegtype locations as a class to core routine
!
!       endif  ! End 2D/3D check
!
!!     case( "tile" )
!     case default
!       write(*,*) "[WARN] Other spatial grid transformations are not currently supported  "
!       write(*,*) "  for the tiled UMD-CROPLAND landcover/crop map type.  Please select either:"
!!       write(*,*) "  -- Mode or Tile "
!       write(*,*) "  -- Mode "
!       write(*,*) " Stopping ..."
!       call LDT_endrun
!   end select  ! End vegtype/cnt aggregation method
!   deallocate( gi, li )

!! ........................................................................

!!- Bring 2-D Vegtype to 3-D Vegcnt tile space:
!   if ( LDT_LSMCrop_struc(n)%crop_gridtransform == "none" .or. &
!        LDT_LSMCrop_struc(n)%crop_gridtransform == "mode" ) then  ! -- NON-TILED SURFACES
!     do r = 1, LDT_rc%lnr(n)
!        do c = 1, LDT_rc%lnc(n)
!           if ( vegtype(c,r) .le. 0 ) then
!              vegtype(c,r) = float(water_class)
!           endif
!           if ( (nint(vegtype(c,r)) .ne. water_class ) .and. &
!                (nint(vegtype(c,r)) .ne. LDT_rc%udef)) then
!              vegcnt(c,r,NINT(vegtype(c,r))) = 1.0
!           endif
!        enddo
!     end do
!   endif   ! End NON-TILED vegetation option
!
!!- Estimate fraction of grid (fgrid) represented by vegetation type::
!   call param_index_fgrdcalc( n, LDT_LSMCrop_struc(n)%crop_proj, &
!        LDT_LSMCrop_struc(n)%crop_gridtransform, &
!        water_class, num_types, vegcnt, fgrd )
!
!   call LDT_releaseUnitNumber(ftn)

   write(LDT_logunit,*) "[INFO] Done reading CDL. "


end subroutine read_CDL
