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
! !ROUTINE: read_PlantingDay
!  \label{read_PlantingDay}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  21 Feb 2025: P.-W. Liu: Modify the reader for planting day
!
! !INTERFACE:
subroutine read_PlantingDay(n, num_bins, fgrd) !PlantingDay Layer

! !USES:
  use ESMF
  use LDT_coreMod,     only : LDT_config, LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_plantingdayMod

#if (defined USE_GDAL)
  use, INTRINSIC:: iso_c_binding
  use fortranc
  use gdal
#endif

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: num_bins
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

!
! !DESCRIPTION:
!  This subroutine reads the planting day map into the input file.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[fgrd]
!     output field with the planting day
!   \end{description}
!EOP
   integer :: water_class
   integer :: ftn
   real    :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real    :: gridcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
   logical :: file_exists

   logical :: file_read_status

   integer :: i, t, c, r, line
   integer :: glpnc, glpnr                          ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr                        ! Parameter subsetted columns and rows
   real      :: dres
   integer :: mi                                    ! Total number of input param grid array points
   integer :: mo                                    ! Total number of output LIS grid array points
   integer*2,allocatable  :: zval(:,:), zval2(:,:)  ! (Input grid)
   real    :: param_gridDesc(20)                    ! Input parameter grid desc fgrd
   real    :: subparam_gridDesc(20)                 ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)

#if (defined USE_GDAL)
   TYPE(gdaldriverh)      :: driver
   TYPE(gdaldataseth)     :: ds
   TYPE(gdalrasterbandh)  :: band
   INTEGER(kind=c_int)    :: xsize, ysize
   REAL(kind=c_double)    :: x1, y1, x2, y2, gt(6)
   INTEGER(kind=c_int)    :: ierr
#endif
   integer               :: x_offset, y_offset
   ! Read input parameter
   integer*2,allocatable  :: read_inputparm(:,:)
   ! input parameter 1d grid
  real,    allocatable  :: gi(:)
  ! input logical mask (to match gi)
  logical*1,allocatable :: li(:)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output logical mask (go1)
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins) ! Output lis 1d grid
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins) ! Output logical mask (go2)
!__________________________________________________________________
!- Initialize:
   temp = LDT_rc%udef
   gridcnt  = 0.
   fgrd    = 0.
   param_gridDesc = 0.

!- Check if planting day file exists:
      inquire( file=trim(LDT_plantingday_struc(n)%plantingdayfile), exist=file_exists )
      if(.not. file_exists) then
         write(LDT_logunit,*) "[ERR] Planting day map: ",trim(LDT_plantingday_struc(n)%plantingdayfile)," does not exist. "
         write(LDT_logunit,*) "program stopping ..."
         call LDT_endrun
      endif
      write(LDT_logunit,*)"[INFO] Reading Palnting day map: ",trim(LDT_plantingday_struc(n)%plantingdayfile)

#if (defined USE_GDAL)
   file_read_status = .false.
   call GDALAllRegister()

   ! Use GDAL routines to open the tiff files:
   driver = gdalgetdriverbyname('Tif'//CHAR(0))
   ds = gdalopen( trim(LDT_plantingday_struc(n)%plantingdayfile)//CHAR(0), GA_ReadOnly)

   if( .not.gdalassociated(ds) ) then
      write(LDT_logunit,*) "[ERR] Opening dataset on file ",&
           trim(LDT_plantingday_struc(n)%plantingdayfile)," failed ..."
      call LDT_endrun
   end if

   ! Determine Tiff-file grid information:
   ierr = gdalgetgeotransform(ds, gt)
   xsize = gdalgetrasterxsize(ds)
   ysize = gdalgetrasterysize(ds)
   CALL gdalapplygeotransform(gt, 0.5_c_double, 0.5_c_double, x2, y1)
   CALL gdalapplygeotransform(gt, real(xsize-0.5_c_double,KIND=C_DOUBLE), real(ysize-0.5_c_double,KIND=C_DOUBLE), x1, y2)
   dres = abs(gt(6))

!- Set parameter grid array inputs:
   LDT_plantingday_struc(n)%plantingday_proj = "latlon"
   param_gridDesc(1)  = 0.        ! Latlon
   param_gridDesc(2)  = xsize
   param_gridDesc(3)  = ysize
   param_gridDesc(4)  = y2        ! LL lat
   param_gridDesc(5)  = x2        ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = y1        !UR lat
   param_gridDesc(8)  = x1        !UR lon
   param_gridDesc(9)  = dres
   param_gridDesc(10) = dres
   param_gridDesc(20) = 64

   !counter = 0 !Used when mutiple layers
   !counter=counter+1

   band = gdalgetrasterband(ds, 1)
   if (.NOT.gdalassociated(band)) THEN
      write(LDT_logunit,*) '[ERR] Failed getting raster band from GeoTIFF dataset on file, ',&
           TRIM(LDT_plantingday_struc(n)%plantingdayfile)
      call LDT_endrun()
   endif

   allocate(zval(xsize, ysize))
   ierr = gdalrasterio_f(band, GF_Read, 0, 0, zval)
   if (ierr /= 0) THEN
      write(LDT_logunit,*) '[ERR] Reading data from GeoTIFF dataset on file ',&
           TRIM(LDT_plantingday_struc(n)%plantingdayfile)
      call LDT_endrun()
   endif
   call gdalclose(ds)

   allocate(zval2(xsize,ysize)) !Rearrange the values's direction
   do r=1,ysize
     do c=1,xsize
        zval2(c,r) = zval(c,ysize-r+1)
     enddo
   enddo
   deallocate(zval)

!- Map Parameter Grid Info to LIS Target Grid/Projection Info --
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_plantingday_struc(n)%plantingday_proj, param_gridDesc, &
                     glpnc, glpnr, subpnc, subpnr,  &
                     subparam_gridDesc, lat_line, lon_line )

! Initialize parameter read-in array:
   allocate( read_inputparm(subpnc, subpnr) )
   read_inputparm = LDT_rc%udef

   x_offset = nint((subparam_griddesc(5)-param_gridDesc(5))/&
        param_gridDesc(9)) + 1
   y_offset = nint((subparam_griddesc(4)-param_gridDesc(4))/&
        param_gridDesc(10)) + 1

   read_inputparm = zval2(x_offset:x_offset+subpnc, &
                          y_offset:y_offset+subpnr)

   deallocate(zval2)
   deallocate(lat_line)
   deallocate(lon_line)

!! -------------------------------------------------------------------
!!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
!! -------------------------------------------------------------------
   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

   if( mi .ne. mo .and. LDT_plantingday_struc(n)%plantingday_gridtransform == "none" ) then
      write(LDT_logunit,*) "[ERR] Spatial transform, 'none', is selected, but number of"
      write(LDT_logunit,*) "  input and output points do not match. Select other spatial"
      write(LDT_logunit,*) "  option (e.g., mode, neighbor, tile, etc.)."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   allocate( gi(mi), li(mi) )
   gi = 0.0
   li = .true.
   lo1 = .false.;  lo2 = .false.

   ! Assign 2-D array to 1-D for grid transformation routines:
   i = 0
   do r = 1, subpnr
      do c = 1, subpnc;  i = i + 1
         gi(i) = read_inputparm(c,r)
      enddo
   enddo
   deallocate(read_inputparm)

! Apply the spatial transform option:
   select case( LDT_plantingday_struc(n)%plantingday_gridtransform)

   ! (1) Single-layer selection:
      case("none", "mode", "neighbor" )

     ! Transform parameter from original grid to LIS output grid:
       call LDT_transform_paramgrid(n, LDT_plantingday_struc(n)%plantingday_gridtransform, &
                          subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

!
!    !- Aggregate:
!
      !- Convert 1D dominant planting day to 2D grid arrays:
          i = 0
          do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                i = i + 1
                temp(c,r) = go1(i)
             enddo
          enddo

       case( "tile" )

       ! Calculate total counts for each planting day in each coarse gridcell
       call LDT_transform_paramgrid( n, LDT_plantingday_struc(n)%plantingday_gridtransform, &
            subparam_gridDesc, mi, num_bins, gi, li, mo, go2, lo2 )

       ! Convert 1D to 2D grid arrays:
           i = 0
           do r = 1, LDT_rc%lnr(n)
              do c = 1, LDT_rc%lnc(n);  i = i + 1
                 do t = 1, num_bins
                       gridcnt(c,r,t) = go2(i,t)
                 end do
              enddo
           enddo

     case default
       write(LDT_logunit,*) "[ERR] Selected grid transform option, "//&
                 LDT_plantingday_struc(n)%plantingday_gridtransform
       write(LDT_logunit,*) "[ERR]  not supported for the planting day reader."
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun

   end select
   deallocate( gi, li )

!! ........................................................................

!- Bring 2-D planting day to 3-D planting day tile space:
   if ( LDT_plantingday_struc(n)%plantingday_gridtransform == "none" .or. &
        LDT_plantingday_struc(n)%plantingday_gridtransform == "mode" .or. &
        LDT_plantingday_struc(n)%plantingday_gridtransform == "neighbor") then  ! -- NON-TILED SURFACES
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if ( temp(c,r) .le. 0 ) then
                temp(c,r) = NINT(temp(c,r))
           endif
           if ( (nint(temp(c,r)) .ne. LDT_rc%udef)) then
                gridcnt(c,r,1) = NINT(temp(c,r))
           endif
        enddo
     end do
   endif   ! End NON-TILED planting day option


  fgrd = gridcnt  ! assign the planting day grid to fraction grid (one dimension only)

   write(LDT_logunit,*) "[INFO] Done reading Planting Day "

#endif
end subroutine read_PlantingDay
