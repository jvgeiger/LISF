! https://github.com/geoschem/gchp/tree/master/ESMF/src/system_tests
module LIS_create_gridMod

   use ESMF
   use LIS_logMod,     only : LIS_verify, LIS_logunit 
   use LIS_coreMod !,    only : LIS_masterproc, LIS_localPet
   !use LIS_coreMod,  only : vm=>LIS_vm, LIS_rc, LIS_Domain
   
   implicit none
   private
   public  :: create_rectilinear_grid
   public  :: create_regular_grid
   public  :: create_regular_grid_2D
   
!EOP
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!BOP
   function create_rectilinear_grid(lon_points, lat_points, grid_name, Nx, Ny, &
                                 coordSys, periodic, staggerloc) result(new_grid)
      integer,          intent(in)    :: Nx ! number of PEs along longitudes
      integer,          intent(in)    :: Ny ! number of PEs along latitudes
      character(len=*), intent(in)    :: grid_name
      real,             intent(in)    :: lat_points(:)
      real,             intent(in)    :: lon_points(:)
      type(ESMF_CoordSys_Flag)        :: coordSys ! Coordinate system: cartesian, spherical
      logical,               optional :: periodic
      type(ESMF_STAGGERLOC), optional :: staggerloc
!
! !RETURNED VALUE:
      type(ESMF_Grid)                 :: new_grid
!
! !DESCRIPTION:
! Create a generic ESMF 2D rectilinear grid given arbitrary latitude and
! longitude grid points.
!
! !LOCAL VARIABLES:
      integer                         :: i, j, num_lats, num_lons
      integer                         :: tlb(1), tub(1)
      integer                         :: rc, STATUS
      real(ESMF_KIND_R8), pointer     :: coordX(:)
      real(ESMF_KIND_R8), pointer     :: coordY(:)
      type(ESMF_STAGGERLOC)           :: staggerloc_
      logical                         :: periodic_
      integer                         :: countsPerDEDimX(Nx)
      integer                         :: countsPerDEDimY(Ny)
      CHARACTER(len=ESMF_MAXSTR)      :: IAm = "create_rectilinear_grid"
!EOP
!------------------------------------------------------------------------------
!BOC
      rc = ESMF_SUCCESS

      write(LIS_logunit,*) '[INFO] Create the ESMF Rectilinear Grid: '// TRIM(grid_name)

      ! Determine the number of grid points to distributed to processors
      !-----------------------------------------------------------------
      num_lons = SIZE(lon_points)
      call decomposeDim(num_lons, countsPerDEDimX, Nx )

      num_lats = SIZE(lat_points)
      call decomposeDim(num_lats, countsPerDEDimY, Ny )

      ! Create the ESMF Grid
      !---------------------
      periodic_ = .FALSE.
      if (PRESENT(periodic)) periodic_ = periodic
      
      if (periodic_) then
         new_grid = ESMF_GridCreate1PeriDim( &
                          ! Define an irregular distribution
                          countsPerDEDim1 = countsPerDEDimX, &
                          countsPerDEDim2 = countsPerDEDimY, &
                          ! Specify mapping of coords dim to Grid dim
                          coordDep1       = (/1/), & ! 1st coord is 1D and depends on 1st Grid dim
                          coordDep2       = (/2/), & ! 2nd coord is 2D and depends on 2nd Grid dim
                          indexflag       = ESMF_INDEX_GLOBAL, &
                          coordSys        = coordSys, &
                          name            = TRIM(grid_name), rc=rc)

         call LIS_verify(rc, 'ESMF_GridCreate1PeriDim failed in '//TRIM(IAm))
      else
         new_grid = ESMF_GridCreateNoPeriDim( &
                          ! Define an irregular distribution
                          countsPerDEDim1 = countsPerDEDimX, &
                          countsPerDEDim2 = countsPerDEDimY, &
                          ! Specify mapping of coords dim to Grid dim
                          coordDep1       = (/1/), & ! 1st coord is 1D and depends on 1st Grid dim
                          coordDep2       = (/2/), & ! 2nd coord is 2D and depends on 2nd Grid dim
                          indexflag       = ESMF_INDEX_GLOBAL, &
                          coordSys        = coordSys, & 
                          name            = TRIM(grid_name), rc=rc)

         call LIS_verify(rc, 'ESMF_GridCreateNoPeriDim failed in '//TRIM(IAm))
      endif

      !-------------------------------------------------------------------
      ! Allocate coordinate storage and associate it with the center
      ! stagger location.  Since no coordinate values are specified in
      ! this call no coordinate values are set yet.
      !-------------------------------------------------------------------
      staggerloc_ = ESMF_STAGGERLOC_CENTER
      if (PRESENT(staggerloc)) staggerloc_ = staggerloc

      call ESMF_GridAddCoord(new_grid, staggerloc=staggerloc_, rc=rc)
      call LIS_verify(rc, 'ESMF_GridAddCoord failed in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Get the pointer to the first coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 1, &
                                       farrayPtr           = coordX, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for x-axis in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Calculate and set coordinates in the first dimension.
      !-------------------------------------------------------------------
      do i = tlb(1), tub(1)
         coordX(i) = lon_points(i)    ! longitude grid points
      enddo

      !-------------------------------------------------------------------
      ! Get the pointer to the second coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 2, &
                                       farrayPtr           = coordY, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for y-axis in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Calculate and set coordinates in the second dimension.
      !-------------------------------------------------------------------
      do j = tlb(1), tub(1)
         coordY(j) = lat_points(j)    ! latitude  grid points
      enddo

      write(LIS_logunit,*) '[INFO] Done creating the ESMF Rectilinear Grid: '// TRIM(grid_name)

   end function create_rectilinear_grid
!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
   function create_regular_grid_2D(grid_desc, grid_name, Nx, Ny, staggerloc) result(new_grid)
      integer, intent(in) :: Nx, Ny
      character(len=*) :: grid_name
      real,    intent(in) :: grid_desc(:)
      type(ESMF_STAGGERLOC), optional :: staggerloc

      type(ESMF_Grid)                 :: new_grid
      integer                         :: i, j, i1, i2, j1, j2
      integer                         :: numGridPoints(2), tlb(2), tub(2)
      integer                         ::  rc, STATUS
      real(ESMF_KIND_R8)              ::  dx, dy
      real(ESMF_KIND_R8), pointer     :: coordX(:,:)
      real(ESMF_KIND_R8), pointer     :: coordY(:,:)
      type(ESMF_STAGGERLOC)           :: staggerloc_

      real, parameter                 :: D2R = 4.0*ATAN(1.0) / 180.

      rc = ESMF_SUCCESS

      write(LIS_logunit,*) '[INFO] Create the ESMF lat/lon Grid: '// TRIM(grid_name)

      !if (LIS_masterproc) print *, "Create the ESMF Grid: ", TRIM(grid_name)

      numGridPoints(1) = int(grid_desc( 2))   ! Along longitudes
      numGridPoints(2) = int(grid_desc( 3))   ! Along latitudes

      dx = grid_desc( 9) ! (grid_desc(8) - grid_desc(5))/(numGridPoints(1)-1)
      dy = grid_desc(10) ! (grid_desc(7) - grid_desc(4))/(numGridPoints(2)-1)

      !if (LIS_masterproc) print *, "KNJR->"//TRIM(grid_name)//"dx/dy:", dx, dy, grid_desc( 9),grid_desc(10)
      !  -------------------------------------------------------
      new_grid = ESMF_GridCreateNoPeriDim(minIndex       = (/1,1/), &
                                          maxIndex       = numGridPoints, &
                                          gridEdgeLWidth = (/0,0/), &
                                          gridEdgeUWidth = (/0,0/), &
                                          coordSys=ESMF_COORDSYS_CART, & ! Cartesian coordinates (x,y not lat,lon)
                                          indexflag      = ESMF_INDEX_GLOBAL, &
                                          regDecomp      = (/NX, NY /), &
                                          name           = TRIM(grid_name), rc=rc)

	call LIS_verify(rc, 'ESMF_GridCreateNoPeriDim failed')

      staggerloc_ = ESMF_STAGGERLOC_CENTER
      if (PRESENT(staggerloc)) staggerloc_ = staggerloc

      call ESMF_GridAddCoord(new_grid, staggerloc=staggerloc_, rc=rc)
      call LIS_verify(rc, 'ESMF_GridAddCoord failed')

      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 1, &
                                       farrayPtr           = coordX, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for x-axis')

      call ESMF_GridGetCoord(new_grid, localDE    = 0, &
                                       staggerLoc = staggerloc_, &
                                       coordDim   = 2, &
                                       farrayPtr  = coordY, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for y-axis')

      do j = tlb(2), tub(2)
         do i = tlb(1), tub(1)
            coordX(i,j) = (i-1)*dx  + grid_desc(5)    ! longitude grid points
            coordY(i,j) = (j-1)*dy  + grid_desc(4)    ! latitude  grid points
         enddo
      enddo
      !if (LIS_masterproc) print *, "Done creating the ESMF Grid: ", TRIM(grid_name)
      write(LIS_logunit,*) '[INFO] Done creating the ESMF lat/lon Grid: '// TRIM(grid_name)

   end function create_regular_grid_2D
!EOC
!------------------------------------------------------------------------------
!BOP
   function create_regular_grid(grid_desc, grid_name, Nx, Ny, staggerloc) result(new_grid)
      integer,          intent(in)    :: Nx ! number of PEs along longitudes
      integer,          intent(in)    :: Ny ! number of PEs along latitudes
      character(len=*), intent(in)    :: grid_name
      real,             intent(in)    :: grid_desc(:)
      type(ESMF_STAGGERLOC), optional :: staggerloc
!
! !RETURNED VALUE:
      type(ESMF_Grid)                 :: new_grid
!
! !DESCRIPTION:
! Create an ESMF 2D lat/lon grid
!
! !LOCAL VARIABLES:
      integer                         :: i, j, i1, i2, j1, j2
      integer                         :: numGridPoints(2)
      integer                         :: tlb(1), tub(1)
      integer                         ::  rc, STATUS
      real(ESMF_KIND_R8)              ::  dx, dy
      real(ESMF_KIND_R8), pointer     :: coordX(:)
      real(ESMF_KIND_R8), pointer     :: coordY(:)
      type(ESMF_STAGGERLOC)           :: staggerloc_
      integer                         :: countsPerDEDimX(Nx)
      integer                         :: countsPerDEDimY(Ny)
      CHARACTER(len=ESMF_MAXSTR)      :: IAm = 'create_regular_grid'
!EOP
!------------------------------------------------------------------------------
!BOC
      rc = ESMF_SUCCESS

      write(LIS_logunit,*) '[INFO] Create the ESMF lat/lon Grid: '// TRIM(grid_name)

      numGridPoints(1) = int(grid_desc(2))   ! Along longitudes
      numGridPoints(2) = int(grid_desc(3))   ! Along latitudes

      dx = grid_desc( 9) 
      dy = grid_desc(10) 

      ! Determine the number of grid points to distributed to processors
      !-----------------------------------------------------------------
      !ALLOCATE(countsPerDEDimX(Nx))
      call decomposeDim(numGridPoints(1), countsPerDEDimX, Nx )

      !ALLOCATE(countsPerDEDimY(Ny))
      call decomposeDim(numGridPoints(2), countsPerDEDimY, Ny )

      ! Create the ESMF Grid
      !---------------------
      new_grid = ESMF_GridCreateNoPeriDim( &
                          ! Define an irregular distribution
                          countsPerDEDim1 = countsPerDEDimX, &
                          countsPerDEDim2 = countsPerDEDimY, &
                          ! Specify mapping of coords dim to Grid dim
                          coordDep1       = (/1/), & ! 1st coord is 1D and depends on 1st Grid dim
                          coordDep2       = (/2/), & ! 2nd coord is 2D and depends on 2nd Grid dim
                          indexflag       = ESMF_INDEX_GLOBAL, &
                          coordSys        = ESMF_COORDSYS_CART, & ! use cartesian coordinates
                          name            = TRIM(grid_name), rc=rc)
      call LIS_verify(rc, 'ESMF_GridCreateNoPeriDim failed in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Allocate coordinate storage and associate it with the center
      ! stagger location.  Since no coordinate values are specified in
      ! this call no coordinate values are set yet.
      !-------------------------------------------------------------------
      staggerloc_ = ESMF_STAGGERLOC_CENTER
      if (PRESENT(staggerloc)) staggerloc_ = staggerloc

      call ESMF_GridAddCoord(new_grid, staggerloc=staggerloc_, rc=rc)
      call LIS_verify(rc, 'ESMF_GridAddCoord failed in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Get the pointer to the first coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 1, &
                                       farrayPtr           = coordX, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for x-axis in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Calculate and set coordinates in the first dimension.
      !-------------------------------------------------------------------
      do i = tlb(1), tub(1)
         coordX(i) = (i-1)*dx  + grid_desc(5)    ! longitude grid points
      enddo

      !-------------------------------------------------------------------
      ! Get the pointer to the second coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 2, &
                                       farrayPtr           = coordY, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for y-axis in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Calculate and set coordinates in the second dimension.
      !-------------------------------------------------------------------
      do j = tlb(1), tub(1)
         coordY(j) = (j-1)*dy  + grid_desc(4)    ! latitude  grid points
      enddo

      write(LIS_logunit,*) '[INFO] Done creating the ESMF lat/lon Grid: '// TRIM(grid_name)

   end function create_regular_grid
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine decomposeDim(dim_world, dim_array, NDEs )
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: dim_world ! total number of grid points
      integer, intent(in)  :: NDEs      ! number of DEs
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: dim_array(0:NDEs-1) 
!
! !DESCRIPTION:
! For a given number of grid points and a number of available processors,
! this subroutines determines the number of grid points assigned to each
! processor.
!
! !LOCAL VARIABLES:
      integer ::   n, im, rm
!EOP
!------------------------------------------------------------------------------
!BOC
      im = dim_world/NDEs
      rm = dim_world-NDEs*im
      do n=0,NDEs-1
                      dim_array(n) = im
      if( n.le.rm-1 ) dim_array(n) = im+1
      enddo
      end subroutine decomposeDim
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: gaussian_find_row
! \label{gaussian_find_row}
!
! !INTERFACE:
function gaussian_find_row(jmax, gaussian_lat_array, lat)
! !USES:
   use LIS_logMod , only : LIS_logunit, LIS_endrun

   implicit none
! !ARGUMENTS:
   integer, intent(in)               :: jmax
   real, dimension(jmax), intent(in) :: gaussian_lat_array
   real, intent(in)                  :: lat
!
! !DESCRIPTION: 
!  This function computes the row number within the global quasi-regular 
!  Gaussian grid corresponding to the latitude given by lat.
!
!   \begin{description}
!   \item [jmax]
!      the total number of latitude circles
!   \item [gaussian\_lat\_array]
!      array of computed latitudes, in degrees
!   \item [lat]
!      latitude to search for
!   \end{description}
!EOP

   integer :: gaussian_find_row

   real    :: eps
   integer :: r

   eps = 180./(2*jmax)

   gaussian_find_row = -9999
   do r = 1, jmax
      if ( abs(gaussian_lat_array(r) - lat) < eps ) then
      !if ( gaussian_latitudes(1,r) == lat ) then
         gaussian_find_row = r
      endif
   enddo

   if ( gaussian_find_row == -9999 ) then
      write(LIS_logunit,fmt=*) '[ERR] gaussian_find_row -- '// &
                               'Could not find lat',lat
      call LIS_endrun
   endif

end function gaussian_find_row
!------------------------------------------------------------------------------

end module LIS_create_gridMod
