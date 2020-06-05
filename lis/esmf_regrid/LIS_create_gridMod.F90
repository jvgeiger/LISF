! https://github.com/geoschem/gchp/tree/master/ESMF/src/system_tests
module LIS_create_gridMod

   use ESMF
   use LIS_logMod,     only : LIS_verify, LIS_logunit 
   use LIS_coreMod !,    only : LIS_masterproc, LIS_localPet
   !use LIS_coreMod,  only : vm=>LIS_vm, LIS_rc, LIS_Domain
   
   implicit none
   private
   public  :: create_regular_grid
   public  :: create_regular_grid_2D
!   public  :: create_arakawa_grid
   public  :: create_gaussian_grid
   
!EOP
!------------------------------------------------------------------------------
contains
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

      print "(a,9i5)", "KNJR->"//TRIM(grid_name),LIS_localPet,tlb(1),tub(1),tlb(2),tub(2), &
                   LIS_ews_halo_ind(1,LIS_localPet+1), LIS_ewe_halo_ind(1,LIS_localPet+1), &
                   LIS_nss_halo_ind(1,LIS_localPet+1), LIS_nse_halo_ind(1,LIS_localPet+1)

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
      PRINT "(a,i3,4f11.4)",TRIM(grid_name)//"-lat/lon:",LIS_localPet,minval(coordX),maxval(coordX),minval(coordY),maxval(coordY)
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

      print "(a,5i5)", "EW->"//TRIM(grid_name),LIS_localPet,tlb(1),tub(1), &
                   LIS_ews_halo_ind(1,LIS_localPet+1), LIS_ewe_halo_ind(1,LIS_localPet+1)

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

      print "(a,5i5)", "NS->"//TRIM(grid_name),LIS_localPet,tlb(1),tub(1), &
                   LIS_nss_halo_ind(1,LIS_localPet+1), LIS_nse_halo_ind(1,LIS_localPet+1)

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
   function create_gaussian_grid(grid_desc, grid_name, Nx, Ny, &
                                 staggerloc, global_domain) result(new_grid)
      integer,          intent(in)    :: Nx ! number of PEs along longitudes
      integer,          intent(in)    :: Ny ! number of PEs along latitudes
      character(len=*), intent(in)    :: grid_name
      real,             intent(in)    :: grid_desc(:)
      type(ESMF_STAGGERLOC), optional :: staggerloc
      logical,               optional :: global_domain
!
! !RETURNED VALUE:
      type(ESMF_Grid)                 :: new_grid
!
! !DESCRIPTION:
! Create an ESMF 2D Gaussian grid
!
! !LOCAL VARIABLES:
      integer                         :: i, j, i1, i2, j1, j2, imin_lat, imax_lat
      integer                         :: numGridPoints(2)
      integer                         :: tlb(1), tub(1)
      integer                         ::  rc, STATUS
      real(ESMF_KIND_R8)              ::  dx, dy
      real(ESMF_KIND_R8), pointer     :: lat_grid(:)
      real(ESMF_KIND_R8), pointer     :: lon_grid(:)
      real(ESMF_KIND_R8), pointer     :: coordX(:)
      real(ESMF_KIND_R8), pointer     :: coordY(:)
      type(ESMF_STAGGERLOC)           :: staggerloc_
      logical                         :: global_domain_
      integer                         :: countsPerDEDimX(Nx)
      integer                         :: countsPerDEDimY(Ny)
      real                            :: min_lat, max_lat
      real                            :: min_lon, max_lon, dlon
      real, pointer                   :: lat_points(:)
      real, pointer                   :: lat_weights(:)
      real, parameter                 :: pi=3.14159265358979
      CHARACTER(len=ESMF_MAXSTR)      :: IAm = "create_gaussian_grid"
!EOP
!------------------------------------------------------------------------------
!BOC
      rc = ESMF_SUCCESS

      write(LIS_logunit,*) '[INFO] Create the ESMF Gaussian Grid: '// TRIM(grid_name)

      write(LIS_logunit,*) grid_desc

      !if (LIS_masterproc) print *, "Create the ESMF Grid: ", TRIM(grid_name)

      numGridPoints(1) = int(grid_desc(2))   ! Along longitudes
      numGridPoints(2) = int(grid_desc(3))   ! Along latitudes

      ! Determine the grid points along:
      !   - latitudes: roots of the Gauss-Legendre polynomial
      !   - longitudes: equally spaced
      !---------------------------------------------------------
      ALLOCATE(lat_points(numGridPoints(2)))
      ALLOCATE(lat_weights(numGridPoints(2)))
      ! grid_desc(7) --> degree South
      ! grid_desc(4) --> degree North
      min_lat = grid_desc(4)
      max_lat = grid_desc(7)
      IF (max_lat < min_lat) THEN
         min_lat = grid_desc(7)
         max_lat = grid_desc(4)
      ENDIF
         
      ! Determine the global Gaussian latitude grid points
      call gausslat(numGridPoints(2), lat_points, lat_weights)
      lat_points(:) = 180.0*asin(lat_points(:))/pi
      !call gaussian_comp_lats(numGridPoints(2), lat_points)

      ! Extract the Gaussian grid points specific to the subdomain of interest
      global_domain_ = .FALSE.
      if (PRESENT(global_domain)) global_domain_ = global_domain

      IF (global_domain_) THEN
         imin_lat = 1
         imax_lat = numGridPoints(2)
      ELSE
         imin_lat = gaussian_find_row(numGridPoints(2), lat_points, min_lat)
         imax_lat = gaussian_find_row(numGridPoints(2), lat_points, max_lat)
      ENDIF

      ALLOCATE(lat_grid(imax_lat-imin_lat+1))
      lat_grid(:) = lat_points(imax_lat:imin_lat:-1)
      !lat_grid(:) = lat_points(imin_lat:imax_lat)
      lat_grid(1) = min_lat
      lat_grid(imax_lat-imin_lat+1) = max_lat

      write(LIS_logunit,*) '[KNJR] Lat Data: '// TRIM(grid_name),&
                            imax_lat-imin_lat+1,numGridPoints(2),&
                            maxval(lat_grid), minval(lat_grid)

      write(LIS_logunit,*) lat_grid

      ! Determine the longitude grid points
      ALLOCATE(lon_grid(numGridPoints(1)))
      
      min_lon = grid_desc(5)
      max_lon = grid_desc(8)
      dlon    = grid_desc(9)
      !IF ((min_lon .GE. 0.0) .AND. (max_lon .LT. 0.0)) THEN
      !   min_lon = grid_desc(5)
      !   max_lon = 360.0 + grid_desc(8)
      !ENDIF
      DO i = 1, numGridPoints(1)
         lon_grid(i) = min_lon + (i-1)*dlon
         ! values must be between -180 to 180
         IF (lon_grid(i) > 180.0) lon_grid(i) = lon_grid(i) - 360.0
      ENDDO

      ! Determine the number of grid points to distributed to processors
      !-----------------------------------------------------------------
      call decomposeDim(numGridPoints(1), countsPerDEDimX, Nx )

      call decomposeDim(imax_lat-imin_lat+1, countsPerDEDimY, Ny )
      !call decomposeDim(numGridPoints(2), countsPerDEDimY, Ny )

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
         coordX(i) = lon_grid(i)    ! longitude grid points
      enddo
      write(LIS_logunit,*) '[<-KNJR->] Lon grid: '// TRIM(grid_name), tlb(1), tub(1)

      IF (TRIM(grid_name) .EQ. 'Model Grid') THEN
         write(LIS_logunit,*) '****Lon size:   ', tub(1)-tlb(1)+1, LIS_rc%gridDesc(1,2)
         write(LIS_logunit,*) '****Lon corners:', coordX(tlb(1)),LIS_rc%gridDesc(1,5), &
                                                  coordX(tub(1)),LIS_rc%gridDesc(1,8)
      ENDIF

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
         coordY(j) = lat_grid(j)    ! latitude  grid points
      enddo
      write(LIS_logunit,*) '[<-KNJR->] Lat grid: '// TRIM(grid_name), tlb(1), tub(1)

      IF (TRIM(grid_name) .EQ. 'Model Grid') THEN
         write(LIS_logunit,*) '****Lat size:   ', tub(1)-tlb(1)+1, LIS_rc%gridDesc(1,3)
         write(LIS_logunit,*) '****Lat corners:', coordY(tlb(1)),LIS_rc%gridDesc(1,4), &
                                                  coordY(tub(1)),LIS_rc%gridDesc(1,7)
      ENDIF

      DEALLOCATE(lat_points, lat_weights)
      DEALLOCATE(lat_grid, lon_grid)

      write(LIS_logunit,*) '[INFO] Done creating the ESMF Gaussian Grid: '// TRIM(grid_name)

   end function create_gaussian_grid
!EOC
!------------------------------------------------------------------------------
!BOP
!   function create_arakawa_grid(vm, grid_desc, grid_name, Nx, Ny, horz_stagger) result(new_grid)
!      type(ESMF_VM) :: vm
!      integer, intent(in) :: Nx, Ny
!      character(len=*) :: grid_name
!      real,    intent(in) :: grid_desc(:)
!      type(ESMF_IGridHorzStagger) :: horz_stagger ! ESMF_IGRID_HORZ_STAGGER_A, ESMF_IGRID_HORZ_STAGGER_D_NE
!
!      type(ESMF_DELayout) :: delayout
!      type(ESMF_Grid)                 :: new_grid
!      integer                         :: i, j, i1, i2, j1, j2
!      integer                         :: numGridPoints(2), tlb(2), tub(2)
!      integer                         ::  rc, STATUS
!      real(ESMF_KIND_R8)              ::  dx, dy
!      real(ESMF_KIND_R8)              :: min_coord(2)
!      real(ESMF_KIND_R8)              :: max_coord(2)
!
!      rc = ESMF_SUCCESS
!
!      !if (LIS_masterproc) print *, "Create the ESMF Grid: ", TRIM(grid_name)
!
!      ! Create a layout with the right breakdown
!      delayout = ESMF_DELayoutCreate(vm, (/ NX, NY /), rc=rc)
!      call LIS_verify(rc, 'ESMF_DELayoutCreate failed')
!
!      numGridPoints(1) = int(grid_desc( 2))   ! Along longitudes
!      numGridPoints(2) = int(grid_desc( 3))   ! Along latitudes
!
!      min_coord(1) = grid_desc(5)
!      max_coord(1) = grid_desc(8)
!      min_coord(2) = grid_desc(4)
!      max_coord(2) = grid_desc(7)
!
!      !dx = grid_desc( 9) ! (grid_desc(8) - grid_desc(5))/(numGridPoints(1)-1)
!      !dy = grid_desc(10) ! (grid_desc(7) - grid_desc(4))/(numGridPoints(2)-1)
!
!      new_grid = ESMF_IGridCreateHorzXYUni(counts               = numGridPoints, &
!                                           minGlobalCoordPerDim = min_coord, &
!                                           maxGlobalCoordPerDim = max_coord, &
!                                           horzStagger          = horz_stagger, &
!                                           name                 = TRIM(grid_name), rc=rc)
!      call LIS_verify(rc, 'ESMF_IGridCreateHorzXYUni failed')
!
!      call ESMF_IGridDistribute(new_grid, delayout=delayout, rc=rc)
!      call LIS_verify(rc, 'ESMF_IGridDistribute failed')
!
!      if (LIS_masterproc) print *, "Done creating the ESMF Grid: ", TRIM(grid_name)
!
!   end function create_arakawa_grid
!!EOC
!------------------------------------------------------------------------------
!BOP

   function create_irregular_grid(vm, grid_desc, grid_name, Nx, Ny, staggerloc) result(new_grid)
      type(ESMF_VM) :: vm
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

      integer                         :: minIndex(2)
      integer, allocatable            :: IMS(:), JMS(:)
      integer                         :: IM_WORLD, JM_WORLD

!EOp
!------------------------------------------------------------------------------
!BOC
      IM_WORLD = int(grid_desc( 2))   ! Along longitudes
      JM_WORLD = int(grid_desc( 3))   ! Along latitudes

      allocate( IMS(0:NX-1), stat=STATUS )
      allocate( JMS(0:NY-1), stat=STATUS )

     call decomposeDim ( IM_WORLD,ims, NX )
     call decomposeDim ( JM_WORLD,jms, NY )

      minIndex=(/1,1/)

      new_grid = ESMF_GridCreate( name            = TRIM(grid_name),     &
                                  minIndex        = minIndex,            &
                                  countsPerDEDim1 = ims,                 &
                                  countsPerDEDim2 = jms,                 &
                                  indexflag       = ESMF_INDEX_GLOBAL,   &
                                  gridMemLBound   = (/1,1/),             &
                                  gridEdgeLWidth  = (/0,0/),             &
                                  gridEdgeUWidth  = (/0,0/),             &
                                  coordDep1       = (/1,2/),             &
                                  coordDep2       = (/1,2/),             &
                                  rc=rc)
      call LIS_verify(rc, 'ESMF_GridCreate failed')

   end function create_irregular_grid
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
      function get_gaussian_points(a, b, n) result(gaussian_points)
          ! Function to determine the Gauss-Legendre points in the 
          ! interval [a, b].
          integer            :: n
          real               :: a  ! minimum value
          real               :: b  ! maximum value
          real(ESMF_KIND_R8) :: gaussian_points(n)
           double precision  :: xroots(n)
          integer            :: i
!EOP
!------------------------------------------------------------------------------
!BOC
          xroots = gauss_legendre_roots(n)
          gaussian_points(:) = a + 0.5d0*(xroots(:) + 1.0d0)*(b-a)

      end function get_gaussian_points
!EOC
!------------------------------------------------------------------------------
!BOP
      function gauss_legendre_roots(n) result(xroots)
!
! !INPUT PARAMETERS:
          integer, intent(in) :: n
!
! !RETURNED VALUE:
          double precision, dimension(n) :: xroots
!
! !DESCRIPTION:
!  Compute the roots of the Guass-Legendre polynomial or order n
!  between -1 and 1.
!  Function adapted from:  https://userpages.umbc.edu/~squire/download/gauleg.f90
!
! !LOCAL VARIABLES:
          integer :: i, j, m
          double precision, dimension(n) :: glweights
          double precision :: p1, p2, p3, pp, z, z1
          double precision, parameter :: eps = 3.0d-15      ! Relative precision
          double precision, parameter :: myPI = 4.0d0*ATAN(1.0d0)
!EOP
!------------------------------------------------------------------------------
!BOC 
          print*, "PI = ", myPI
          ! Roots are symmetric in the interval - so only need to find half of them
          m = (n+1)/2

          ! Loop over the desired roots
          do i = 1, m
             z = cos( myPI * (i-0.25d0) / (n+0.5d0) )
             z1 = 0.0

             ! Starting with the above approximation to the ith root,
             ! we enter the main loop of refinement by Newton's method
             do while(abs(z-z1) .gt. eps)
                p1 = 1.0d0
                p2 = 0.0d0
                ! Loop up the recurrence relation to get the Legendre polynomial evaluated at z 
                do j = 1, n
                   p3 = p2
                   p2 = p1
                   p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
                end do
                ! p1 is now the desired Legendre polynomial. 
                ! We compute pp, its derivative, by a standard relation 
                ! involving also p2, the polynomial of one lower order.      
                pp = n*(z*p1-p2)/(z*z-1.0d0)
                z1 = z
                ! Newton's method
                z = z1 - p1/pp
             end do

             ! Roots will be bewteen -1.0 & 1.0
             ! and symmetric about the origin 
             xroots(    i) = -z
             xroots(n+1-i) = +z

             ! Compute the weight and its symmetric counterpart 
             glweights(    i) = 2.0d0/((1.0d0-z*z)*pp*pp)
             glweights(n+1-i) = glweights(i)
          end do
    end function gauss_legendre_roots
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
