module LIS_create_gridMod

   use ESMF
   use LIS_logMod,     only : LIS_verify 
   use LIS_coreMod !,    only : LIS_masterproc, LIS_localPet
   !use LIS_coreMod,  only : vm=>LIS_vm, LIS_rc, LIS_Domain
   
   implicit none
   private
   public  :: create_regular_grid
   
!EOP
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!BOP
   function create_regular_grid(vm, grid_desc, grid_name, Nx, Ny, staggerloc) result(new_grid)
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

      rc = ESMF_SUCCESS

      !if (LIS_masterproc) print *, "Create the ESMF Grid: ", TRIM(grid_name)

      numGridPoints(1) = int(grid_desc( 2))   ! Along longitudes
      numGridPoints(2) = int(grid_desc( 3))   ! Along latitudes

      dx = grid_desc( 9) ! (grid_desc(8) - grid_desc(5))/(numGridPoints(1)-1)
      dy = grid_desc(10) ! (grid_desc(7) - grid_desc(4))/(numGridPoints(2)-1)

      if (LIS_masterproc) print *, "KNJR->"//TRIM(grid_name)//"dx/dy:", dx, dy, grid_desc( 9),grid_desc(10)
      !  -------------------------------------------------------
      new_grid = ESMF_GridCreateNoPeriDim(minIndex       = (/1,1/), &
                                          maxIndex       = numGridPoints, &
                                          gridEdgeLWidth = (/0,0/), &
                                          gridEdgeUWidth = (/0,0/), &
                                          coordSys=ESMF_COORDSYS_CART, & ! Cartesian coordinates (x,y not lat,lon)
                                          indexflag      = ESMF_INDEX_GLOBAL, &
                                          regDecomp      = (/NX, NY /), &
                                          name           = TRIM(grid_name), rc=rc)

!      new_grid = ESMF_GridCreate( name           = TRIM(grid_name),     &
!                                  regDecomp      = (/NX, NY /),         &
!                                  minIndex       = (/1,1/),             &
!                                  maxIndex       = numGridPoints,       &
!                                  indexflag      = ESMF_INDEX_GLOBAL,   &
!                                  gridMemLBound  = (/1,1/),             &
!                                  gridEdgeLWidth = (/0,0/),             &
!                                  gridEdgeUWidth = (/0,0/),             &
!                                  coordDep1      = (/1,2/),             &
!                                  coordDep2      = (/1,2/),             &
!                                  rc=rc)

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
      if (LIS_masterproc) print *, "Done creating the ESMF Grid: ", TRIM(grid_name)

   end function create_regular_grid
!EOC
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
      subroutine decomposeDim ( dim_world,dim,NDEs )
      implicit   none
      integer    dim_world, NDEs
      integer    dim(0:NDEs-1)
      integer    n,im,rm,nbeg,nend
      im = dim_world/NDEs
      rm = dim_world-NDEs*im
      do n=0,NDEs-1
                      dim(n) = im
      if( n.le.rm-1 ) dim(n) = im+1
      enddo
      end subroutine decomposeDim

end module LIS_create_gridMod


