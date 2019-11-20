!
! Example/test code which shows User Component calls.

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!BOP
!
! !DESCRIPTION:
!  User-supplied Component, most recent interface revision.
!
!
!\begin{verbatim}

    module LIS_metforcing_componentMod

      ! ESMF Framework module
      use ESMF
      use LIS_coreMod,              only : LIS_rc, LIS_vm, LIS_masterproc, LIS_localPet
      use LIS_logMod,               only : LIS_verify
      use LIS_FORC_AttributesMod
      use nldas2_forcingMod,        only : nldas2_struc
      use LIS_field_bundleMod
      use LIS_create_gridMod,       only : create_regular_grid
      USE LIS_metforcing_nldas2Mod

      USE LIS_timeMgrMod,           ONLY : LIS_date2time

      IMPLICIT NONE
    
      public met_forcing_register
      public forcing_bundle_bname
      public  :: nldas2_file_info
      public  :: set_nldas2_file_name
      public setFor_forcast_member_index
      public setFor_nest_index


      real, save :: forcing_gridDesc(10)
      character(len=20), parameter ::                  Iam = "Forcing Data Comp:"
      character(len=20), parameter :: forcing_bundle_bname = "Forcing_Data_Bundle_"

      type(ESMF_STAGGERLOC), save :: staggerloc

      type(ESMF_TypeKind_Flag) :: type_kind = ESMF_TYPEKIND_R4
      integer, allocatable :: lbnd(:,:), ubnd(:,:)

      type nldas2_file_info
          character(len=100) ::  file_nameA
          character(len=100) ::  file_nameB
          character(len=  1) ::  file_type   ! 'A' or 'B'
      end type nldas2_file_info
 
      type(nldas2_file_info), save :: nldas2_file_name

      integer :: forcast_member_index
      integer :: nest_index

        
!-------------------------------------------------------------------------
    contains
!-------------------------------------------------------------------------
!   !  The Register routine sets the subroutines to be called
!   !   as the init, run, and finalize routines.  Note that these are
!   !   private to the module.
 
    subroutine met_forcing_register(comp, rc)
        type(ESMF_GridComp) :: comp
        integer, intent(out) :: rc

        if (LIS_masterproc) print *, TRIM(Iam) // ": Register routine"

        ! Register the callback routines.

        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, met_forcing_init, rc=rc)
        call LIS_verify(rc, 'ESMF_GridAddCoord failed')

        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, met_forcing_run, rc=rc)
        call LIS_verify(rc, 'ESMF_GridCompSetEntryPoint failed for run')

        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, met_forcing_final, rc=rc)
        call LIS_verify(rc, 'ESMF_GridCompSetEntryPoint failed for finalize')

        if (LIS_masterproc) &
           print *, TRIM(Iam) // ": Registered Initialize, Run, and Finalize routines"

    end subroutine met_forcing_register

!-------------------------------------------------------------------------
!   !  User Comp Component created by higher level calls, here is the
!   !   Initialization routine.
 
    
    subroutine met_forcing_init(comp, importState, exportState, clock, rc)
        type(ESMF_GridComp) :: comp
        type(ESMF_State) :: importState, exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc

       ! Local variables
        type(ESMF_Field)                    :: dummy
        type(ESMF_FieldBundle), allocatable :: met_forcing_bundle(:)
        type(ESMF_Grid),        allocatable :: esmfGrid(:)
        type(ESMF_ArraySpec)                :: arrayspec
        real(ESMF_KIND_R8),         pointer :: ptr2D(:,:)
        real(ESMF_KIND_R8),         pointer :: coordX(:,:), coordY(:,:)
        integer                             :: i, j
        real(ESMF_KIND_R8)                  :: min(2), max(2), dx, dy
        type(ESMF_Config)  :: config
        integer :: fieldcount
        integer :: nx, ny, N, iv
        character(len=50) :: staggerloc_var
        character(len=100) :: var_name
        real(kind=8) :: temp_var
        character(len=2) :: num_st
        real                                :: dummy_array(1,1) = 0.0

!EOP
!------------------------------------------------------------------------------
!BOC
        if (LIS_masterproc) print *, TRIM(Iam) // ": Init starting"
   
        forcing_gridDesc(:) = 0.0

        Nx = LIS_rc%npesx ! Number of processors along x
        Ny = LIS_rc%npesy ! Number of processors along y

        staggerloc_var = "center"

        IF (TRIM(staggerloc_var) == 'center') THEN
           staggerloc = ESMF_STAGGERLOC_CENTER
        ELSEIF (TRIM(staggerloc_var) == 'edge1') THEN
           staggerloc = ESMF_STAGGERLOC_EDGE1
        ELSEIF (TRIM(staggerloc_var) == 'edge2') THEN
           staggerloc = ESMF_STAGGERLOC_EDGE2
        ENDIF

        allocate(met_forcing_bundle(LIS_rc%nnest))
        allocate(esmfGrid(LIS_rc%nnest))
        allocate(lbnd(LIS_rc%nnest, 2))
        allocate(ubnd(LIS_rc%nnest, 2))

        if (LIS_masterproc) then
           print*, "<>------ Forcing data Grid -------<>"
           print*, "Number of nest:     ", LIS_rc%nnest
           print*, "Nx:                 ", Nx
           print*, "Ny:                 ", Ny
           print*, "num points along x: ", nldas2_struc(1)%gridDesc(2)
           print*, "num points along y: ", nldas2_struc(1)%gridDesc(3)
           print*, "lower lat:          ", nldas2_struc(1)%gridDesc(4)
           print*, "lower lon:          ", nldas2_struc(1)%gridDesc(5)
           print*, "upper lat:          ", nldas2_struc(1)%gridDesc(7)
           print*, "upper lon:          ", nldas2_struc(1)%gridDesc(8)
           print*, "<>------ Forcing data Grid -------<>"
        endif

        ! Set the list of fields to be regridded
        CALL set_list_nldas2_fields()

        do n = 1, LIS_rc%nnest
           forcing_gridDesc(2) = nldas2_struc(n)%gridDesc(2) ! num points along x
           forcing_gridDesc(3) = nldas2_struc(n)%gridDesc(3) ! num points along y
           forcing_gridDesc(4) = nldas2_struc(n)%gridDesc(4) ! lower lat
           forcing_gridDesc(5) = nldas2_struc(n)%gridDesc(5) ! lower lon
           forcing_gridDesc(7) = nldas2_struc(n)%gridDesc(7) ! upper lat
           forcing_gridDesc(8) = nldas2_struc(n)%gridDesc(8) ! upper lon
           forcing_gridDesc(9) = nldas2_struc(n)%gridDesc(9) ! x-grid size
           forcing_gridDesc(10)= nldas2_struc(n)%gridDesc(10) ! y-grid size

           ! Create the bundle
           write(num_st, '(i2.2)') n
           met_forcing_bundle(n) = ESMF_FieldBundleCreate(name = forcing_bundle_bname//num_st, rc=rc)
           call LIS_verify(rc, 'ESMF_FieldBundleCreate failed for Met Forcing Data')
        
           ! Create the grid
           esmfGrid(n) = create_regular_grid(LIS_vm, forcing_gridDesc, "Forcing Grid", &
                             Nx, Ny, staggerloc = staggerloc)

           ! Get the computational bounds
           call ESMF_GridGetCoord(esmfGrid(n),                           &
                                  localDE             = 0,               &
                                  staggerLoc          = staggerloc,      &
                                  coordDim            = 1,               &
                                  farrayPtr           = coordX,          &
                                  computationalLBound = lbnd(n,:),       &
                                  computationalUBound = ubnd(n,:),  rc=rc)
           call LIS_verify(rc, 'ESMF_GridGetCoord failed for Met Forcing Data')

           !print*,"Grid 2:",LIS_localPet,lbnd,ubnd

           ! Add fields to the bundle
           DO iv = 1, number_nldas2_fields
              call addTracerToBundle(met_forcing_bundle(n), esmfGrid(n), & 
                                     TRIM(list_nldas2_fields(iv)), type_kind, dummy_array)
           ENDDO

           ! Add bundle to the state
           call ESMF_StateAdd(exportState, (/met_forcing_bundle(n)/), rc=rc)
           call LIS_verify(rc, 'ESMF_StateAddBundle failed for bundle')
        enddo

        if (LIS_masterproc) print *, ": Forcing Comp Init returning"

    end subroutine met_forcing_init
!EOC
!-------------------------------------------------------------------------
!   !  The Run routine where data is computed.
!   !
 
    subroutine met_forcing_run(comp, importState, exportState, clock, rc)
        type(ESMF_GridComp) :: comp
        type(ESMF_State) :: importState, exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc

       ! Local variables
        type(ESMF_Field)                    :: returnedfield
        type(ESMF_ARRAY)                    :: array
        type(ESMF_FieldBundle), allocatable :: met_forcing_bundle(:)
        type(ESMF_grid),        allocatable :: esmfGrid(:)
        real(ESMF_KIND_R8),         pointer :: ptr2D(:,:), ptr(:,:)
        integer                             :: n, fieldcount
        character(len=2)                    :: num_st
        real(ESMF_KIND_R8), dimension(:,:), pointer :: coordX, coordY
        integer :: i, j, imin, imax, jmin, jmax, kk
        character(len=100) :: var_name, file_name
        REAL(kind=8)       :: dummy_time
        REAL(kind=4)       :: dummy_gmt
        INTEGER            :: curDayOfYear

        rc = ESMF_SUCCESS
        if (LIS_masterproc) print *, TRIM(Iam) // ": Run starting"

        allocate(met_forcing_bundle(LIS_rc%nnest))
        allocate(esmfGrid(LIS_rc%nnest))

        kk = forcast_member_index
        n  = nest_index

        !do n = 1, LIS_rc%nnest

           ! Get the FieldBundle data from the State
           write(num_st, '(i2.2)') n
           call ESMF_StateGet(exportState, forcing_bundle_bname//num_st, met_forcing_bundle(n), rc=rc)

           ! Get date/Time data
!           CALL LIS_date2time( dummy_time, curDayOfYear, dummy_gmt, &
!                               LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
!                               LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
!
!           ! Get the file name for reading
!           call get_nldas2_grib_file(file_name, "A", n, LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
!                                      curDayOfYear, LIS_rc%hr)

           CALL get_nldas2_file_name(file_name, 'A')

           if (LIS_masterproc) print *, "ForcingFile: "//TRIM(file_name)

           ! Read the grib file and populate the bundle
           CALL read_nldas2a_grib_file(file_name, met_forcing_bundle(n), n, lbnd(n,:), ubnd(n,:))

        !enddo

      if (LIS_masterproc) print *, TRIM(Iam) // ": Run ending"

    end subroutine met_forcing_run
!EOC
!-------------------------------------------------------------------------
!   !  The Finalization routine where things are deleted and cleaned up.
!   !
 
    subroutine met_forcing_final(comp, importState, exportState, clock, rc)
        type(ESMF_GridComp) :: comp
        type(ESMF_State) :: importState, exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc

        ! Local variables
        type(ESMF_FieldBundle), allocatable :: met_forcing_bundle(:)
        type(ESMF_grid),        allocatable :: esmfGrid(:)
        integer                             :: n
        character(len=2)                    :: num_st

        rc = ESMF_SUCCESS
        if (LIS_masterproc) print *, TRIM(Iam) // ": Final starting"

        ! garbage collection   
        allocate(met_forcing_bundle(LIS_rc%nnest))

        do n = 1, LIS_rc%nnest

           ! Get the FieldBundle data from the State
           write(num_st, '(i2.2)') n
           call ESMF_StateGet(exportState, forcing_bundle_bname//num_st, met_forcing_bundle(n), rc=rc)

           call ESMF_FieldBundleDestroy(met_forcing_bundle(n), rc=rc)

        enddo

        if (LIS_masterproc) print *, TRIM(Iam) // ": Final returning"

    end subroutine met_forcing_final
!-------------------------------------------------------------------------
      subroutine get_nldas2_file_name(file_name, file_type)
          character(len=1), intent(in) ::  file_type   ! 'A' or 'B'
          character(len=*), intent(out) ::  file_name
          IF (file_type == 'A') THEN
             file_name = nldas2_file_name%file_nameA
          ELSE
             file_name = nldas2_file_name%file_nameB
          ENDIF
      end subroutine get_nldas2_file_name
!-------------------------------------------------------------------------
      subroutine set_nldas2_file_name(file_type, n, curYear, curMonth, curDay, &
                                      curDayOfYear, curHour)
          CHARACTER(len=1), intent(in) :: file_type ! 'A' or 'B'
          INTEGER, intent(in) :: n
          INTEGER, intent(in) :: curYear      ! current year
          INTEGER, intent(in) :: curMonth     ! current month
          INTEGER, intent(in) :: curDay       ! current day
          INTEGER, intent(in) :: curDayOfYear ! current day of the year
          INTEGER, intent(in) :: curHour      ! current hour

          IF (file_type == 'A') THEN
             CALL get_nldas2_grib_file(nldas2_file_name%file_nameA, 'A', n, &
                                     curYear, curMonth, curDay, &
                                     curDayOfYear, curHour)
          ELSE
             CALL get_nldas2_grib_file(nldas2_file_name%file_nameB, 'B', n, &
                                     curYear, curMonth, curDay, &
                                     curDayOfYear, curHour)
          ENDIF
      end subroutine set_nldas2_file_name
!-------------------------------------------------------------------------

      subroutine setFor_forcast_member_index(kk)
          integer, intent(in) :: kk
          forcast_member_index = kk
      end subroutine setFor_forcast_member_index

      subroutine setFor_nest_index(n)
          integer, intent(in) :: n
          nest_index = n
      end subroutine setFor_nest_index

    end module LIS_metforcing_componentMod
    
!\end{verbatim}
    
