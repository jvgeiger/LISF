!
! Example/test code which shows User Component calls.

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

!BOP
!
! !DESCRIPTION:
!  User-supplied Component
!
!
!\begin{verbatim}

    module  LIS_model_componentMod

      ! ESMF Framework module
      use ESMF
      use LIS_coreMod,    only : LIS_rc, LIS_vm, LIS_masterproc, LIS_domain
      use LIS_logMod,     only : LIS_verify
      use LIS_FORC_AttributesMod
      use nldas2_forcingMod
      use LIS_field_bundleMod
      use LIS_create_gridMod, only : create_regular_grid
      use LIS_metforcing_nldas2Mod, only : number_nldas2_fields, list_nldas2_fields

      implicit none
    
      public model_register
      public model_bundle_bname
      public setMod_forcast_member_index
      public setMod_nest_index
      public setMod_findex

      real,  save :: model_gridDesc(10)
      character(len=20), parameter :: Iam = 'Model Data Comp:'
      character(len=18), parameter :: model_bundle_bname = "Model_Data_Bundle_"

      !type(ESMF_STAGGERLOC)        :: staggerloc = ESMF_STAGGERLOC_CENTER
      !type(ESMF_STAGGERLOC)        :: staggerloc = ESMF_STAGGERLOC_EDGE1
      !type(ESMF_STAGGERLOC)        :: staggerloc = ESMF_STAGGERLOC_EDGE2
       type(ESMF_STAGGERLOC), save :: staggerloc

      type(ESMF_TypeKind_Flag) :: type_kind = ESMF_TYPEKIND_R4
      integer, parameter       ::      prec = ESMF_KIND_R4
      type(ESMF_recordBundle), pointer :: recordBundle(:)

      integer :: forcast_member_index
      integer :: nest_index
      integer :: findex

!--------------------------------------------------------------------------------
    contains
!--------------------------------------------------------------------------------
!   !  The Register routine sets the subroutines to be called
!   !   as the init, run, and finalize routines.  Note that these are
!   !   private to the module.
 
    subroutine model_register(comp, rc)
        type(ESMF_GridComp) :: comp
        integer, intent(out) :: rc

        if (LIS_masterproc) print *, TRIM(Iam) // " Register routine"

        ! Register the callback routines.

        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, model_init, rc=rc)
        call LIS_verify(rc, 'ESMF_GridCompSetEntryPoint failed for init')

        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, model_run, rc=rc)
        call LIS_verify(rc, 'ESMF_GridCompSetEntryPoint failed for run')

        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, model_final, rc=rc)
        call LIS_verify(rc, 'ESMF_GridCompSetEntryPoint failed for finalize')

        if (LIS_masterproc) &
           print *, TRIM(Iam) // " Registered Initialize, Run, and Finalize routines"

    end subroutine model_register

!--------------------------------------------------------------------------------
!   !  User Comp Component created by higher level calls, here is the
!   !   Initialization routine.
 
    subroutine model_init(comp, importState, exportState, clock, rc)
      type(ESMF_GridComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

!   ! Local variables
        type(ESMF_Field)                    :: dummy
        type(ESMF_FieldBundle), allocatable :: model_bundle(:)
        type(ESMF_Grid),        allocatable :: esmfGrid(:)
        type(ESMF_ArraySpec)                :: arrayspec
        real(prec),         pointer :: ptr2D(:,:)
        real(ESMF_KIND_R8),         pointer :: coordX(:,:), coordY(:,:)
        integer                             :: i, j, counts(2), tlb(2), tub(2)
        real(ESMF_KIND_R8)                  :: min(2), max(2), dx, dy
        integer :: fieldcount
        type(ESMF_Config)  :: config
        integer :: nx, ny, N, iv
        character(len=50) :: staggerloc_var
        character(len=100) :: var_name
        character(len=2)                    :: num_st
        real                                :: dummy_array(1,1) = 0.0

      ! Initially import state contains a field with a grid but no data.

      if (LIS_masterproc) print *, TRIM(Iam) // " Model Comp Init starting"

       model_gridDesc(:) = 0.0

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

        if (LIS_masterproc) then
           print*, "<>------ Model data Grid -------<>"
           print*, "Number of nest:     ", LIS_rc%nnest
           print*, "Nx ", Nx
           print*, "Ny ", Ny
           print*, "num points along x: ", LIS_rc%gnc(1) ! LIS_rc%gridDesc(1,2)
           print*, "num points along y: ", LIS_rc%gnr(1) ! LIS_rc%gridDesc(1,3)
           print*, "lower lat:          ", LIS_rc%gridDesc(1,34), LIS_rc%minLat(n)
           print*, "lower lon:          ", LIS_rc%gridDesc(1,35), LIS_rc%minLon(n)
           print*, "upper lat:          ", LIS_rc%gridDesc(1,37)
           print*, "upper lon:          ", LIS_rc%gridDesc(1,38)
           print*, "<>------ Model data Grid -------<>"
        endif

        allocate(recordBundle(LIS_rc%nnest))
        allocate(model_bundle(LIS_rc%nnest))
        allocate(esmfGrid(LIS_rc%nnest))

        do n = 1, LIS_rc%nnest

           model_gridDesc(2) = LIS_rc%gnc(n) ! LIS_rc%gridDesc(n,2) ! Global num points along x
           model_gridDesc(3) = LIS_rc%gnr(n) ! LIS_rc%gridDesc(n,3) ! Global num points along y
           model_gridDesc(4) = LIS_rc%gridDesc(n,34) ! lower lat
           model_gridDesc(5) = LIS_rc%gridDesc(n,35) ! lower lon
           model_gridDesc(7) = LIS_rc%gridDesc(n,37)! 7) ! upper lat
           model_gridDesc(8) = LIS_rc%gridDesc(n,38) !8) ! upper lon
           model_gridDesc(9) = LIS_rc%gridDesc(n,39) ! x-grid size
           model_gridDesc(10)= LIS_rc%gridDesc(n,40) ! y-grid size

           ! Create the bundle
           write(num_st, '(i2.2)') n
           model_bundle(n) = ESMF_FieldBundleCreate(name = model_bundle_bname//num_st, rc=rc)
           recordBundle(n)%prevBundle = ESMF_FieldBundleCreate(name = 'prev_rec'//num_st, rc=rc)
           recordBundle(n)%nextBundle = ESMF_FieldBundleCreate(name = 'next_rec'//num_st, rc=rc)

           esmfGrid(n) = create_regular_grid(LIS_vm, model_gridDesc, "Model Grid", &
                             Nx, Ny, staggerloc = staggerloc)

           ! Add fields to bundles
           DO iv = 1, number_nldas2_fields
              call addTracerToBundle(model_bundle(n), esmfGrid(n), &
                                     TRIM(list_nldas2_fields(iv)), type_kind, dummy_array)
              call addTracerToBundle(recordBundle(n)%prevBundle, esmfGrid(n), &
                                     TRIM(list_nldas2_fields(iv)), type_kind, dummy_array)
              call addTracerToBundle(recordBundle(n)%nextBundle, esmfGrid(n), &
                                     TRIM(list_nldas2_fields(iv)), type_kind, dummy_array)
          ENDDO

           ! Add bundle to the state
           call ESMF_StateAdd(importState, (/model_bundle(n)/), rc=rc)
           call LIS_verify(rc, 'ESMF_StateAdd failed for bundle')
        enddo
  
      if (LIS_masterproc) print *, TRIM(Iam) // " Init returning"
   
    end subroutine model_init
!EOC
!--------------------------------------------------------------------------------
!   !  The Run routine where data is computed.
!   !
 
    subroutine model_run(comp, importState, exportState, clock, rc)
      type(ESMF_GridComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

!   ! Local variables
        type(ESMF_Field)                    :: returnedfield
        type(ESMF_FieldBundle), allocatable :: model_bundle(:)
        type(ESMF_grid),        allocatable :: esmfGrid(:)
        real(prec),         pointer :: ptr2D(:,:)
        integer                             :: n, fieldcount
        integer                             :: i, kk
        character(len=2)                    :: num_st
        character(len=100) :: var_name
        integer, save :: rec_number = 0

        if (LIS_masterproc) print *, TRIM(Iam) // " Run starting"

        allocate(model_bundle(LIS_rc%nnest))
        allocate(esmfGrid(LIS_rc%nnest))

        kk = forcast_member_index
        n  = nest_index

        !do n = 1, LIS_rc%nnest

           ! Get the FieldBundle data from the State
           write(num_st, '(i2.2)') n
           call ESMF_StateGet(importState, model_bundle_bname//num_st, model_bundle(n), rc=rc)
           call LIS_verify(rc, 'ESMF_StateGet failed to retreive model bundle')

          IF (rec_number == 0) THEN
             ! Get the fisrt record
             CALL copy_esmfBundle(recordBundle(n)%prevBundle, model_bundle(n))
             
             !var_name = TRIM(LIS_FORC_Tair%varname(1))
             !call getTracerFromBundle(recordBundle(n)%prevBundle, ptr2D, TRIM(var_name))
             !if (LIS_masterproc) print*,"Prev-"//TRIM(var_name)//": ",minval(ptr2D),maxval(ptr2D)

             call fromBundle_to_metdata(recordBundle(n)%prevBundle, n, kk, 1)
          ELSE IF (rec_number == 1) THEN
             ! Get the second record
             CALL copy_esmfBundle(recordBundle(n)%nextBundle, model_bundle(n))

             call fromBundle_to_metdata(recordBundle(n)%nextBundle, n, kk, 2)
          ELSE
             ! Because of a new record, update the previous and next record
             CALL copy_esmfBundle(recordBundle(n)%prevBundle, recordBundle(n)%nextBundle)
             CALL copy_esmfBundle(recordBundle(n)%nextBundle, model_bundle(n))

             call fromBundle_to_metdata(recordBundle(n)%prevBundle, n, kk, 1)
             call fromBundle_to_metdata(recordBundle(n)%nextBundle, n, kk, 2)

             !var_name = TRIM(LIS_FORC_Tair%varname(1))
             !call getTracerFromBundle(recordBundle(n)%prevBundle, ptr2D, TRIM(var_name))
             !if (LIS_masterproc) print*,"Prev-"//TRIM(var_name)//": ",minval(ptr2D),maxval(ptr2D)

             !call getTracerFromBundle(recordBundle(n)%nextBundle, ptr2D, TRIM(var_name))
             !if (LIS_masterproc) print*,"Next-"//TRIM(var_name)//": ",minval(ptr2D),maxval(ptr2D)

             !call getTracerFromBundle(model_bundle(n), ptr2D, TRIM(var_name))
             !if (LIS_masterproc) print*,"Curr-"//TRIM(var_name)//": ",minval(ptr2D),maxval(ptr2D)
          ENDIF
          rec_number = rec_number + 1

          !if (LIS_masterproc) then
          !   call ESMF_FieldBundleGet(model_bundle(n), fieldCount=fieldcount, rc=rc)
          !   call LIS_verify(rc, 'ESMF_FieldBundleGet failed for model bundle')
          !   print *, "***** "//TRIM(Iam) // "Number of fields in Model Bundle =", fieldcount
          !   do i = 1, fieldcount
          !      call ESMF_FieldBundleGet(model_bundle(n), i, returnedfield, rc=rc)
          !      call LIS_verify(rc, 'ESMF_FieldBundleGet failed to return field')
          !      call ESMF_FieldGet(returnedfield, name=var_name, rc=rc)
          !      call LIS_verify(rc, 'ESMF_FieldGet failed to return field name')
          !      print*, '**********> ', i, TRIM(var_name)
          !   enddo
          !endif

      if (LIS_masterproc) print *, TRIM(Iam) // " Run ending"

    end subroutine model_run


!--------------------------------------------------------------------------------
!   !  The Finalization routine where things are deleted and cleaned up.
!   !
 
    subroutine model_final(comp, importState, exportState, clock, rc)
      type(ESMF_GridComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      ! Local variables
        type(ESMF_FieldBundle), allocatable :: model_bundle(:)
        type(ESMF_grid),        allocatable :: esmfGrid(:)
        integer                             :: n
        character(len=2)                    :: num_st

      type(ESMF_Field) :: field
      real(prec), dimension(:,:), pointer :: exact_data, interp_data
      real(prec), dimension(:,:), pointer :: coordX, coordY
      integer :: haloWidth, haloUWidth(2,1), tlb(2), tub(2)

      rc = ESMF_SUCCESS

      if (LIS_masterproc) print *, TRIM(Iam) // "  Final starting"

        ! garbage collection   
        allocate(model_bundle(LIS_rc%nnest))

        do n = 1, LIS_rc%nnest

           ! Get the FieldBundle data from the State
           write(num_st, '(i2.2)') n
           call ESMF_StateGet(importState, model_bundle_bname//num_st, model_bundle(n), rc=rc)

           call ESMF_FieldBundleDestroy(model_bundle(n), rc=rc)

        enddo

      if (LIS_masterproc) print *, TRIM(Iam) // "  Final returning"
   
    end subroutine model_final

!--------------------------------------------------------------------------------
!BOP

    subroutine fromBundle_to_metdata(myBundle, n, kk, order)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
      integer, intent(in) :: kk
      integer, intent(in) :: order
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_FieldBundle),  intent(inOut) :: myBundle
!
! !DESCRIPTION:
!  This subroutine extracts data from an ESMF bundle of 2D fields
!  to populate the model array metdata (1D).
!
! !LOCAL VARIABLES:
      integer :: fieldcount, iv, r, c, rc
      integer :: i_min, j_min
      real(kind=4), pointer :: varfield(:,:)
!EOP
!--------------------------------------------------------------------------------
!BOC
      ! Get the number of fields in the bundle
      call ESMF_FieldBundleGet(myBundle, fieldCount=fieldcount, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for model bundle')

      ! Loop over all the fields in the bundle to populate the metdata array
      ! It is assumed the ordering of the fields in the bundle is the same
      ! ordering in the metdata array.
      do iv = 1, fieldcount
         call getTracerFromBundle(myBundle, varfield, iv)

         i_min = lbound(varfield, 1) ! lower bound of the first  dimension
         j_min = lbound(varfield, 2) ! lower bound of the second dimension

         !PRINT'(3i5,2f10.4)',iv,i_min,j_min,maxval(varfield),minval(varfield)
         !PRINT'(4i5)',n,kk,LIS_rc%lnr(n),LIS_rc%lnc(n)

         do r=1,LIS_rc%lnr(n)
            do c=1,LIS_rc%lnc(n)
               if (LIS_domain(n)%gindex(c,r).ne.-1) then
                  if (order.eq.1) then
                     nldas2_struc(n)%metdata1(kk,iv,LIS_domain(n)%gindex(c,r)) &
                          = varfield(i_min+c-1,j_min+r-1)
                  elseif (order.eq.2) then
                     nldas2_struc(n)%metdata2(kk,iv,LIS_domain(n)%gindex(c,r))&
                          = varfield(i_min+c-1,j_min+r-1)
                  endif
               endif
            end do
         enddo

         PRINT'(2i5,4f13.4)',kk,iv,maxval(nldas2_struc(n)%metdata1(kk,iv,:)),minval(nldas2_struc(n)%metdata1(kk,iv,:)), maxval(nldas2_struc(n)%metdata2(kk,iv,:)),minval(nldas2_struc(n)%metdata2(kk,iv,:))
         !PRINT'(3i5,2f10.4)',iv,i_min,j_min,maxval(nldas2_struc(n)%metdata1(kk,iv,:)),minval(nldas2_struc(n)%metdata1(kk,iv,:))
      enddo

    end subroutine fromBundle_to_metdata
!EOC
!--------------------------------------------------------------------------------

      subroutine setMod_forcast_member_index(kk)
          integer, intent(in) :: kk
          forcast_member_index = kk
      end subroutine setMod_forcast_member_index
     
      subroutine setMod_nest_index(n)
          integer, intent(in) :: n
          nest_index = n
      end subroutine setMod_nest_index

      subroutine setMod_findex(n)
          integer, intent(in) :: n
          findex = n
      end subroutine setMod_findex
     
    end module  LIS_model_componentMod
    
!\end{verbatim}
    
