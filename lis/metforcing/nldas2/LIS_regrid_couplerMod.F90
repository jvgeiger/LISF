!-------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION:
!  User-supplied Coupler
!
      module LIS_regrid_couplerMod

      use ESMF
      use LIS_coreMod,    only : LIS_rc, LIS_vm, LIS_masterproc
      use LIS_logMod,     only : LIS_verify
      use LIS_metforcing_componentMod, only : forcing_bundle_bname
      use LIS_model_componentMod,      only : model_bundle_bname
    
      implicit none
    
      public regridCoupler_register
        
      type esmf_regrid_params
           ! global data
           type(ESMF_RouteHandle),      allocatable :: routehandle(:)
           type(ESMF_RegridMethod_Flag)             :: regridMethod 
      end type

      type(esmf_regrid_params), save :: regrid_params
      character(len=25), parameter :: Iam = 'ESMF Regridding Coupler:'

!-------------------------------------------------------------------------
    contains
!-------------------------------------------------------------------------
!   !  The Register routine sets the subroutines to be called
!   !   as the init, run, and finalize routines.  Note that these are
!   !   private to the module.
 
    subroutine regridCoupler_register(comp, rc)
      type(ESMF_CplComp) :: comp
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS
      if (LIS_masterproc) print *, TRIM(Iam) // ": in user setservices routine"

      ! Register the callback routines.
      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, regridCoupler_init, rc=rc)
      call LIS_verify(rc, 'ESMF_GridCompSetEntryPoint failed for regridCoupler_init')

      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_RUN, regridCoupler_run, rc=rc)
      call LIS_verify(rc, 'ESMF_GridCompSetEntryPoint failed for regridCoupler_run')

      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, regridCoupler_final, rc=rc)
      call LIS_verify(rc, 'ESMF_GridCompSetEntryPoint failed for regridCoupler_finalize')

      if (LIS_masterproc) print *, TRIM(Iam) // ": Registered Initialize, Run, and Finalize routines"

    end subroutine regridCoupler_register

!-------------------------------------------------------------------------
!   !User Comp Component created by higher level calls, here is the
!   ! Initialization routine.
 
    
    subroutine regridCoupler_init(comp, importState, exportState, clock, rc)
      type(ESMF_CplComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      ! Local variables
      integer :: itemcount, fieldcount, n
      type(ESMF_FieldBundle), allocatable :: met_forcing_bundle(:)
      type(ESMF_FieldBundle), allocatable :: model_bundle(:)

      type(ESMF_Grid) :: srcGrid, dstGrid
      integer, dimension(2) :: clbnd, cubnd
      character(len=20) :: regridMethod_var
      character(len=2) :: num_st
      integer                    :: index_tracer
      type(ESMF_FIELD)           :: field_tracer
      type(ESMF_FIELD)           :: field_model
      type(ESMF_FIELD)           :: field_forcing
      character(len=100)         :: name_tracer

      rc = ESMF_SUCCESS

      if (LIS_masterproc) print *, TRIM(Iam) // ": User Coupler Init starting"

      regrid_params%regridMethod = ESMF_REGRIDMETHOD_BILINEAR
      ! regrid_params%regridMethod = ESMF_REGRIDMETHOD_NEAREST_STOD

      allocate(met_forcing_bundle(LIS_rc%nnest), STAT=rc)
      call LIS_verify(rc, 'Failed to allocate met_forcing_bundle')

      allocate(model_bundle(LIS_rc%nnest), STAT=rc)
      call LIS_verify(rc, 'Failed to allocate model_bundle')

      allocate(regrid_params%routehandle(LIS_rc%nnest), STAT=rc)
      call LIS_verify(rc, 'Failed to allocate routehandle')

      do n = 1, LIS_rc%nnest
         write(num_st, '(i2.2)') n
         call ESMF_StateGet(importState, forcing_bundle_bname//num_st, &
                            met_forcing_bundle(n), rc=rc)
         call LIS_verify(rc, 'ESMF_StateGet failed for forcing bundle')

         call ESMF_FieldBundleGet(met_forcing_bundle(n), &
                                  fieldIndex = 1, &
                                  field      = field_forcing, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldBundleGet failed for forcing field')

         if (LIS_masterproc) then
            call ESMF_FieldBundleGet(met_forcing_bundle(n), fieldCount=fieldcount, rc=rc)
            call LIS_verify(rc, 'ESMF_FieldBundleGet failed for forcing bundle')
            print *, "---> Number of fields in Forcing data Bundle =", fieldcount

            DO index_tracer = 1, fieldcount
               call ESMF_FieldBundleGet(met_forcing_bundle(n), &
                                        fieldIndex = index_tracer, &
                                        field      = field_tracer, rc=rc)
               call LIS_verify(rc, 'ESMF_FieldBundleGet failed for forcing bundle')

               call ESMF_FieldGet(field_tracer, name = name_tracer, rc=rc)
               call LIS_verify(rc, 'ESMF_FieldGet failed to get the tracer name')
               print *, "    *",index_tracer,TRIM(name_tracer)
            ENDDO
         endif

         call ESMF_StateGet(exportState, model_bundle_bname//num_st, &
                            model_bundle(n), rc=rc)
         call LIS_verify(rc, 'ESMF_StateGet failed for model bundle')

         call ESMF_FieldBundleGet(model_bundle(n), &
                                  fieldIndex = 1, &
                                  field      = field_model, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldBundleGet failed for forcing field')

         if (LIS_masterproc) then
            call ESMF_FieldBundleGet(model_bundle(n), fieldCount=fieldcount, rc=rc)
            call LIS_verify(rc, 'ESMF_FieldBundleGet failed for model bundle')
            print *, "---> Number of fields in Model data Bundle =", fieldcount
            DO index_tracer = 1, fieldcount
               call ESMF_FieldBundleGet(model_bundle(n), &
                                        fieldIndex = index_tracer, &
                                        field      = field_tracer, rc=rc)
               call LIS_verify(rc, 'ESMF_FieldBundleGet failed for forcing bundle')
            
               call ESMF_FieldGet(field_tracer, name = name_tracer, rc=rc)
               call LIS_verify(rc, 'ESMF_FieldGet failed to get the tracer name')
               print *, "    *",index_tracer,TRIM(name_tracer)
            ENDDO
         endif

         call ESMF_FieldBundleRegridStore(srcFieldBundle = met_forcing_bundle(n), &
                                          dstFieldBundle = model_bundle(n), &
         !call ESMF_FieldRegridStore(srcField     = field_forcing, &
         !                           dstField     = field_model, &
                                    routeHandle  = regrid_params%routehandle(n), &
                                    regridmethod = regrid_params%regridMethod, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldBundleRegridStore failed')
      enddo

      if (LIS_masterproc) print *, TRIM(Iam) // ": User Coupler Init returning"

      deallocate(model_bundle)
      deallocate(met_forcing_bundle)

    end subroutine regridCoupler_init


!-------------------------------------------------------------------------
!   !  The Run routine where data is coupled.
!   !
 
    subroutine regridCoupler_run(comp, importState, exportState, clock, rc)
      type(ESMF_CplComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_FieldBundle), allocatable :: met_forcing_bundle(:)
      type(ESMF_FieldBundle), allocatable :: model_bundle(:)
      integer :: n
      character(len=2) :: num_st

      if (LIS_masterproc) print *, TRIM(Iam) // ": User Coupler Run starting"

      allocate(met_forcing_bundle(LIS_rc%nnest))
      allocate(model_bundle(LIS_rc%nnest))

       do n = 1, LIS_rc%nnest
          write(num_st, '(i2.2)') n
          call ESMF_StateGet(importState, forcing_bundle_bname//num_st, &
                             met_forcing_bundle(n), rc=rc)
          call LIS_verify(rc, 'ESMF_StateGet failed for forcing bundle')

          call ESMF_StateGet(exportState, model_bundle_bname//num_st, &
                             model_bundle(n), rc=rc)
          call LIS_verify(rc, 'ESMF_StateGet failed for model bundle')

          call ESMF_FieldBundleRegrid(srcFieldBundle = met_forcing_bundle(n), &
                                      dstFieldBundle = model_bundle(n), &
                                      routehandle    = regrid_params%routehandle(n),  &
                                      checkflag      = .TRUE., rc=rc)
          call LIS_verify(rc, 'ESMF_FieldBundleRegrid failed')
      enddo

 
      if (LIS_masterproc) print *, TRIM(Iam) // ": User Coupler Run returning"

    end subroutine regridCoupler_run

!EOC
!-------------------------------------------------------------------------
!   !  The Finalization routine where things are deleted and cleaned up.
!   !
 
    subroutine regridCoupler_final(comp, importState, exportState, clock, rc)
      type(ESMF_CplComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      ! Local variables
      integer :: n

      if (LIS_masterproc) print *, TRIM(Iam) // ": User Coupler Final starting"
   
      do n = 1, LIS_rc%nnest
         ! Release resources stored for the Regridding.
         call ESMF_FieldRegridRelease(regrid_params%routehandle(n), rc=rc)
         call LIS_verify(rc, 'ESMF_FieldRegridRelease failed')
      enddo

      if (LIS_masterproc) print *, TRIM(Iam) // ": User Coupler Final returning"

    end subroutine regridCoupler_final

    end module LIS_regrid_couplerMod
    
!\end{verbatim}
    
