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
      use LIS_ESMF_Regrid_Utils, only : runESMF_Regridding, createESMF_RouteHandle
    
      implicit none
    
      public regridCoupler_register
        
      type esmf_regrid_params
           ! global data
           type(ESMF_RouteHandle),      allocatable :: routehandle(:)
           type(ESMF_DynamicMask),      allocatable :: dynamicMask(:)
           type(ESMF_RegridMethod_Flag)             :: regridMethod 
      end type

      type(esmf_regrid_params), save :: regrid_params
      character(len=25), parameter :: Iam = 'ESMF Regridding Coupler:'

      type(ESMF_DynamicMask) :: dynamicMask
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
      type(ESMF_FIELD)           :: srcField
      type(ESMF_FIELD)           :: dstField
      character(len=100)         :: name_tracer

      rc = ESMF_SUCCESS

      if (LIS_masterproc) print *, TRIM(Iam) // ": User Coupler Init starting"

      !regrid_params%regridMethod = ESMF_REGRIDMETHOD_BILINEAR
      regrid_params%regridMethod = ESMF_REGRIDMETHOD_NEAREST_STOD

      allocate(met_forcing_bundle(LIS_rc%nnest), STAT=rc)
      call LIS_verify(rc, 'Failed to allocate met_forcing_bundle')

      allocate(model_bundle(LIS_rc%nnest), STAT=rc)
      call LIS_verify(rc, 'Failed to allocate model_bundle')

      allocate(regrid_params%routehandle(LIS_rc%nnest), STAT=rc)
      call LIS_verify(rc, 'Failed to allocate routehandle')

      allocate(regrid_params%dynamicMask(LIS_rc%nnest), STAT=rc)
      call LIS_verify(rc, 'Failed to allocate dynamicMask')

      ! Option for using an ESMF Bundle
      !--------------------------------
!      do n = 1, LIS_rc%nnest
!         write(num_st, '(i2.2)') n
!         call ESMF_StateGet(importState, forcing_bundle_bname//num_st, &
!                            met_forcing_bundle(n), rc=rc)
!         call LIS_verify(rc, 'ESMF_StateGet failed for forcing bundle')
!
!
!         call ESMF_StateGet(exportState, model_bundle_bname//num_st, &
!                            model_bundle(n), rc=rc)
!         call LIS_verify(rc, 'ESMF_StateGet failed for model bundle')
!
!         call createESMF_RouteHandle(met_forcing_bundle(n), model_bundle(n), &
!                                     regrid_params%regridMethod, &
!                                     regrid_params%routehandle(n), rc)
!         call LIS_verify(rc, 'createESMF_RouteHandle failed')
!
!      enddo

      ! Option for using an ESMF Field
      !-------------------------------
      do n = 1, LIS_rc%nnest
         write(num_st, '(i2.2)') n
         ! Get first field from met forcing bundle
         call ESMF_StateGet(importState, forcing_bundle_bname//num_st, &
                            met_forcing_bundle(n), rc=rc)
         call LIS_verify(rc, 'ESMF_StateGet failed for forcing bundle')

         call ESMF_FieldBundleGet(met_forcing_bundle(n), &
                                  fieldIndex = 1, &
                                  field      = srcField, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldBundleGet failed for forcing field')

         ! Get first field from model bundle
         call ESMF_StateGet(exportState, model_bundle_bname//num_st, &
                            model_bundle(n), rc=rc)
         call LIS_verify(rc, 'ESMF_StateGet failed for model bundle')

         call ESMF_FieldBundleGet(model_bundle(n), &
                                  fieldIndex = 1, &
                                  field      = dstField, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldBundleGet failed for model field')

         ! Determine the routehandle
         call createESMF_RouteHandle(srcField, dstField,  regrid_params%regridMethod, &
                                     regrid_params%routehandle(n), dynamicMask, rc)
         call LIS_verify(rc, 'createESMF_RouteHandle failed')
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
      integer :: n, fieldcount, i
      type(ESMF_FIELD)           :: srcField
      type(ESMF_FIELD)           :: dstField
      character(len=2) :: num_st

      if (LIS_masterproc) print *, TRIM(Iam) // ": User Coupler Run starting"

      allocate(met_forcing_bundle(LIS_rc%nnest))
      allocate(model_bundle(LIS_rc%nnest))

      ! Regridding with the ESMF Bundle option
      !---------------------------------------
!      do n = 1, LIS_rc%nnest
!         write(num_st, '(i2.2)') n
!         call ESMF_StateGet(importState, forcing_bundle_bname//num_st, &
!                            met_forcing_bundle(n), rc=rc)
!         call LIS_verify(rc, 'ESMF_StateGet failed for forcing bundle')
!
!         call ESMF_StateGet(exportState, model_bundle_bname//num_st, &
!                            model_bundle(n), rc=rc)
!         call LIS_verify(rc, 'ESMF_StateGet failed for model bundle')
!
!         call runESMF_Regridding(met_forcing_bundle(n), model_bundle(n), &
!                                 regrid_params%routehandle(n), rc)
!
!         call LIS_verify(rc, 'runESMF_Regridding failed')
!      enddo

      ! Regridding with the ESMF Field option
      !---------------------------------------
      do n = 1, LIS_rc%nnest
         write(num_st, '(i2.2)') n
         call ESMF_StateGet(importState, forcing_bundle_bname//num_st, &
                            met_forcing_bundle(n), rc=rc)
         call LIS_verify(rc, 'ESMF_StateGet failed for forcing bundle')

         call ESMF_StateGet(exportState, model_bundle_bname//num_st, &
                            model_bundle(n), rc=rc)
         call LIS_verify(rc, 'ESMF_StateGet failed for model bundle')

         ! Get the number of fields to regrid
         call ESMF_FieldBundleGet(met_forcing_bundle(n), fieldCount=fieldcount, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldBundleGet failed for forcing bundle')

         DO i = 1, fieldcount
            call ESMF_FieldBundleGet(met_forcing_bundle(n), &
                                     fieldIndex = i, &
                                     field      = srcField, rc=rc)
            call LIS_verify(rc, 'ESMF_FieldBundleGet failed for forcing field')

            call ESMF_FieldBundleGet(model_bundle(n), &
                                     fieldIndex = i, &
                                     field      = dstField, rc=rc)
            call LIS_verify(rc, 'ESMF_FieldBundleGet failed for model field')

            call runESMF_Regridding(srcField, dstField, regrid_params%routehandle(n), &
                                    dynamicMask, rc)
            call LIS_verify(rc, 'runESMF_Regridding failed')
         END DO
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
!EOC
!-------------------------------------------------------------------------
!BOP
   subroutine simpleDynMaskProc(dynamicMaskList, dynamicSrcMaskValue, &
      dynamicDstMaskValue, rc)
      type(ESMF_DynamicMaskElementR4R8R4), pointer        :: dynamicMaskList(:)
      real(ESMF_KIND_R4),            intent(in), optional :: dynamicSrcMaskValue
      real(ESMF_KIND_R4),            intent(in), optional :: dynamicDstMaskValue
      integer,                       intent(out)          :: rc
      integer :: i, j
      real(ESMF_KIND_R8)  :: renorm
!EOP
!-------------------------------------------------------------------------
!BOC
      if (associated(dynamicMaskList)) then
         do i=1, size(dynamicMaskList)
            dynamicMaskList(i)%dstElement = 0.d0 ! set to zero
            renorm = 0.d0 ! reset
            do j=1, size(dynamicMaskList(i)%factor)
               !if (.not. &
                    !match(dynamicSrcMaskValue,dynamicMaskList(i)%srcElement(j))) then
                   !dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement &
                   !+ dynamicMaskList(i)%factor(j) &
                   !* dynamicMaskList(i)%srcElement(j)
                   !renorm = renorm + dynamicMaskList(i)%factor(j)
               !endif
               if (dynamicSrcMaskValue /= dynamicMaskList(i)%srcElement(j)) then
                  dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement &
                         + dynamicMaskList(i)%factor(j) &
                         * dynamicMaskList(i)%srcElement(j)
                  renorm = renorm + dynamicMaskList(i)%factor(j)
               endif
            enddo
            if (renorm > 0.d0) then
               dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement / renorm
            else if (present(dynamicSrcMaskValue)) then
               dynamicMaskList(i)%dstElement = dynamicSrcMaskValue
            else
               rc = ESMF_RC_ARG_BAD  ! error detected
               return
            endif
         enddo
      endif
      ! return successfully
      rc = ESMF_SUCCESS
    end subroutine

!EOC
!-------------------------------------------------------------------------

    end module LIS_regrid_couplerMod
    
!\end{verbatim}
    
