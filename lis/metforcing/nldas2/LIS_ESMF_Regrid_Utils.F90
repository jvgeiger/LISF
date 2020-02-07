!-------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION:
!
      module LIS_ESMF_Regrid_Utils

      use ESMF
      use LIS_coreMod,    only : LIS_rc, LIS_vm, LIS_masterproc
      use LIS_logMod,     only : LIS_verify
    
      implicit none
    
      private
      public runESMF_Regridding
      public createESMF_RouteHandle
        
      interface createESMF_RouteHandle
         module procedure createESMF_RouteHandleField
         module procedure createESMF_RouteHandleBundle
      end interface

      interface runESMF_Regridding
         module procedure runESMF_RegriddingField
         module procedure runESMF_RegriddingBundle
      end interface

      REAL, PARAMETER :: UNDEF_VALUE = 9999.000 !1.0e15
      character(len=ESMF_MAXSTR), parameter :: Iam = 'ESMF Regridding Utility:'
!EOP
!-------------------------------------------------------------------------
    contains
!-------------------------------------------------------------------------
!BOP    
    subroutine createESMF_RouteHandleBundle(srcBundle, dstBundle, regridMethod, routehandle, rc)

      type(ESMF_RegridMethod_Flag), intent(in)    :: regridMethod
      type(ESMF_FieldBundle),       intent(inOut) :: srcBundle
      type(ESMF_FieldBundle),       intent(inOut) :: dstBundle
      type(ESMF_RouteHandle),       intent(out)   :: routehandle
      integer,                      intent(out)   :: rc
!
! !DESCRIPTION:
! Given the source and destination ESMF bundle, create the routehandle
!
! !LOACAL VARIABLES:
      integer :: itemcount, n
      integer :: src_fieldcount, dst_fieldcount
      integer                    :: index_tracer
      type(ESMF_FIELD)           :: src_field
      type(ESMF_FIELD)           :: dst_field
      character(len=ESMF_MAXSTR) :: src_fieldname, dst_fieldname
!EOP
!-------------------------------------------------------------------------
!BOC

      rc = ESMF_SUCCESS

      ! Verify that the two bundles have the same fields
      !-------------------------------------------------

      call ESMF_FieldBundleGet(srcBundle, fieldCount=src_fieldcount, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for source bundle')

      call ESMF_FieldBundleGet(dstBundle, fieldCount=dst_fieldcount, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleGet failed for destination bundle')

      IF (src_fieldcount .NE. dst_fieldcount) THEN
         rc = -1
         call LIS_verify(rc, 'The two bundles do not have the same number of fields')
      END IF
       
      DO index_tracer = 1, src_fieldcount
         call ESMF_FieldBundleGet(srcBundle, &
                                  fieldIndex = index_tracer, &
                                  field      = src_field, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldBundleGet failed for source bundle')

         call ESMF_FieldGet(src_field, name = src_fieldname, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldGet failed to get the source field name')

         call ESMF_FieldBundleGet(dstBundle, &
                                  fieldIndex = index_tracer, &
                                  field      = dst_field, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldBundleGet failed for destination bundle')

         call ESMF_FieldGet(dst_field, name = dst_fieldname, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldGet failed to get the destination field name')

         IF (TRIM(src_fieldname) .NE. TRIM(dst_fieldname)) THEN
            rc = -1
            call LIS_verify(rc, 'Source field '//TRIM(src_fieldname)//' is not the same as Destination field '//TRIM(dst_fieldname))
         END IF
         IF (LIS_masterproc) then
            print *, "    *",index_tracer,TRIM(src_fieldname)
         END IF

      ENDDO

      ! Determine the routehandle
      !--------------------------
      call ESMF_FieldBundleRegridStore(srcFieldBundle = srcBundle, &
                                       dstFieldBundle = dstBundle, &
                                       routeHandle    = routehandle, &
                                       regridmethod   = regridMethod, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldBundleRegridStore failed')


    end subroutine createESMF_RouteHandleBundle
!EOC
!-------------------------------------------------------------------------
!BOP    
    subroutine createESMF_RouteHandleField(srcField, dstField, regridMethod, routehandle, &
                                           dynamicMask, rc)

      type(ESMF_RegridMethod_Flag), intent(in)    :: regridMethod
      type(ESMF_Field),             intent(inOut) :: srcField
      type(ESMF_Field),             intent(inOut) :: dstField
      type(ESMF_DynamicMask),       intent(inOut) :: dynamicMask
      type(ESMF_RouteHandle),       intent(out)   :: routehandle
      integer,                      intent(out)   :: rc
!
! !DESCRIPTION:
! Given the source and destination ESMF fields, create the routehandle.

      integer :: srcTerm
!
!EOP
!-------------------------------------------------------------------------
!BOC

      rc = ESMF_SUCCESS

      srcTerm = 0
      call ESMF_FieldRegridStore(srcField          = srcfield, &
                                 dstField          = dstfield, &
                                 srcTermProcessing = srcTerm, &
                                 srcMaskValues     = [0], &
                                 unmappedAction    = ESMF_UNMAPPEDACTION_IGNORE, &
                                 routeHandle       = routehandle, &
                                 regridmethod      = regridMethod, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldRegridStore failed')

      call ESMF_DynamicMaskSetR4R8R4(dynamicMask, &
                                     simpleDynMaskProc, &
                                     dynamicSrcMaskValue=UNDEF_VALUE, rc=rc)

    end subroutine createESMF_RouteHandleField
!EOC
!-------------------------------------------------------------------------
!BOP

    subroutine runESMF_RegriddingBundle(srcBundle, dstBundle, routehandle, rc)

      type(ESMF_FieldBundle),       intent(in) :: srcBundle
      type(ESMF_FieldBundle),       intent(inOut) :: dstBundle
      type(ESMF_RouteHandle),       intent(inOut) :: routehandle
      integer,                      intent(out)   :: rc

!EOP
!-------------------------------------------------------------------------
!BOC
       call ESMF_FieldBundleRegrid(srcFieldBundle = srcBundle, &
                                   dstFieldBundle = dstBundle, &
                                   routehandle    = routehandle,  &
                                   checkflag      = .TRUE., rc=rc)
       call LIS_verify(rc, 'ESMF_FieldBundleRegrid failed')

    end subroutine runESMF_RegriddingBundle
!EOC
!-------------------------------------------------------------------------
!BOP

    subroutine runESMF_RegriddingField(srcField, dstField, routehandle, &
                                       dynamicMask, rc)

      type(ESMF_RouteHandle), intent(inOut) :: routehandle
      type(ESMF_Field),       intent(in) :: srcField
      type(ESMF_Field),       intent(inOut) :: dstField
      type(ESMF_DynamicMask), intent(inOut) :: dynamicMask
      integer,                intent(out)   :: rc

!EOP
!-------------------------------------------------------------------------
!BOC
       call ESMF_FieldRegrid(srcField      = srcField,        &
                             dstField      = dstField,        &
                             routehandle   = routehandle,     &  
                             dynamicMask   = dynamicMask,     &
                             zeroregion    = ESMF_REGION_SELECT, &
                             termorderflag = ESMF_TERMORDER_SRCSEQ, &
                             checkflag   = .TRUE., rc=rc)
       call LIS_verify(rc, 'ESMF_FieldRegrid failed')

    end subroutine runESMF_RegriddingField
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

    end module LIS_ESMF_Regrid_Utils
    
!\end{verbatim}
    
