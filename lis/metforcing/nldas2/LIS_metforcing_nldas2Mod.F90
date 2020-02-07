!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!BOP
!
! !DESCRIPTION:
!  Module for getting access to nldas2 files
!
!
!\begin{verbatim}

    module LIS_metforcing_nldas2Mod

! !USES:
       use LIS_metforcingMod, only : LIS_forc
       use grib_api
      ! ESMF Framework module
      use ESMF
      use LIS_coreMod,            only : LIS_rc, LIS_vm, LIS_masterproc
      use LIS_logMod,             only : LIS_logunit, LIS_verify, LIS_warning
      use LIS_FORC_AttributesMod
      use nldas2_forcingMod,      only : nldas2_struc
      use LIS_field_bundleMod

      implicit none
    
      private
      public get_nldas2_grib_file
      public read_nldas2b_grib_file
      public read_nldas2a_grib_file
      public   :: number_nldas2_fields
      public   :: list_nldas2_fields
      public   :: set_list_nldas2_fields

      character(len=20), parameter :: Iam = 'LIS_metforcing_nldas2Mod:'

      integer, parameter :: number_nldas2_fields = 9 ! 11
      character(len=100) :: list_nldas2_fields(number_nldas2_fields)
!
! 01 -> LIS_FORC_Tair%varname(1)
! 02 -> LIS_FORC_Qair%varname(1)
! 03 -> LIS_FORC_SWdown%varname(1)
! 04 -> LIS_FORC_LWdown%varname(1)
! 05 -> LIS_FORC_Wind_E%varname(1)
! 06 -> LIS_FORC_Wind_N%varname(1)
! 07 -> LIS_FORC_Psurf%varname(1)
! 08 -> LIS_FORC_Rainf%varname(1)
! 09 -> LIS_FORC_CRainf%varname(1)
! 10 -> LIS_FORC_PET%varname(1)
! 11 -> LIS_FORC_CAPE%varname(1)
! 12 -> LIS_FORC_Forc_Hgt%varname(1)
! 13 -> LIS_FORC_Ch%varname(1)

!EOP        
!-------------------------------------------------------------------------
    contains
!-------------------------------------------------------------------------
!BOP
      SUBROUTINE set_list_nldas2_fields()
!
! !DESCRIPTION:
! Set the list of fields to be regridded.
!EOP
!-------------------------------------------------------------------------
!BOC
       list_nldas2_fields( 1) = TRIM(LIS_FORC_Tair%varname(1))
       list_nldas2_fields( 2) = TRIM(LIS_FORC_Qair%varname(1))
       list_nldas2_fields( 3) = TRIM(LIS_FORC_LWdown%varname(1))
       list_nldas2_fields( 4) = TRIM(LIS_FORC_SWdown%varname(1))
       list_nldas2_fields( 5) = TRIM(LIS_FORC_Wind_E%varname(1))
       list_nldas2_fields( 6) = TRIM(LIS_FORC_Wind_N%varname(1))
       list_nldas2_fields( 7) = TRIM(LIS_FORC_Psurf%varname(1))
       list_nldas2_fields( 8) = TRIM(LIS_FORC_Rainf%varname(1))
       list_nldas2_fields( 9) = TRIM(LIS_FORC_CRainf%varname(1))
!       list_nldas2_fields(10) = TRIM(LIS_FORC_PET%varname(1))
!       list_nldas2_fields(11) = TRIM(LIS_FORC_CAPE%varname(1))

      END SUBROUTINE set_list_nldas2_fields
!EOC
!-------------------------------------------------------------------------
!BOP
      SUBROUTINE get_nldas2_grib_file(file_name, ftype, n, &
                                     curYear, curMonth, curDay, &
                                     curDayOfYear, curHour)
!
! !INPUT PARAMETERS:
      CHARACTER(len=1), intent(in) :: ftype ! 'A' or 'B'
      INTEGER, intent(in) :: n
      INTEGER, intent(in) :: curYear      ! current year
      INTEGER, intent(in) :: curMonth     ! current month
      INTEGER, intent(in) :: curDay       ! current day
      INTEGER, intent(in) :: curDayOfYear ! current day of the year
      INTEGER, intent(in) :: curHour      ! current hour
!
! !INPUT/OUTOUT PARAMETERS:
      CHARACTER(len=*), intent(out) :: file_name
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      character*4  :: fyr
      character*3  :: fdoy
      character*2  :: fmo, fda, fhr
      integer      :: doy2
!
!EOP
!-------------------------------------------------------------------------
!BOC
      write(unit=fyr, fmt='(i4.4)')  curYear
      write(unit=fdoy,fmt='(i3.3)')  curDayOfYear
      write(unit=fmo, fmt='(i2.2)')  curMonth
      write(unit=fda, fmt='(i2.2)')  curDay
      write(unit=fhr, fmt='(i2.2)')  curHour

      file_name = trim(nldas2_struc(n)%nldas2dir)//"/"//fyr//"/"//fdoy//&
                     "/NLDAS_FOR"//ftype//"0125_H.A"//fyr//fmo//fda//"."//fhr// &
                     "00.002.grb"

      END SUBROUTINE get_nldas2_grib_file
!EOC
!-------------------------------------------------------------------------
!BOP


    SUBROUTINE read_nldas2a_grib_file(file_name, field_bundle, n, tlb, tub)

!
! !INPUT PARAMETERS:
       integer,          intent(in) :: n           ! index of the nest
       integer,          intent(in) :: tlb(2)      ! lower bound grid points
       integer,          intent(in) :: tub(2)      ! upper bound grid points
       character(len=*), intent(in) :: file_name   ! name of the hourly NLDAS2 forecast file
!
! !INPUT/OUTPUT PARAMETERS:
       type(ESMF_FieldBundle), intent(inOut) :: field_bundle
!
! DESCRIPTION:
!   Read data from NLDAS met forcing grib file and populate and ESMF bundle.

       integer, parameter :: number_fields = 11
!
! ! LOCAL VARIABLES:
       integer                     :: iv, ftn
       integer                     :: nldas2
       integer                     :: k,t,c,r,iret,rc
       integer                     :: var_index
       real                        :: missingValue
       integer                     :: num_infile_vars
       integer                     :: igrib
       logical                     :: pcp_flag, var_found
       logical                     :: var_status(number_fields)
       logical                     :: file_exists
       integer                     :: pds5(number_fields)
       integer                     :: pds7(number_fields)
       integer                     :: pds2(number_fields)
       integer                     :: pds5_val, pds7_val, pds2_val
       real, allocatable           :: var1D(:)
       character(len=100)          :: var_name
       real(ESMF_KIND_R4), pointer :: ptr2Dglob(:,:)
       real(ESMF_KIND_R4), pointer :: ptr2Dloc(:,:)

!EOP
!-------------------------------------------------------------------------
       ! Total number of grid points
       nldas2 = nldas2_struc(n)%ncold * nldas2_struc(n)%nrold

       ! Order of variables to be assigned in the NLDAS-2 Forcing "B" files.
       ! Note that this is NOT the order of the fields in the actual files,
       ! but the order in which they are assigned to "metdata". - dmm
       pds5 = (/ 011,051,204,205,033,034,001,061,153,228,157/) !parameter
       pds7 = (/ 002,002,000,000,010,010,000,000,000,000,180/) !htlev2
       pds2 = (/  84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84/)


       var_status = .false.

       INQUIRE (file=file_name, exist=file_exists)
       if (.NOT. file_exists) then
          write(LIS_logunit,*) '[ERR] File does not exist: ',trim(file_name)
          return
       else 
          call grib_open_file(ftn,trim(file_name),'r',iret)
          if(iret.ne.0) then
             write(LIS_logunit,*) &
                  '[ERR] Could not open file: ',trim(file_name)
             return
          endif

          IF ( LIS_masterproc) print*,"---<> Opened file: ", TRIM(file_name)

          call grib_count_in_file(ftn,num_infile_vars,iret)
          call LIS_verify(iret, 'error in grib_count_in_file in read_nldas2_grib_file')

          allocate(var1D(nldas2_struc(n)%ncold*nldas2_struc(n)%nrold), stat=iret)
          call LIS_verify(iret, 'Cannot allocate var1D')

          IF ( LIS_masterproc) print*,"---<> Loop over the variables in the file "

          do k=1,num_infile_vars
             call grib_new_from_file(ftn, igrib, iret)
             call LIS_warning(iret, 'error in grib_new_from_file in read_nldas2_grib_file')
             if(iret.ne.0) then
                write(LIS_logunit,*) &
                     '[ERR] Could not retrieve entries in file: ',trim(file_name)
                deallocate(var1D)
                return
             endif

             call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
             call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_nldas2_grib_file')

             call grib_get(igrib,'level',pds7_val,rc)
             call LIS_verify(rc, 'error in grib_get: level in read_nldas2_grib_file')

             call grib_get(igrib,'generatingProcessIdentifier',pds2_val,rc)
             call LIS_verify(rc, 'error in grib_get: generatingProcessIdentifier in read_nldas2_grib_file')

             var_found = .false.
             do iv=1,number_fields
                if((pds5_val.eq.pds5(iv)).and.&
                     (pds7_val.eq.pds7(iv)).and.&
                     (pds2_val.eq.pds2(iv))) then
                   var_found = .true.
                   var_index = iv
                   var_status(iv) = .true.
                   exit
                endif
             enddo

             var1D = LIS_rc%udef
             call grib_get(igrib, 'values', var1D, rc)
             call LIS_warning(rc, 'error in grib_get:values in read_nldas2_grib_file')

             if(rc.ne.0) then
                write(LIS_logunit,*) &
                     '[ERR] Could not retrieve entries in file: ',trim(file_name)
                deallocate(var1D)
                return
             endif

             call grib_get(igrib,'missingValue',missingValue,rc)
             call LIS_verify(rc, 'error in grib_get:missingValue in read_nldas2_grib_file')

             call grib_release(igrib,rc)
             call LIS_verify(rc, 'error in grib_release in read_nldas2_grib_file')

             IF ((var_index >= 1) .AND. (var_index <= number_nldas2_fields)) THEN
                var_name = list_nldas2_fields(var_index)
             ELSE
                var_name = ""
             ENDIF
   
             IF (var_name .NE. "") THEN

                !write(unit=LIS_logunit,fmt=*)'[INFO] dealing with .. ',trim(var_name)

                allocate(ptr2Dglob(nldas2_struc(n)%ncold, nldas2_struc(n)%nrold), stat=iret)
                call LIS_verify(iret, 'Cannot allocate ptr2Dglob for '//TRIM(var_name))
                allocate(ptr2Dloc(tlb(1):tub(1), tlb(2):tub(2)), stat=iret)
                call LIS_verify(iret, 'Cannot allocate ptr2Dloc for '//TRIM(var_name))
                ptr2Dglob = reshape(var1D, (/ nldas2_struc(n)%ncold, nldas2_struc(n)%nrold /) )
                ptr2Dloc(:,:) = ptr2Dglob(tlb(1):tub(1), tlb(2):tub(2))
                IF ( LIS_masterproc)  print*,"FORCING-"//TRIM(var_name)//": ",minval(ptr2Dloc),maxval(ptr2Dloc)
                call updateTracerToBundle(field_bundle, ptr2Dloc, TRIM(var_name))
                DEALLOCATE(ptr2Dglob, ptr2Dloc)
             ENDIF
          enddo

          DEALLOCATE(var1D)
       endif

    END SUBROUTINE read_nldas2a_grib_file
!EOC
!-------------------------------------------------------------------------

!BOP


    SUBROUTINE read_nldas2b_grib_file(file_name, field_bundle, n, tlb, tub)

!
! !INPUT PARAMETERS:
       integer,          intent(in) :: n           ! index of the nest
       integer,          intent(in) :: tlb(2)      ! lower bound grid points
       integer,          intent(in) :: tub(2)      ! upper bound grid points
       character(len=*), intent(in) :: file_name   ! name of the hourly NLDAS2 forecast file
!
! !INPUT/OUTPUT PARAMETERS:
       type(ESMF_FieldBundle), intent(inOut) :: field_bundle
!
! DESCRIPTION:
!   Read data from NLDAS met forcing grib file and populate and ESMF bundle.

       integer, parameter :: number_fields = 10
!
! ! LOCAL VARIABLES:
       integer                     :: iv, ftn
       integer                     :: nldas2
       integer                     :: k,t,c,r,iret,rc
       integer                     :: var_index
       real                        :: missingValue
       integer                     :: num_infile_vars
       integer                     :: igrib
       logical                     :: pcp_flag, var_found
       logical                     :: var_status(10)
       logical                     :: file_exists
       integer                     :: pds5(number_fields)
       integer                     :: pds7(number_fields)
       integer                     :: pds2(number_fields)
       integer                     :: pds5_val, pds7_val, pds2_val
       real, allocatable           :: var1D(:)
       character(len=100)          :: var_name
       real(ESMF_KIND_R4), pointer :: ptr2Dglob(:,:)
       real(ESMF_KIND_R4), pointer :: ptr2Dloc(:,:)

!EOP
!-------------------------------------------------------------------------
       ! Total number of grid points
       nldas2 = nldas2_struc(n)%ncold * nldas2_struc(n)%nrold

       ! Order of variables to be assigned in the NLDAS-2 Forcing "B" files.
       ! Note that this is NOT the order of the fields in the actual files,
       ! but the order in which they are assigned to "metdata". - dmm
       pds5 = (/ 011,051,204,179,033,034,001,061,063,007/) !parameter
       pds7 = (/ 001,001,000,000,001,001,001,000,000,001/) !htlev2
       pds2 = (/  84, 84, 84, 84, 84, 84, 84, 84, 84, 84/)

       var_status = .false.

       INQUIRE (file=file_name, exist=file_exists)
       if (.NOT. file_exists) then
          write(LIS_logunit,*) '[ERR] File does not exist: ',trim(file_name)
          return
       else 
          call grib_open_file(ftn,trim(file_name),'r',iret)
          if(iret.ne.0) then
             write(LIS_logunit,*) &
                  '[ERR] Could not open file: ',trim(file_name)
             return
          endif

          !IF ( LIS_masterproc) print*,"---<> Opened file: ", TRIM(file_name)

          call grib_count_in_file(ftn,num_infile_vars,iret)
          call LIS_verify(iret, 'error in grib_count_in_file in read_nldas2_grib_file')

          allocate(var1D(nldas2_struc(n)%ncold*nldas2_struc(n)%nrold))

          do k=1,num_infile_vars
             call grib_new_from_file(ftn, igrib, iret)
             call LIS_warning(iret, 'error in grib_new_from_file in read_nldas2_grib_file')
             if(iret.ne.0) then
                write(LIS_logunit,*) &
                     '[ERR] Could not retrieve entries in file: ',trim(file_name)
                deallocate(var1D)
                return
             endif

             call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
             call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_nldas2_grib_file')

             call grib_get(igrib,'level',pds7_val,rc)
             call LIS_verify(rc, 'error in grib_get: level in read_nldas2_grib_file')

             call grib_get(igrib,'generatingProcessIdentifier',pds2_val,rc)
             call LIS_verify(rc, 'error in grib_get: generatingProcessIdentifier in read_nldas2_grib_file')

             var_found = .false.
             do iv=1,number_fields
                if((pds5_val.eq.pds5(iv)).and.&
                     (pds7_val.eq.pds7(iv)).and.&
                     (pds2_val.eq.pds2(iv))) then
                   var_found = .true.
                   var_index = iv
                   var_status(iv) = .true.
                   exit
                endif
             enddo

             var1D = LIS_rc%udef
             call grib_get(igrib, 'values', var1D, rc)
             call LIS_warning(rc, 'error in grib_get:values in read_nldas2_grib_file')

             if(rc.ne.0) then
                write(LIS_logunit,*) &
                     '[ERR] Could not retrieve entries in file: ',trim(file_name)
                deallocate(var1D)
                return
             endif

             call grib_get(igrib,'missingValue',missingValue,rc)
             call LIS_verify(rc, 'error in grib_get:missingValue in read_nldas2_grib_file')

             call grib_release(igrib,rc)
             call LIS_verify(rc, 'error in grib_release in read_nldas2_grib_file')

             var_name = ""
             SELECT CASE (var_index)
                CASE(1)
                   var_name = TRIM(LIS_FORC_Tair%varname(1))
                CASE(2)
                   var_name = TRIM(LIS_FORC_Qair%varname(1))
                CASE(3)
                   var_name = TRIM(LIS_FORC_LWdown%varname(1))
                CASE(4)
                   var_name = TRIM(LIS_FORC_SWdown%varname(1))
                CASE(5)
                   var_name = TRIM(LIS_FORC_Wind_E%varname(1))
                CASE(6)
                   var_name = TRIM(LIS_FORC_Wind_N%varname(1))
                CASE(7)
                   var_name = TRIM(LIS_FORC_Psurf%varname(1))
                CASE(8)
                   var_name = TRIM(LIS_FORC_Rainf%varname(1))
                CASE(9)
                   var_name = TRIM(LIS_FORC_CRainf%varname(1))
             END SELECT
   
             IF (var_name .NE. "") THEN
                allocate(ptr2Dglob(nldas2_struc(n)%ncold, nldas2_struc(n)%nrold))
                allocate(ptr2Dloc(tlb(1):tub(1), tlb(2):tub(2)))
                ptr2Dglob = reshape(var1D, (/ nldas2_struc(n)%ncold, nldas2_struc(n)%nrold /) )
                ptr2Dloc(:,:) = ptr2Dglob(tlb(1):tub(1), tlb(2):tub(2))
                !IF ( LIS_masterproc)  print*,"FORCING-"//TRIM(var_name)//": ",minval(ptr2Dloc),maxval(ptr2Dloc)
                call updateTracerToBundle(field_bundle, ptr2Dloc, TRIM(var_name))
             ENDIF
          enddo

          DEALLOCATE(var1D)
       endif

    END SUBROUTINE read_nldas2b_grib_file
!EOC
!-------------------------------------------------------------------------

    end module LIS_metforcing_nldas2Mod
    
!\end{verbatim}
    
