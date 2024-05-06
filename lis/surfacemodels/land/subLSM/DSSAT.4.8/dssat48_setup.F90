!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: dssat48_setup
! \label{dssat48_setup}
!
! !REVISION HISTORY:
!  09 Apr 2024: Pang-Wei Liu 
! 
! !INTERFACE:
subroutine dssat48_setup()
! !USES:
   use netcdf
   use LIS_logMod,    only : LIS_verify, LIS_logunit, LIS_endrun

   use LIS_coreMod,   only : LIS_rc, LIS_surface, &
          LIS_localPet, LIS_ews_ind, LIS_ewe_ind,&
          LIS_nss_ind, LIS_nse_ind, & 
          LIS_ews_halo_ind, LIS_ewe_halo_ind,&
          LIS_nss_halo_ind, LIS_nse_halo_ind

   use LIS_fileIOMod, only : LIS_read_param
   use dssat48_lsmMod
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for DSSAT48.
!  
! The routines invoked are: 
! \begin{description}
! \item[dssat48\_setvegparms](\ref{dssat48_setvegparms}) \newline
!   initializes the vegetation-related parameters in dssat48
! \end{description}
!EOP

  implicit none
  integer           :: n
  integer           :: mtype
  integer           :: t
  integer           :: col, row
!  integer           :: nc, nr
!  integer           :: ews, ewe, nss, nse 
!  double precision  :: xmn_part  ! center x of local LL starting point
!  double precision  :: ymn_part  ! center y of local LL starting point
!  integer, allocatable :: global_kstn(:,:,:)
!  integer, allocatable :: level_data(:, :)
  integer, allocatable :: placeholder(:,:)
  CHARACTER*10   :: mukey
! _______________________________________________________________

  mtype = LIS_rc%lsm_index

  write(LIS_logunit,*)"[INFO] Reading in Crop Type and topography (dssat48_setup)"

  ! Read in spatial SnowModel maps from LDT:
  do n=1, LIS_rc%nnest

     ! Allocate memory for place holder for #n nest

     allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
     ! Read MUKEY from Input File
        write(LIS_logunit,*) "DSSAT48: reading MUKEY from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, 'MUKEY', placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            WRITE( mukey, '(I10)') placeholder(col, row) !Turn mukey number to a string
            dssat48_struc(n)%dssat48(t)%SLNO = 'LC'//repeat('0',8-len(trim(adjustl(mukey))))//adjustl(mukey)
        enddo

     !! Get Cropyears
     !  ! open NetCDF parameter file
     !   ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
     !   call LIS_verify(ios, 'Error in nf90_open in DSSAT48_read_cropyears')

     !   ! inquire the ID of east-west dimension
     !   ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
     !   call LIS_verify(ios, 'Error in nf90_inq_dimid in DSSAT48_read_MULTILEVEL_param')

     !   ! inquire the ID of north-south dimension
     !   ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
     !   call LIS_verify(ios, 'Error in nf90_inq_dimid in DSSAT48_read_MULTILEVEL_param')



        !write(LIS_logunit,*) "DSSAT48: reading CROPTYPE from ", trim(LIS_rc%paramfile(n))
        !do k = 1, 12 %1 to how many years
        !   call DSSAT48_read_MULTILEVEL_param(n, 'CROPTYPE', k, placeholder)
        !do t = 1, LIS_rc%npatch(n, mtype)
        !    col = LIS_surface(n, mtype)%tile(t)%col
        !    row = LIS_surface(n, mtype)%tile(t)%row
        !    dssat48_struc(n)%dssat48(t)%CDL(k) = placeholder(col, row)
        !enddo


     ! Initialize crop data layer and soil propertiy parameters:
     !snowmodel_struc(n)%sm(:)%smtopo = 0.
     !snowmodel_struc(n)%sm(:)%smvege = 0.

     ! Assign LDT-based CDL and soil maps :
     !if( dssat48_struc(n)%sm_params_opt == "LDT" ) then
     !   write(LIS_logunit,*) "[INFO] Reading in DSSAT48 parameters from LDT"
     !   ascii_topoveg = 2.0    ! No file read in by SnowModel; use LDT input
     !else
     !   write(LIS_logunit,*) "[INFO] Reading in SnowModel LSM parameters from snowmodel.par file "
     !endif

     !if( ascii_topoveg == 2.0 ) then
     !  write(LIS_logunit,*) "SnowModel: Reading parameter CDL from: ",trim(LIS_rc%paramfile(n))
     !  call LIS_read_param(n, "CDL", placeholder)
     !  do t = 1, LIS_rc%npatch(n, mtype)
     !     col = LIS_surface(n, mtype)%tile(t)%col
     !     row = LIS_surface(n, mtype)%tile(t)%row
     !     dssat48_struc(n)%sm(t)%smCDL = placeholder(col, row)
     !    ! Assign LDT-version of SnowModel topo map to topo_land ...
     !     cdl(col,row) = placeholder(col, row)
     !  enddo
  enddo

  ! ------------------------------ 
end subroutine dssat48_setup
 !BOP
  !
  ! !ROUTINE: CROCUS81_read_MULTILEVEL_param
  !  \label{read_MULTILEVEL_param}
  !
  ! !REVISION HISTORY:
  !  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
  !  30 Oct  2013: Shugong Wang; Generalization for reading MULTILEVEL spatial parameter
  !  22 Apr  2024: Pang-Wei Liu; reading multi year CDL

  ! !INTERFACE:
  subroutine DSSAT48_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
! !USES:
    use netcdf
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_localPet,   &
                            LIS_ews_halo_ind, LIS_ewe_halo_ind, &
                            LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    implicit none
! !ARGUMENTS: 
    integer, intent(in)          :: n
    integer, intent(in)          :: level
    character(len=*), intent(in) :: ncvar_name
    real, intent(out)            :: placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n))
! !DESCRIPTION:
!  This subroutine reads MULTILEVEL parameters from the LIS
!  NetCDF parameter data file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[level]
!    level index (month, quarter, soil layer, snow layer) of the data to be read
!   \item[array]
!    array containing returned values
!   \end{description}
!
!EOP      

    integer       :: ios1
    integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
    integer       :: nc, nr, t, nlevel, k
    real, pointer :: level_data(:, :, :)
    real, pointer :: level_data1(:, :, :)
    logical       :: file_exists

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) 'Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in DSSAT48_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in DSSAT48_read_MULTILEVEL_param')

        !! inquire the ID of north-south dimension !PL: This lines are redundant
        !ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        !call LIS_verify(ios, 'Error in nf90_inq_dimid in DSSAT48_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in DSSAT48_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in DSSAT48_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in DSSAT48_read_MULTILEVEL_param')

        ! inquire the ID of parameter. 
        ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! inquire the IDs of all dimensions. The third dimension is the level dimension
        ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')

        ! inquire the length of the level dimension
        ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire the length of the 3rd dimension')

        ! allocate memory
        allocate(level_data(LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in RUC37_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in RUC37_read_MULTILEVEL_param')

        ! grab parameter at specific level
        placeholder(:, :) = &
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        ! free memory 
        deallocate(level_data)

    else
        write(LIS_logunit, *) 'MULTILEVEL parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit, *) 'program stopping ...'
        call LIS_endrun
    endif
  end subroutine DSSAT48_read_MULTILEVEL_param
