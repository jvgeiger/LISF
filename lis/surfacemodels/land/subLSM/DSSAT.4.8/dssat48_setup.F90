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
!  26 Jun 2023: Pang-Wei Liu 
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
  integer       :: ios, nid, param_ID, nc_ID, nr_ID
  integer       :: nc, nr
!  integer           :: ews, ewe, nss, nse 
!  double precision  :: xmn_part  ! center x of local LL starting point
!  double precision  :: ymn_part  ! center y of local LL starting point
!  integer, allocatable :: global_kstn(:,:,:)
  integer, allocatable :: level_data(:, :)
  integer, allocatable :: placeholder(:,:)
  CHARACTER*10   :: mukey
! _______________________________________________________________

  mtype = LIS_rc%lsm_index

  write(LIS_logunit,*)"[INFO] Reading in Crop Type and topography (dssat48_setup)"

  ! Read in spatial SnowModel maps from LDT:
  do n=1, LIS_rc%nnest

     ! Allocate memory for place holder for #n nest
     allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))


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

!       write(LIS_logunit,*) "SnowModel: Reading parameter SOIL from: ",&
!          trim(LIS_rc%paramfile(n))
!       call LIS_read_param(n, "SOIL", placeholder)
!       do t = 1, LIS_rc%npatch(n, mtype)
!          col = LIS_surface(n, mtype)%tile(t)%col
!          row = LIS_surface(n, mtype)%tile(t)%row
!          dssat48_struc(n)%sm(t)%smsoil = placeholder(col, row)
!         ! Assign LDT-version of dssat48 soil types ...
!          soiltype(col,row) = placeholder(col, row)
!       enddo
!
!     endif


       !----------------------------------------------!
        ! MULTILEVEL reading spatial spatial parameters Pang 2024.02.09!
        !----------------------------------------------!
          ! read: Mukey from external file
          ! write(LIS_logunit,*) "DSSAT48: reading parameter MUKEY  from ", trim('lis_input_IAcounty_merged.nc')
          !CALL DSSAT48_read_EXTERNAL_param(n, 'MUKEY', placeholder)
          ! open NetCDF parameter file
          ios = nf90_open(path=trim('lis_input_IAcounty_merged.nc'), mode=NF90_NOWRITE, ncid=nid)
          call LIS_verify(ios, 'Error in nf90_open in dssat48_setup')

          ! inquire the ID of east-west dimension
          ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
          call LIS_verify(ios, 'Error in nf90_inq_dimid in dssat48_setup')

          ! inquire the ID of north-south dimension
          ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
          call LIS_verify(ios, 'Error in nf90_inq_dimid in dssat48_setup')
          ! inquire the length of east-west dimension

          ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
          call LIS_verify(ios, 'Error in nf90_inquire_dimension in dssat48_setup')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in dssat48_setup')

        ! inquire the ID of parameter. 
        ios = nf90_inq_varid(nid, Trim('MUKEY'), param_ID)
        call LIS_verify(ios, trim('MUKEY')//' field not found in the external file')

        ! allocate memory
        allocate(level_data(LIS_rc%gnc(n), LIS_rc%gnr(n)))

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in dssat48_setup')
      
        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in dssat48_setup')

        ! grab parameter at specific level
        placeholder(:, :) = &
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1))
        ! free memory 
        deallocate(level_data)
        do t = 1, LIS_rc%npatch(n, mtype)
                  col = LIS_surface(n, mtype)%tile(t)%col
                  row = LIS_surface(n, mtype)%tile(t)%row
                  WRITE( mukey, '(I10)') placeholder(col, row) !Turn mukey number to a string
                  !PRINT*, 'In dsst48_setup: t, mukey', t, trim(mukey)
                  dssat48_struc(n)%dssat48(t)%SLNO = 'LC'//repeat('0',8-len(trim(adjustl(mukey))))//adjustl(mukey)
                  !PRINT*, 'dssat48_struc(n)%dssat48(t)%SLNO ',dssat48_struc(n)%dssat48(t)%SLNO
        enddo

  enddo

  ! ------------------------------ 
end subroutine dssat48_setup
 
