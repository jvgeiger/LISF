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
   use LIS_logMod,    only : LIS_verify, LIS_logunit, LIS_endrun

   use LIS_coreMod,   only : LIS_rc, LIS_surface
!         LIS_localPet, LIS_ews_ind, LIS_ewe_ind,&
!         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind, LIS_ewe_halo_ind,&
!         LIS_nss_halo_ind, LIS_nse_halo_ind

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
!  integer           :: ews, ewe, nss, nse 
!  double precision  :: xmn_part  ! center x of local LL starting point
!  double precision  :: ymn_part  ! center y of local LL starting point
!  integer, allocatable :: global_kstn(:,:,:)
  real, allocatable :: placeholder(:,:)

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

  enddo

  ! ------------------------------ 
end subroutine dssat48_setup
 
