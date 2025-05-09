!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: SMscaling_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle soil moisture bias correction b/w LIS and DSSAT
! 
! !REVISION HISTORY: 
!  22 Aug 2016    Sujay Kumar; initial specification
!  09 Oct 2024: Pang-Wei Liu; bias correction for coupling LIS and DSSAT
!
module SMscaling_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_LMLCMod
!---
  use netcdf
  use LIS_coreMod
  use LIS_logmod
  use LIS_histDataMod
  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: SMscaling_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: SMscaling_struc
!EOP
  type, public:: SMscaling_dec
     
     real,    allocatable       :: dssat_xrange(:,:,:,:)
     real,    allocatable       :: lsm_xrange(:,:,:,:)
     real,    allocatable       :: dssat_cdf(:,:,:,:)
     real,    allocatable       :: lsm_cdf(:,:,:,:)
     real,    allocatable       :: dssat_mu(:,:,:)
     real,    allocatable       :: lsm_mu(:,:,:) 
     real,    allocatable       :: dssat_sigma(:,:,:)
     real,    allocatable       :: lsm_sigma(:,:,:)
     real,    allocatable       :: landmaskind(:,:)
     real,    allocatable       :: bc_glbngrid_red(:)
     integer                :: ntimes
     logical                :: cdf_read_mon  !(for reading monthly CDF when
                                             !LIS_rc%da > 1 but the first model time step,
                                             !e.g., 4/29 13:00:00)
     integer                :: cdf_read_opt  ! 0: read all months at one time
                                             ! 1: read only the current month
     character*20             :: bcmethod
     integer                :: bclayer
     integer                :: nbins
     character(len=LIS_CONST_PATH_LEN) :: dssatcdffiles(4)
     character(len=LIS_CONST_PATH_LEN) :: lsmcdffiles(4)

  end type SMscaling_dec
  
  type(SMscaling_dec),allocatable :: SMscaling_struc(:)
  
contains

!BOP
! 
! !ROUTINE: SMscaling_setup
! \label{SMscaling_setup}
! 
! !INTERFACE: 
  subroutine SMscaling_setup
! !USES: 
    use LIS_dataAssimMod, only : LIS_getCDFattributes 
    implicit none 
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data strctures required for handling LSM and DSSAT
!   soil moisture data. 
!  
!EOP
    integer                ::  n, k, i, l !l: layer of soil
    integer                ::  rc, ngrid
    integer                ::  r, c

    allocate(SMscaling_struc(LIS_rc%nnest))

    call ESMF_ConfigFindLabel(LIS_config,"Bias correction scaling strategy:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, SMscaling_struc(n)%bcmethod, rc = rc)
       call LIS_verify(rc,'Bias correction scaling strategy: not defined')
       write(LIS_logunit,*) "[INFO] DSSAT SM strategy ", SMscaling_struc(n)%bcmethod
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"Numbers of soil layers for bias correction:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, SMscaling_struc(n)%bclayer, rc = rc)
       call LIS_verify(rc,'Numbers of soil layers for bias correction: not defined')
       write(LIS_logunit,*) "[INFO] DSSAT SM number of layer ", SMscaling_struc(n)%bclayer
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"DSSAT model CDF files:",rc=rc)
    do n=1,LIS_rc%nnest
       do i=1,SMscaling_struc(n)%bclayer
          call ESMF_ConfigGetAttribute(LIS_config, SMscaling_struc(n)%dssatcdffiles(i), rc = rc)
          call LIS_verify(rc,'DSSAT model CDF files: not defined')
          write(LIS_logunit,*) "[INFO] DSSAT SM file ", trim(SMscaling_struc(n)%dssatcdffiles(i))
       enddo
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"LSM model CDF files:",rc=rc)
    do n=1,LIS_rc%nnest
       do i=1,SMscaling_struc(n)%bclayer
          call ESMF_ConfigGetAttribute(LIS_config, SMscaling_struc(n)%lsmcdffiles(i), rc = rc)
          call LIS_verify(rc,'LSM model CDF files: not defined')
             write(LIS_logunit,*) "[INFO] LSM SM file ", trim(SMscaling_struc(n)%lsmcdffiles(i))
       enddo
    enddo

     call ESMF_ConfigFindLabel(LIS_config,"LIS-DSSAT number of bins in the CDF:",rc=rc)
     do n=1,LIS_rc%nnest
           call ESMF_ConfigGetAttribute(LIS_config, SMscaling_struc(n)%nbins, rc = rc)
           call LIS_verify(rc,'LIS-DSSAT number of bins in the CDF: not defined')
           write(LIS_logunit,*) "[INFO] Bin number ", SMscaling_struc(n)%nbins
     enddo

    write(LIS_logunit,*)&
         '[INFO] read LIS-DSSAT bias correction specifications'       

!----------------------------------------------------------------------------
!   Create the array containers that will contain the soil moisture statistics 
!   of LSM and DSSAT in the grid space. Since there is only upto four layers
!   being bias corrected, the array size is LIS_rc%ngrid(k). 
!   
!----------------------------------------------------------------------------
    k = 1 !k is always 1 in this case
    do n=1,LIS_rc%nnest
       if(SMscaling_struc(n)%bcmethod.eq."CDF matching" .or. &
               SMscaling_struc(n)%bcmethod.eq."Linear scaling" .or. &
               SMscaling_struc(n)%bcmethod.eq."Anomaly scaling") then

          call LIS_getCDFattributes(k,SMscaling_struc(n)%dssatcdffiles(1),&
               SMscaling_struc(n)%ntimes,ngrid)                 

!          if (SMscaling_struc(n)%cdf_read_opt.eq.0) then   !Leave it open for future read method option
!------------------- TURN CDF OFF ---------------------------------------------
!             allocate(SMscaling_struc(n)%dssat_xrange(&
!                  SMscaling_struc(n)%bclayer,              &        
!                  LIS_rc%ngrid(1), SMscaling_struc(n)%ntimes, &
!                  SMscaling_struc(n)%nbins))

!             allocate(SMscaling_struc(n)%lsm_xrange(&
!                   SMscaling_struc(n)%bclayer,              &
!                  LIS_rc%ngrid(1), SMscaling_struc(n)%ntimes, &
!                  SMscaling_struc(n)%nbins))

!             allocate(SMscaling_struc(n)%dssat_cdf(&  
!                   SMscaling_struc(n)%bclayer,              &    
!                  LIS_rc%ngrid(1), SMscaling_struc(n)%ntimes, &
!                  SMscaling_struc(n)%nbins))

!             allocate(SMscaling_struc(n)%lsm_cdf(&
!                   SMscaling_struc(n)%bclayer,              &
!                  LIS_rc%ngrid(1), SMscaling_struc(n)%ntimes, &
!                  SMscaling_struc(n)%nbins))
!---------------------------------------------------------------------------
             !Dimension of mu & sigma: layers(upto 4 layers), ngrid(number of local pixels), ntimes(12)
             allocate(SMscaling_struc(n)%dssat_mu(&
                     SMscaling_struc(n)%bclayer,              &
                     LIS_rc%ngrid(k), SMscaling_struc(n)%ntimes))
             
             allocate(SMscaling_struc(n)%dssat_sigma(&
                  SMscaling_struc(n)%bclayer,              &
                  LIS_rc%ngrid(k),SMscaling_struc(n)%ntimes))

             allocate(SMscaling_struc(n)%lsm_mu(&
                  SMscaling_struc(n)%bclayer,              &
                  LIS_rc%ngrid(k), SMscaling_struc(n)%ntimes))

             allocate(SMscaling_struc(n)%lsm_sigma(&
                  SMscaling_struc(n)%bclayer,              &
                  LIS_rc%ngrid(k),SMscaling_struc(n)%ntimes))

             allocate(SMscaling_struc(n)%bc_glbngrid_red(k))
             allocate(SMscaling_struc(n)%landmaskind(LIS_rc%gnc(n),&
                  LIS_rc%gnr(n))) !Allocate Space to landmaskind

             SMscaling_struc(n)%landmaskind = LIS_LMLC(n)%glandmask

             SMscaling_struc(n)%bc_glbngrid_red(k) = 0
             do r=1,LIS_rc%gnr(k)
                do c=1,LIS_rc%gnc(k)
                   if(SMscaling_struc(n)%landmaskind(c,r).gt.0) then
                      SMscaling_struc(n)%bc_glbngrid_red(k) =&
                           SMscaling_struc(n)%bc_glbngrid_red(k) + 1
                      SMscaling_struc(n)%landmaskind(c,r) = &
                           SMscaling_struc(n)%bc_glbngrid_red(k)   
                   endif
                enddo 
             enddo

!---------------------------------------------------------------------------
!          else  
!             allocate(SMscaling_struc(n)%model_xrange(&
!                  LIS_rc%obs_ngrid(k), 1, &
!                  SMscaling_struc(n)%nbins))
!             allocate(SMscaling_struc(n)%obs_xrange(&
!                  LIS_rc%obs_ngrid(k), 1, &
!                  SMscaling_struc(n)%nbins))
!             allocate(SMscaling_struc(n)%model_cdf(&
!                  LIS_rc%obs_ngrid(k), 1, &
!                  SMscaling_struc(n)%nbins))
!             allocate(SMscaling_struc(n)%obs_cdf(&
!                  LIS_rc%obs_ngrid(k), 1, &
!                  SMscaling_struc(n)%nbins))
!             allocate(SMscaling_struc(n)%model_mu(LIS_rc%obs_ngrid(k),1))
!             allocate(SMscaling_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),1))
!             allocate(SMscaling_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),1))
!             allocate(SMscaling_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),1)) 
!          endif 
!----------------------------------------------------------------------------
! Read the model and observation CDF data
!----------------------------------------------------------------------------
!         if (SMscaling_struc(n)%cdf_read_opt.eq.0) then   

        DO l=1, SMscaling_struc(n)%bclayer

          call read_MeanSigmaData_all_dssat(n,k,&      !Read mu and sigma
               SMscaling_struc(n)%ntimes,&
               LIS_rc%ngrid(k), &
               SMscaling_struc(n)%dssatcdffiles(l), &
               "SoilMoist",&
               SMscaling_struc(n)%dssat_mu(l,:,:),&
               SMscaling_struc(n)%dssat_sigma(l,:,:))

          call read_MeanSigmaData_all_dssat(n,k,&
               SMscaling_struc(n)%ntimes,&
               LIS_rc%ngrid(k), &
               SMscaling_struc(n)%lsmcdffiles(l), &
               "SoilMoist",&
               SMscaling_struc(n)%lsm_mu(l,:,:),&
               SMscaling_struc(n)%lsm_sigma(l,:,:))
        ENDDO
!-----------TURN CDF OFF------------------------------------------------        
       !DO l=1, SMscaling_struc(n)%bclayer
          !call LIS_readCDFdata(n,k,&         !Soil layer 1
          !     SMscaling_struc(n)%nbins,&
          !     SMscaling_struc(n)%ntimes,&
          !     LIS_rc%ngrid, &
          !     SMscaling_struc(n)%dssatcdffiles(1), &
          !     "SoilMoist",&
          !     SMscaling_struc(n)%dssat_xrange,&
          !     SMscaling_struc(n)%dssat_cdf)    
          
          !call LIS_readCDFdata(n,k,&
          !     SMscaling_struc(n)%nbins,&
          !     SMscaling_struc(n)%ntimes,&
          !     LIS_rc%ngrid, &
          !     SMscaling_struc(n)%lsmcdffiles(1), &
          !     "SoilMoist",&
          !     SMscaling_struc(n)%lsm_xrange,&
          !     SMscaling_struc(n)%lsm_cdf)
       !ENDDO      
       endif
    enddo
!------------------------------------------------------------------------------
  end subroutine SMscaling_setup


!BOP
! !ROUTINE: read_MeanSigmaData_all_dssat
! \label{read_MeanSigmaData_all}
!
! !INTERFACE: 
   subroutine read_MeanSigmaData_all_dssat(n, k, ntimes, ngrid, filename, varname, mu, sigma)

     implicit none
! !ARGUMENTS:      
     integer,   intent(in)    :: n
     integer,   intent(in)    :: k
     integer,   intent(in)    :: ntimes
     integer,   intent(in)    :: ngrid
     character(len=*)         :: filename
     character(len=*)         :: varname
     real                     :: mu(ngrid, ntimes)
     real                     :: sigma(ngrid, ntimes)
! 
! !DESCRIPTION: 
!  This routine reads the input CDF file (generated by LDT in NETCDF format)
!  and reads the mean (mu) and standard deviation (sigma) values for each
!  grid point. 
!  Both these fields are expected to be in the 1-d grid vector dimension. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]             index of the nest
!  \item[filename]      name of the CDF file
!  \item[varname]       name of the variable being extracted.
!  \item[mu]            mean values
!  \item[sigma]         standard deviation values
! \end{description}
!EOP
     integer                  :: j
     integer                  :: nlevsId, gId
     integer                  :: ngrid_file, nlevs_file
     integer                  :: muid, sigmaid
     real, allocatable        :: mu_file(:,:,:)
     real, allocatable        :: sigma_file(:,:,:)
     integer                  :: nid

     write(LIS_logunit,*) '[INFO] Reading mean and sigma from CDF file ',trim(filename)
     call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
          ncid=nid),'failed to open file '//trim(filename))

     call LIS_verify(nf90_inq_dimid(nid,trim(varname)//"_levels",nlevsId), &
          'nf90_inq_dimid failed for '//trim(varname)//"_levels")

     call LIS_verify(nf90_inquire_dimension(nid,nlevsId, len=nlevs_file),&
          'nf90_inquire_dimension failed for nbins')

     call LIS_verify(nf90_inq_dimid(nid, 'ngrid',gId), &
          'Error nf90_inq_dimid: ngrid')

     call LIS_verify(nf90_inquire_dimension(nid, gId, len=ngrid_file), &
          'Error nf90_inquire_dimension:ngrid')

     allocate(mu_file(ngrid_file,ntimes, nlevs_file))
     allocate(sigma_file(ngrid_file,ntimes, nlevs_file))

     call LIS_verify(nf90_inq_varid(nid,trim(varname)//'_mu',muid),&
          'nf90_inq_varid failed for for '//trim(varname)//'_mu')
     call LIS_verify(nf90_inq_varid(nid,trim(varname)//'_sigma',sigmaid),&
          'nf90_inq_varid failed for '//trim(varname)//'_sigma')
     call LIS_verify(nf90_get_var(nid,muid,mu_file),&
          'nf90_get_var failed for '//trim(varname)//'_mu')
     call LIS_verify(nf90_get_var(nid,sigmaid,sigma_file),&
          'nf90_get_var failed for '//trim(varname)//'_sigma')

     if (LIS_rc%ngrid(k).gt.0) then !Conver the mean and sigma data from global domain to local
        do j=1,ntimes
          call LIS_convertBCVarToLocalSpace(n,k,mu_file(:,j,1),mu(:,j))
          call LIS_convertBCVarToLocalSpace(n,k,sigma_file(:,j,1),sigma(:,j))
        enddo
     endif

     !mu(:,:)    = mu_file(:,:,1)
     !sigma(:,:) = sigma_file(:,:,1)

     deallocate(mu_file)
     deallocate(sigma_file)

     call LIS_verify(nf90_close(nid),&
          'failed to close file '//trim(filename))
     write(LIS_logunit,*) '[INFO] Successfully read mean and sigma from CDF file ',trim(filename)

   end subroutine read_MeanSigmaData_all_dssat

!BOP
! !ROUTINE: LIS_convertBCVarToLocalSpace
! \label{LIS_convertBCVarToLocalSpace}
!  
! !INTERFACE:
  subroutine LIS_convertBCVarToLocalSpace(n,k, gvar,lvar)
! !USES:
      
    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: n
    integer, intent(in)   :: k
    real                  :: gvar(LIS_rc%glbngrid_red(n)) !Size should be the global "land" pixel size
    real                  :: lvar(LIS_rc%ngrid(k))
! !DESCRIPTION:
!  Converts a variable in the global 1-d space to local 1-d DAspace. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \item[k]
!    index of the instance (LSM), this is always 1 here
!   \item [gvar]
!    input variable dimensioned in the 1-d global space
!   \item [lvar]
!    output variable dimensioned in the 1-d local space
!  \end{description}
!EOP
    integer             :: count1
    integer             :: c,r,t
    !integer             :: gid, tid, ntiles, stid
    integer             :: lmaskind
    !----- Code similar to LIS_DAobservationdsMod.F90 ----------------------------------
    !do t=1, LIS_rc%ngrid(k)
    !   c = LIS_domain(n,k)%col(t) + LIS_ews_halo_ind(k,LIS_localPet+1)-1
       !r = LIS_domain(n,k)%row(t) + LIS_nss_halo_ind(k,LIS_localPet+1)-1
       !lmask = int(SMscaling_struc(n)%landmaskind(c,r))
       !lmask = int(LIS_obs_domain(n,k)%landmask(c,r)) !Remove this line
       !if(lmask.gt.0) then
       !   lvar(t) = gvar(lmask)
       !endif
    !enddo
    !-------------------------------------------------------------------------------------
    !----- Using codes combining concepts of that from LIS_DAobservationdsMod.F90 & LIS_historyMod.F90
    count1 = 1
    do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
          lmaskind = int(SMscaling_struc(n)%landmaskind(c,r))
          if(lmaskind.gt.0) then
             lvar(count1) = gvar(lmaskind)
             count1 = count1+1
          endif
       enddo
    enddo

    !------ Similar Code similar to LIS_historyMod.F90 -----------------------------------------
    !------ Can we do this way? Check with Sujay -----------------------------------------
    !count1=1
    !do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
    !   do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)

    !      gid = c+(r-1)*LIS_rc%gnc(n)
    !      ntiles = LIS_domain(n)%ntiles_pergrid(gid)
    !      stid = LIS_domain(n)%str_tind(gid)
    !      do t=1,ntiles
    !         tid = stid + t-1
    !         var(count1) = gtmp(tid)
    !         count1 = count1 + 1
    !      enddo
    !   enddo
    !enddo
    !PRINT*,'count1: ', count1
    !------------------------------------------------------------------------------------
end subroutine LIS_convertBCVarToLocalSpace
end module SMscaling_Mod
