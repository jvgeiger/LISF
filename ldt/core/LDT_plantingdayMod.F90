!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_plantingdayMod
!BOP
!
! !MODULE: LDT_plantingdayMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read planting day map
!  for DSSAT model.
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the
!  planting day.
!
! !REVISION HISTORY:
!
!  24 Feb 2025: Pang-Wei Liu; Create for planting day data
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_plantingday_readParamSpecs
  public :: LDT_plantingday_init    !allocates memory for required structures
  public :: LDT_plantingday_writeHeader
  public :: LDT_plantingday_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_plantingday_struc

  type, public :: plantingday_type_dec

     character(len=LDT_CONST_PATH_LEN)    :: plantingdayfile
     real, allocatable :: plantingday_gridDesc(:)
     character*50     :: plantingday_proj
     character*50     :: plantingday_gridtransform
     ! Planting day parameters
     type(LDT_paramEntry) :: plantingday

  end type plantingday_type_dec

  type(plantingday_type_dec), allocatable :: LDT_plantingday_struc(:)

contains

  subroutine LDT_plantingday_readParamSpecs

    character*100     :: source
    integer           :: rc
    integer           :: n

    allocate(LDT_plantingday_struc(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config,"Planting day data source:",rc=rc)
    do n=1,LDT_rc%nnest
       allocate(LDT_plantingday_struc(n)%plantingday_gridDesc(20))
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_plantingday_struc(n)%plantingday,&
            "PLANTINGDAY",source)
    enddo

  end subroutine LDT_plantingday_readParamSpecs


!BOP
!
! !ROUTINE: LDT_plantingday_init
! \label{LDT_plantingday_init}
!
! !INTERFACE:
  subroutine LDT_plantingday_init
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_gridOptChecks
!    use LDT_logMod,    only : LDT_verify

! !DESCRIPTION:
!
! Allocates memory for data structures for reading
! the planting day datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[plantingdaysetup](\ref{plantingdaysetup}) \newline
!    calls the registry to invoke the drainage tile setup methods.
!  \end{description}
!
!EOP
    implicit none
    integer  :: n
    integer  :: k
    integer  :: rc
    type(LDT_fillopts)  :: plantingday

    logical  :: plantingday_select

! _____________________________________________________________

   plantingday_select = .false.

   do n = 1, LDT_rc%nnest
      if( LDT_plantingday_struc(n)%plantingday%selectOpt.gt.0 ) then
        plantingday_select = .true.
      endif
   enddo

   if(plantingday_select) then
     write(LDT_logunit,*)" - - - - - - - - - - Planting Day - - - - - - - - - - - - -"
   endif

   do n=1,LDT_rc%nnest

      if( plantingday_select ) then

         call set_plantingday_attribs( n, LDT_plantingday_struc(n)%plantingday%source )

         LDT_plantingday_struc(n)%plantingday%vlevels = LDT_plantingday_struc(n)%plantingday%num_bins
         allocate(LDT_plantingday_struc(n)%plantingday%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_plantingday_struc(n)%plantingday%vlevels))
      endif

   enddo

  ! Planting Day ldt.config entries:
    if( plantingday_select ) then

       call ESMF_ConfigFindLabel(LDT_config,"Planting day map:",rc=rc) !Directory of data
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_plantingday_struc(n)%plantingdayfile,rc=rc)
          call LDT_verify(rc,'Planting day map: not specified')
       enddo

       call ESMF_ConfigFindLabel(LDT_config,"Planting day spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_plantingday_struc(n)%plantingday_gridtransform,&
               rc=rc)
          call LDT_verify(rc,'Planting day spatial transform: option not specified in the config file')
       enddo

     ! Set units and full names:
       do n=1,LDT_rc%nnest
          LDT_plantingday_struc(n)%plantingday%units="DOY"
          call setPlantingParmsFullnames( n, "Planting Day", &
                  LDT_plantingday_struc(n)%plantingday%source )
       enddo

       ! Read in lat/lon extents and res for sources that are binary:
       do n=1,LDT_rc%nnest
           call ESMF_ConfigGetAttribute(LDT_config,LDT_plantingday_struc(n)%plantingday_proj,&
                 label="Planting day map projection:",rc=rc)
           call LDT_verify(rc,'Planting day map projection: option not specified in the config file')

           if( LDT_plantingday_struc(n)%plantingday_proj.NE."latlon" ) then
                 !Handle unsupported map projection.
                 write(LDT_logunit,*) &
                      '[ERR] Planting day map projection only supports latlon'
                 call LDT_endrun()
           endif
       enddo

    end if

!-- Read Planting day maps:
    do n = 1,LDT_rc%nnest
    !- Planting day:
       if( LDT_plantingday_struc(n)%plantingday%selectOpt == 1 ) then
          write(LDT_logunit,*) "Reading planting day map: "//&
               trim(LDT_plantingday_struc(n)%plantingdayfile)
          call readplantingday(&
               trim(LDT_plantingday_struc(n)%plantingday%source)//char(0),&
               n, LDT_plantingday_struc(n)%plantingday%num_bins, &
               LDT_plantingday_struc(n)%plantingday%value(:,:,1))
       end if
    enddo

  end subroutine LDT_plantingday_init

  subroutine LDT_plantingday_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer    :: n
    integer    :: ftn
    integer    :: dimID(3)
    integer    :: tdimID(3)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

    if( LDT_plantingday_struc(n)%plantingday%selectOpt.gt.0 ) then
       call LDT_verify(nf90_def_dim(ftn,'plantingday',&
            LDT_plantingday_struc(n)%plantingday%vlevels,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
                LDT_plantingday_struc(n)%plantingday)
    endif

#endif
  end subroutine LDT_plantingday_writeHeader

  subroutine LDT_plantingday_writeData(n,ftn)

    integer     :: n
    integer     :: ftn

    if( LDT_plantingday_struc(n)%plantingday%selectOpt.gt.0 ) then
       call LDT_writeNETCDFdata(n,ftn,LDT_plantingday_struc(n)%plantingday)
    endif

  end subroutine LDT_plantingday_writeData

end module LDT_plantingdayMod
