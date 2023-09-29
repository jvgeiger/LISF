!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module dssat48_module
!BOP
!
! !MODULE: dssat48_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the DSSAT model 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[n]
!     nest id. unit: -
!
! !REVISION HISTORY:
!  11 May 2023: Pang-Wei Liu
!  31 Aug 2023: J. Erlingis; Add forcing
!EOP
    implicit none
    private
    type, public :: dssat48dec
        !-------------------------------------------------------------------------
        ! forcing
        !-------------------------------------------------------------------------
        real               :: tair, tmax, tmin
        real               :: qair, swdown, lwdown
        real               :: uwind, vwind, wndspd
        real               :: psurf, rainf, snowf
        real               :: tdew, totprc
      
        !-------------------------------------------------------------------------
        ! spatial parameter
        !-------------------------------------------------------------------------
        !real               :: CDL
        !real               :: NUKEY
        !-------------------------------------------------------------------------
        ! multilevel spatial parameter
        !-------------------------------------------------------------------------
        !real, pointer      :: ALB(:) 
        !-------------------------------------------------------------------------
        ! state
        !-------------------------------------------------------------------------
        !real, pointer      :: SM(:)
        !real, pointer      :: SNOWRHO(:)
        !-------------------------------------------------------------------------
        ! output
        !-------------------------------------------------------------------------
        !real               :: SM
        !real               :: YIELD
        !real               :: LAI
        !-------------------------------------------------------------------------
        ! TIME CONTROL
        !-------------------------------------------------------------------------
         INTEGER            :: YREND, EXPNO, TRTALL
         INTEGER            :: NYRS, NREPS, YRSIM, ENDYRS, YRDOY_END
         CHARACTER*1        :: RNMODE
         CHARACTER*80       :: PATHEX
         CHARACTER*92       :: FILEX_P
         CHARACTER*120      :: FILECTL
    end type dssat48dec
          !INTEGER, PARAMETER :: &
          !Dynamic variable values
          !  RUNINIT  = 1, &
          !  INIT     = 2, & !Will take the place of RUNINIT & SEASINIT
          !               !     (not fully implemented)
          !  SEASINIT = 2, &
          !  RATE     = 3, &
          !  EMERG    = 3, & !Used for some plant processes.  
          !  INTEGR   = 4, & 
          !  OUTPUT   = 5, & 
          !  SEASEND  = 6, &
          !  ENDRUN   = 7    
end module dssat48_module
