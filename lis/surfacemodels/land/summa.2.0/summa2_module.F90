!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: summa2_module.F90
module summa2_module
!
! !DESCRIPTION:
!  Declare forcing-only option (summa2) variables
!
!  \begin{description}
!   \item[forcing]
!     Array of meteorological forcing
!   \item{tair}
!     2m air temperature forcing
!   \item{qair}
!     2m specific humidity forcing
!   \item{swdown}
!     downward shortwave forcing
!   \item{lwdown}
!     downward longwave forcing
!   \item{uwind}
!     u-wind component forcing
!   \item{vwind}
!     v-wind component forcing
!   \item{psurf}
!     surface pressure forcing
!   \item{rainf}
!     total rainfall forcing
!   \item{rainf\_c}
!     convective rainfall forcing
!   \item{snowf}
!     total snowfall forcing
!   \end{description}
!
!EOP
   use nrtype, only : i4b
   use data_types

   implicit none


   type summa2dec

      ! define indices
      integer(i4b)     :: iVar           ! index of a model variable
      integer(i4b)     :: iStruct        ! loop through data structures
      integer(i4b)     :: iGRU
      integer(i4b)     :: iHRU,jHRU,kHRU ! index of the hydrologic response unit

      ! summa2-Forcing Variables
      real :: tair
      real :: qair
      real :: swdown
      real :: lwdown
      real :: uwind
      real :: vwind
      real :: psurf
      real :: rainf
      real :: rainf_c
      real :: snowf

   end type summa2dec

end module summa2_module
