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
! !ROUTINE: Crocus81_descale_dhdt
! \label{Crocus81_descale_dhdt}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 8 Jan  2024 : Mahdi Navari ; Modified for ICESat2 dhdt da 
!
! !INTERFACE:
subroutine Crocus81_descale_dhdt(n, LSM_State, LSM_Incr_State, LIS_LSM_particle_weight)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noahmp401_lsmMod
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
  type(ESMF_State)       :: LIS_LSM_particle_weight
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  !type(ESMF_Field)       :: snodField 
  !type(ESMF_Field)       :: ParticleWeightField

  !integer                :: t
  !integer                :: status
  !real, pointer          :: snod(:)
  !real, pointer          :: pw(:)
end subroutine Crocus81_descale_dhdt

