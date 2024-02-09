!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module Crocus81_da_dhdt_Mod
!BOP
!
! !MODULE: Crocus81_da_dhdt_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! !USES:        

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: Crocus81_da_dhdt_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
! 
! !ROUTINE: Crocus81_da_dhdt_init
! \label{Crocus81_da_dhdt_init}
! 
! !INTERFACE:
  subroutine Crocus81_da_dhdt_init()
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none
  end subroutine Crocus81_da_dhdt_init
end module Crocus81_da_dhdt_Mod
