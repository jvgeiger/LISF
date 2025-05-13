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
! !ROUTINE: setPlantingParmsFullnames
!  \label{setPlantingParmsFullnames}
!
! !REVISION HISTORY:
!  26 Feb 2025: P.W. Liu; Support for planting day map
!
! !INTERFACE:
subroutine setPlantingParmsFullnames(n,datatype,source)

! !USES:
  use LDT_plantingdayMod
  use LDT_logMod
  use LDT_paramDataMod

  implicit none
  integer,         intent(in) :: n
  character(len=*),intent(in) :: datatype
  character(len=*),intent(in) :: source

! !ARGUMENTS:

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[datatype]
!     Planting day map data type
!   \item[source]
!     Planting day map dataset source
!   \end{description}
!EOP
!

   select case( datatype )

    case( "Planting Day" )
    select case( source )
        case( "CREATE" )
          LDT_plantingday_struc(n)%plantingday%standard_name =&
             "Planting day by create"
      end select

    case default
      write(LDT_logunit,*) "[ERR] Planting day type not recognized: ",trim(source)
      write(LDT_logunit,*) " Program stopping ..."
      stop
   end select

end subroutine setPlantingParmsFullnames
