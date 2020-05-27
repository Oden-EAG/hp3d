!--------------------------------------------------------------------
!
!     routine name      - elem_residual
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 2019
!
!     purpose:          - routine returns element residual (squared)
!                         for the Primal Poisson and UW Time Harmonic
!                         Maxwell equation
!
!     arguments:
!        in:
!             Mdle      - element middle node number
!        out:
!             Resid     - element residual (squared)
!             Nref_flag - suggested h-refinement flag
!
!---------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem_residual(Mdle, Resid,Nref_flag)
!..modules used
   use commonParam
   use control
   use data_structure3D
   use element_data
   use parametersDPG
!..no implicit statements
   implicit none
!..declare input/output variables
   integer, intent(in)  :: Mdle
   real(8), intent(out) :: Resid
   integer, intent(out) :: Nref_flag
!
!..number of test and trial degrees of freedom
   integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ
   integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
   integer :: nrTest
!
!..(enriched) order of element nodes
   integer :: norder(19),norderP(19),nordP
!
!..element type
   character(len=4) :: etype
!
!..fld_flag refers to either pump (0) or signal (1) field
   integer :: fld_flag
!
!---------------------------------------------------------------------
!
   norder (1:19) = 0
   norderP(1:19) = 0
!
   etype = NODES(Mdle)%type
!..determine order of approximation
   call find_order(Mdle, norder)
!..set the enriched order of appoximation
   select case(etype)
      case('mdlb')
         nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdln','mdld')
         nordP = NODES(Mdle)%order+NORD_ADD
      case('mdlp')
         nordP = NODES(Mdle)%order+NORD_ADD*11
      case default
         write(*,*) 'elem_residual: invalid etype param. stop.'
         stop
   end select
!..note: compute_enriched_order works only for hexa and prism currently
   call compute_enriched_order(etype,nordP, norderP)
!..compute nrdof for trial
   call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!..compute nrdof for test
   call celndof(etype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
   select case(NO_PROBLEM)
      case(1,2)
         nrTest = nrdofHH
         call elem_residual_heat(Mdle,             &
            nrTest,nrdofHH,nrdofH,nrdofV,          &
            Resid,Nref_flag)
      case(3,4)
         nrTest = 2*nrdofEE
         if (NO_PROBLEM .eq. 3) then
            fld_flag = 1
         else
            fld_flag = 0
         endif
         call elem_residual_maxwell(Mdle,fld_flag,    &
            nrTest,nrdofEE,nrdofH,nrdofE,nrdofQ,      &
            Resid,Nref_flag)
      case default
         write(*,*) 'elem_residual: invalid NO_PROBLEM param. stop.'
         stop
   end select
!
end subroutine elem_residual
