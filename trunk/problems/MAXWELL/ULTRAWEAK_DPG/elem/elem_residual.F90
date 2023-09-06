!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!> @brief       Compute (squared) norm of element residual
!!
!> @param[in]   mdle       - middle node number
!> @param[out]  Resid      - element residual (norm squared)
!> @param[out]  Nref_flag  - suggested  h-refinement refinement
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine elem_residual(Mdle, Resid,Nref_flag)
!
      use commonParam
      use control
      use data_structure3D
      use element_data
      use parametersDPG
!
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(out) :: Resid
      integer, intent(out) :: Nref_flag
!
!  ...number of test and trial degrees of freedom
      integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ
      integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
      integer :: nrTest
!
!  ...(enriched) order of element nodes
      integer :: norder(19),norderP(19),nordP
      integer :: ntype
!
!  ...fld_flag refers to either pump (0) or signal (1) field
      integer :: fld_flag
!
!---------------------------------------------------------------------
!
      norder (1:19) = 0
      norderP(1:19) = 0
!
!  ...get mdle node type
      ntype = NODES(Mdle)%ntype
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
!
!  ...set the enriched order of appoximation
      select case(ntype)
         case(MDLB)
            nordP = NODES(Mdle)%order+NORD_ADD*111
         case(MDLP)
            nordP = NODES(Mdle)%order+NORD_ADD*11
         case(MDLN,MDLD)
            nordP = NODES(Mdle)%order+NORD_ADD
         case default
            write(*,*) 'elem_residual: invalid ntype param. stop.'
            stop
      end select
!
!  ...note: compute_enriched_order works only for hexa and prism currently
      call compute_enriched_order(ntype,nordP, norderP)
!  ...compute nrdof for trial
      call celndof(ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!  ...compute nrdof for test
      call celndof(ntype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
      nrTest = 2*nrdofEE
!
      call elem_residual_maxwell(Mdle,nrTest,                       &
                                 nrdofEE,nrdofH,nrdofE,nrdofQ,      &
                                 Resid,Nref_flag)
!
end subroutine elem_residual
