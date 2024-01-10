!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief Routine determines whether for a node case, all supported
!!        variable components of a specific variable type with
!!        Dirichlet flags have homogeneous BC data
!! @param[in]  Dtype   - discretization type
!! @param[in]  Icase   - encoded nodal case number
!! @param[in]  Bcond   - encoded nodal BC flag
!! @param[out] Is_homD - true/false
!! @param[out] Ncase   - decoded nodal case numbers
!! @param[out] Ibcnd   - decoded nodal BC flags
!
!> @date Mar 2023
!-----------------------------------------------------------------------
!
subroutine homogenD(Dtype,Icase,Bcond, Is_homD,Ncase,Ibcnd)
!
   use physics
!
   implicit none
!
!..input arguments
   integer, intent(in)  :: Dtype,Icase,Bcond
!
!..output arguments
   logical, intent(out) :: Is_homD
   integer, intent(out) :: Ncase(NR_PHYSA)
   integer, intent(out) :: Ibcnd(NRINDEX_HEV)
!
   logical :: is_Dirichlet
!
   integer :: iphys,ic,ivar
!
!-----------------------------------------------------------------------
!
   call decod(Icase,2,NR_PHYSA, Ncase)
   call decod(Bcond,2,NRINDEX_HEV,  Ibcnd)
   Is_homD = .true.
   ic=0
   do iphys = 1,NR_PHYSA
!
!  ...skip if the variable is not supported or the wrong type
      if ((Ncase(iphys).eq.0) .or. (D_TYPE(iphys).ne.Dtype)) then
         ic = ic + NR_COMP(iphys)
         cycle
      endif
!
!  ...determine if a Dirichlet variable
      is_Dirichlet = .false.
      do ivar=1,NR_COMP(iphys)
         ic = ic+1
         if (Ibcnd(ic).eq.1) is_Dirichlet = .true.
      enddo
!
!  ...check if a homogeneous Dirichlet BC variable
      if (is_Dirichlet .and. (.not.PHYSAd(iphys))) Is_homD = .false.
   enddo
!
end subroutine homogenD
