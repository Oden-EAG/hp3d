!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> Purpose :  routine determines whether for a node case, all supported
!             variable components of a specific variable type with 
!             Dirichlet flags have homogeneous BC data
!! @param[in]  D_type  - discretization type
!! @param[in]  Icase   - encoded case number 
!! @param[in]  Bcond   - encoded nodal BC flag
!! @param[out] Is_homD = true/false
!! @param[out] Ncase   - decoded nodal case
!! @param[out] Ibcnd   - decoded BC flags

!-----------------------------------------------------------------------
!
      subroutine homogenD(D_type,Icase,Bcond, Is_homD,Ncase,Ibcnd)
!
      use physics
      implicit none
!
!  ...Arguments
      logical          :: Is_homD
      character(len=6) :: D_type
      integer          :: Icase, Bcond
!
!  ...decoded Icase oand Bcond
      integer, dimension(NR_PHYSA)  :: ncase
      integer, dimension(NRINDEX)   :: ibcnd
      integer :: iphys,ic,ivar
      logical          :: is_Dirichlet
!
!
      call decod(Icase,2,NR_PHYSA, Ncase)
      call decod(Bcond,2,NRINDEX,  Ibcnd)
      Is_homD = .true.
      ic=0
      do iphys = 1,NR_PHYSA
!
!  .....skip if the variable is not supported or the wrong type
        if ((Ncase(iphys).eq.0).or.(DTYPE(iphys).ne.D_type)) then
          ic = ic + NR_COMP(iphys)
          cycle
        endif
!
!  .....determine if a Dirichlet variable
        is_Dirichlet = .false.
        do ivar=1,NR_COMP(iphys)
          ic = ic+1
          if (ibcnd(ic).eq.1) is_Dirichlet = .true.
        enddo
!
!  .....check if a homogeneous Dirichlet BC variable
        if (is_Dirichlet.and..not.PHYSAd(iphys)) Is_homD = .false.
      enddo
!
      end subroutine homogenD