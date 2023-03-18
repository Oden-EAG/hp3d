!======================================================================
! REMARK: unused routine!!!
c----------------------------------------------------------------------
c
c   routine name       - curve_2SurfrCur
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 07
c
c   purpose            - routine supports parametrization of a curve
c                        lying on intersection of two reconstructed
c                        surfaces
c
c   arguments :
c     in:
c               No     - the curve number
c               Eta    - reference coordinate  (between 0 and 1)
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        reference coordinate
c
c----------------------------------------------------------------------
c
      subroutine curve_2SurfrCur(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      integer :: No
      real(8) :: X(NDIM),Dxdeta(NDIM)
c
c  ...endpoints coordinates
      real(8) :: xv(NDIM,2)
c
      integer :: i
c
c  ...get the curve endpoints
      do i=1,2
        np = CURVES(No)%EndPoNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
c
c  ...use segment of straight line for now...
      do j=1,NDIM
        X(j) = Eta*xv(j,2) + (1.d0-Eta)*xv(j,1)
        Dxdeta(j) = xv(j,2) - xv(j,1)
      enddo
c
c
      end subroutine curve_2SurfrCur
