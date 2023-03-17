!======================================================================
! REMARK: unused routine!!!
c----------------------------------------------------------------------
c
c   routine name       - curve_1SurfrCur
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine supports parametrization of a curve
c                        lying on a single reconstructed surface
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
      subroutine curve_1SurfrCur(No,Eta, X,Dxdeta)
c
      use GMP
      implicit none
c
      integer :: No
      real(8) :: Eta,X(NDIM),Dxdeta(NDIM)
c
c  ...endpoints coordinates
      real(8) :: xv(NDIM,2)
c
      integer :: i,j,np
c
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
      end subroutine curve_1SurfrCur
