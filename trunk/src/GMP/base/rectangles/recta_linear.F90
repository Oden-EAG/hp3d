!------------------------------------------------------------------------
!
!   routine name       - recta
!
!------------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine evaluates physical coordinates
!                        and its derivatives wrt to reference
!                        coordinates for a point in a GMP rectangle
!                        USING A BILINEAR APPROXIMATION
!
!   arguments :
!     in:
!               No     - the rectangle number
!               Eta    - reference coordinates of a point
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the reference coordinates
!
!------------------------------------------------------------------------
!
   subroutine recta_linear(No,Eta, X,Dxdeta)
!
      use GMP
!
      implicit none
!
      integer :: No
      real(8) :: Eta(2),X(3),Dxdeta(NDIM,2)
!
!  ...vertex coordinates, work space
      real(8) :: xv(NDIM,4),aux(NDIM,3)
!
      integer :: i,j,np
!
      integer :: iprint
!------------------------------------------------------------------------
!
      iprint=0
!
!  ...evaluate vertex coordinates
      do i=1,4
        np = RECTANGLES(No)%VertNo(i)
        call pointr(np, xv(1:NDIM,i))
      enddo
!
!  ...evaluate coefficients defining the bilinear map
      do j=1,NDIM
        aux(j,1) = xv(j,2)-xv(j,1)
        aux(j,2) = xv(j,4)-xv(j,1)
        aux(j,3) = xv(j,1)-xv(j,2)+xv(j,3)-xv(j,4)
      enddo
!
!  ...evaluate physical coordinates and their derivatives
      do j=1,NDIM
        Dxdeta(j,1) = aux(j,1) + Eta(2)*aux(j,3)
        Dxdeta(j,2) = aux(j,2) + Eta(1)*aux(j,3)
        X(j) = xv(j,1) + Eta(1)*aux(j,1) + Eta(2)*aux(j,2) &
             + Eta(1)*Eta(2)*aux(j,3)
      enddo
!
!
   end subroutine recta_linear
