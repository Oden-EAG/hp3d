c----------------------------------------------------------------------
c
c     routine name      - shape_triaV
c
c----------------------------------------------------------------------
c
c     latest revision:  - Apr 14
c
c     purpose:          - routine returns values of triangle H(div)
c                         bubble shape functions and their divergence
c
c     arguments:
c
c     in:
c             Vshap     - triangle vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             Nord      - polynomial order
c
c     out:
c             NrdofV    - number of dof
c             ShapV     - values of the shape functions at the point
c             DivV      - div of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape_triaV(Vshap,Vgrad,Nord, NrdofV,ShapV,DivV)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:2),Vgrad(3,0:2),
     .          ShapV(3,MAXmdltQ),DivV(MAXmdltQ)
c
c  ...work space for routines Legendre and Jacobi
      dimension p1(0:MAXP),p2(0:MAXP+1,0:MAXP)
c
c  ...work space
      dimension vecW(3),cross1(3)
c
c  ...Jacobi weight
      integer alpha
c
      iprint=0
c
c  ...determine face determinant
      call cross_product(Vgrad(1:3,0),Vgrad(1:3,1), cross1)
      call dot_product(cross1,Vgrad(1:3,2), detf)
c
c  ...determine the Whitney shape function
      vecW(1:3) = 0.d0
      do i=0,2
        i1 = mod(i+1,3); i2=mod(i+2,3)
        call cross_product(Vgrad(1:3,i1),Vgrad(1:3,i2), cross1)
        vecW(1:3) = vecW(1:3) + Vshap(i)*cross1(1:3)
      enddo
c
c  ...determine shifted Legendre polynomials in the first variable
      x = Vshap(1); t = Vshap(0) + Vshap(1)
      call LegendreNEW(x,t,Nord-1, p1)
c
c  ...determine shifted Jacobi polynomials in the second variable
      x = Vshap(2); t = Vshap(0) + Vshap(1) + Vshap(2)
      maxalpha = 2*Nord-1
      call Jacobi(x,t,maxalpha,Nord-1, p2)
c
c  ...shape function counter
      k=0
c
c  ...watch for the structure of the loop to enforce hierarchy
c     in polynomial order, np = total order
      do np=0,Nord-1
c
c  .....order in y
        do j=0,np
c
c  .......order in x
          i=np-j
c
c  .......weight in y
c          alpha = 2*i+1
          k=k+1
          ShapV(1:3,k) = p1(i)*p2(i,j)*vecW(1:3)
          DivV(k)      = (np+3)*p1(i)*p2(i,j)*detf
        enddo
c
c  ...end of loop through order of polynomials
      enddo
      NrdofV=k
c
      if (iprint.eq.1) then
        write(*,7001) Vshap(0:2),Nord
 7001   format('shap_triaV: Vshap = ',3f8.3,' Nord = ',i2)
        do k=1,NrdofV
          write(*,7002) k,ShapV(1:3,k),DivV(k)
 7002     format('k = ',i3,' ShapV, DivV = ',3e12.5, 3x, e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_triaV
c
c--------------------------------------------------------------------
c
c     routine name      - shape_triaV_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of face H(div) bubble
c                         shape functions and their divergence
c                         accounting for orientation
c
c     arguments:
c
c     in:
c             Vshap     - triangle vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             Nord      - polynomial order
c             Norient   - face orientation
c
c     out:
c             NrdofV    - number of dof
c             ShapV     - values of the shape functions at the point
c             DivV     - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape_triaV_orient(Vshap,Vgrad,Nord,Norient,
     .                              NrdofV,ShapV,DivV)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:2),Vgrad(3,0:2),
     .          ShapV(3,MAXmdltQ),DivV(MAXmdltQ)
c
c  ...affine coordinates and their gradient reordered according to
c     orientation
      dimension vshap_loc(0:2),vgrad_loc(3,0:2)
c
      iprint=0
c
c  ...rotate affine coordinates acoording to orientation
      select case(Norient)
      case(0)
        n0 = 0; n1 = 1; n2 = 2
      case(1)
        n0 = 1; n1 = 2; n2 = 0
      case(2)
        n0 = 2; n1 = 0; n2 = 1
      case(3)
        n0 = 0; n1 = 2; n2 = 1
      case(4)
        n0 = 1; n1 = 0; n2 = 2
      case(5)
        n0 = 2; n1 = 1; n2 = 0
      end select
      vshap_loc(0) = Vshap(n0); vgrad_loc(1:3,0) = Vgrad(1:3,n0)
      vshap_loc(1) = Vshap(n1); vgrad_loc(1:3,1) = Vgrad(1:3,n1)
      vshap_loc(2) = Vshap(n2); vgrad_loc(1:3,2) = Vgrad(1:3,n2)
c
c  ...determine shape functions in the face coordinates
      call shape_triaV(vshap_loc,vgrad_loc,Nord, NrdofV,ShapV,DivV)
c
      if ((iprint.eq.1).and.(NrdofV.gt.0)) then
        write(*,7001) Vshap(0:2),Nord, Norient
 7001   format('shap_triaV_orient: Vshap = ',3f8.3,
     .         ' Nord = ',i2,' Norient = ',i1)
        do k=1,NrdofV
          write(*,7002) k,ShapV(1:3,k),DivV(k)
 7002     format('k = ',i3,' ShapV, DivV = ',3e12.5,2x,e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_triaV_orient


