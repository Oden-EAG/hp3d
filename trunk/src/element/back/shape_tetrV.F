c----------------------------------------------------------------------
c
c     routine name      - shape_tetrV
c
c----------------------------------------------------------------------
c
c     latest revision:  - Apr 14
c
c     purpose:          - routine returns values of tetrahedron H(div)
c                         bubble shape functions and their divergence
c
c     arguments:
c
c     in:
c             Vshap     - tetrngle vertex shape functions (affine
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
      subroutine shape_tetrV(Vshap,Vgrad,Nord, NrdofV,ShapV,DivV)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:3),Vgrad(3,0:3),
     .          ShapV(3,MAXmdlnV),DivV(MAXmdlnV)
c
c  ...local affine coordinates
      dimension vshapl(0:3),vgradl(3,0:3)
c
c  ...work space for routine Legendre and Jacobi
      dimension p1(0:MAXP),p2(0:MAXP+1,0:MAXP)
c
c  ...work space for routine IJacobi
      dimension u3(1:MAXP+1,1:MAXP),du3dx(1:MAXP+1,1:MAXP)
c
c  ...work space
      dimension vecW(3),cross1(3)
c
c  ...Jacobi weight
      integer alpha
c
      double precision z
c
      iprint=0
      if (Nord.lt.2) then
        NrdofV=0; return
      endif
c
c  ...loop through the three families of functions
      do iter=0,2
c
c  .....set up local affine coordinates
        do l=0,3
          l1 = mod(l+iter,4)
          vshapl(l) = Vshap(l1); vgradl(1:3,l) = Vgrad(1:3,l1)
        enddo
c
c  .....determine the Whitney shape function
        vecW(1:3) = 0.d0
        do i=0,2
          i1 = mod(i+1,3); i2=mod(i+2,3)
          call cross_product(vgradl(1:3,i1),vgradl(1:3,i2), cross1)
          vecW(1:3) = vecW(1:3) + vshapl(i)*cross1(1:3)
        enddo
c
c  .....determine determinant
        call cross_product(vgradl(1:3,0),vgradl(1:3,1), cross1)
        call dot_product(cross1,vgradl(1:3,2), detf)
c
c  .....determine dot product of Whitney shape function and grad lambda_3
        call dot_product(vecW, vgradl(1:3,3), dotp)
c
c  .....determine shifted scaled Legendre polynomials in lambda_1
        x = vshapl(1); t = vshapl(0) + vshapl(1)
        call LegendreNEW(x,t,Nord-1, p1)
c
c  .....determine shifted Jacobi polynomials in the second variable
        y = vshapl(2); t = vshapl(0) + vshapl(1) + vshapl(2)
        maxalpha = 2*Nord-1
        call Jacobi(y,t,maxalpha,Nord-1, p2)
c
c  .....determine shifted integrated polynomials in lambda_3
        z = vshapl(3); t = 1.d0; imax = Nord-1
        call IJacobi(z,t,imax,Nord,0, u3,du3dx,void)
c
c  .....shape function counter
        k=iter-2
c
c  .....watch for the structure of the loop to enforce hierarchy
c       in polynomial order, np = total order
        do np=0,Nord-1
c
c  .......order in lambda_3
          do l=1,np
            ij = np-l+1 ! Jacobi weight for the integrated Jacobi polynomial
c
c  .........order in lambda_2
            do j=0,np-l
c
c  ...........order in lambda_1
              i=np-l-j
c
c  ...........weight in lambda2
c              alpha = 2*i+1
              k=k+3
              ShapV(1:3,k) = p1(i)*p2(i,j)*u3(ij,l)*vecW(1:3)
              DivV(k)      = (np-l+3)*p1(i)*p2(i,j)*detf*u3(ij,l)
     .                     + p1(i)*p2(i,j)*du3dx(ij,l)*dotp
            enddo
          enddo
c
c  .....end of loop through total order of polynomials
        enddo
c
c  ...end of loop through families
      enddo
      NrdofV=k
c
c  ...check number of dof
      ndof = (Nord-1)*Nord*(Nord+1)/2
      if (ndof.ne.NrdofV) then
        write(*,*) 'shape_tetrV: ndof,NrdofV = ',ndof,NrdofV
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7001) Vshap(1:3),Nord
 7001   format('shap_tetrV: Xi = ',3f8.3,' Nord = ',i2)
        do k=1,NrdofV
          write(*,7002) k,ShapV(1:N,k),DivV(k)
 7002     format('k = ',i3,' ShapV, DivV = ',3e12.5,2x,e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_tetrV
