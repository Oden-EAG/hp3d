c----------------------------------------------------------------------
c
c     routine name      - shape_tetrE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of tetrahedron H(curl)
c                         bubble shape functions and their curls
c
c     arguments:
c
c     in:
c             Vshap     - tetrahedron vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             Nord      - polynomial order
c
c     out:
c             NrdofE    - number of dof
c             ShapE     - values of the shape functions at the point
c             CurlE     - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape_tetrE(Vshap,Vgrad,Nord, NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:3),Vgrad(3,0:3),
     .          ShapE(3,MAXmdlnE),CurlE(3,MAXmdlnE)
c
c  ...local affine coordinates
      dimension vshapl(0:3),vgradl(3,0:3)
c
c  ...work space for routine shape_edgeE
      dimension shapE1(3,0:MAXP-1),curlE1(3,0:MAXP-1)
c
c  ...work space for routine IJacobi
      dimension u2(1:MAXP+1,1:MAXP-1),
     .          du2dx(1:MAXP+1,1:MAXP),du2dt(1:MAXP+1,1:MAXP),
     .          u3(1:MAXP+1,1:MAXP),du3dx(1:MAXP+1,1:MAXP)
c
c  ...work space
      dimension crossE2(3),crossE3(3),tmpgrad(3)
c
      double precision z
c
      iprint=0
c
      if (Nord.le.2) then
        NrdofE=0
        return
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
c  .....determine H(curl) shape functions for the first edge
        call shape_edgeE(vshapl,vgradl,3,Nord, nrdofE1,shapE1,curlE1)
c
c  .....determine integrated shifted scaled Jacobi polynomials in y
        y = vshapl(2); z = vshapl(3)
        call IJacobi(y,1.d0-z,Nord-2,Nord-2,1, u2,du2dx,du2dt)
c
c  .....determine integrated shifted Jacobi polynomials in z
        call IJacobi(z,1.d0,Nord-2,Nord-2,0, u3,du3dx,void)
c
        k=iter-2
c
c  .....watch for the structure of the loop to enforce hierarchy
c       in polynomial order, np = total order
        do np=2,Nord-1
c
cccccc  FF Oct 14: Consistent ordering (like in H1)
          do nxy=1,np-1
c
            do i=0,nxy-1
c
              j = nxy-i
c
              l = np-nxy
c  .......order in z
c          do l=1,np-1
c
c  .........combined order in x and y
c            nxy=np-l
c
c  .........order in y
c            do j=1,nxy
c
c  ...........order in x
c             i=nxy-j
c
              k=k+3
              ShapE(1:3,k) = shapE1(1:3,i)*u2(i+1,j)*u3(nxy,l)
ccccc  .......FF Oct 14: We had to use du2dt in some way.
              tmpgrad(1:3) = du2dx(i+1,j)*vgradl(1:3,2)
     .                     - du2dt(i+1,j)*vgradl(1:3,3)
              call cross(3,shapE1(1:3,i),tmpgrad(1:3), crossE2)
              call cross(3,shapE1(1:3,i),vgradl(1:3,3), crossE3)
              CurlE(1:3,k) = curlE1(1:3,i)*u2(i+1,j)*u3(nxy,l)
     .                     - crossE2(1:3)*u3(nxy,l)
     .                     - u2(i+1,j)*du3dx(nxy,l)*crossE3(1:3)
            enddo
          enddo
c
c  .....end of loop through total polynomial order
        enddo
c
c  ...end of loop through families
      enddo
      NrdofE=k
c
c  ...check number of dof
      ndof = (Nord-2)*(Nord-1)*Nord/2
      if (NrdofE.ne.ndof) then
        write(*,*) 'shape_tetrE: NrdofE,ndof = ',NrdofE,ndof
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7001) Vshap(0:3),Nord
 7001   format('shape_tetrE: Vshap = ',4f8.3,' Nord = ',i2)
        do k=1,NrdofE
          write(*,7002) k,ShapE(1:3,k),CurlE(1:3,k)
 7002     format('k = ',i3,' ShapE, CurlE = ',3e12.5,2x,3e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_tetrE
