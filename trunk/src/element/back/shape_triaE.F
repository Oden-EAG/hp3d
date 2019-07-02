c----------------------------------------------------------------------
c
c     routine name      - shape_triaE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of triangle H(curl)
c                         bubble shape functions and their curls
c
c     arguments:
c
c     in:
c             Vshap     - triangle vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             N         - space dimension
c             Nord      - polynomial order
c
c     out:
c             NrdofE    - number of dof
c             ShapE     - values of the shape functions at the point
c             CurlE     - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape_triaE(Vshap,Vgrad,N,Nord, NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:2),Vgrad(N,0:2),
     .          ShapE(N,MAXmdltE),CurlE(2*N-3,MAXmdltE)
c
c  ...local affine coordinates
      dimension vshapl(0:2),vgradl(N,0:2)
c
c  ...work space for routine shape_edgeE
      dimension shapE1(N,0:MAXP-1),curlE1(2*N-3,0:MAXP-1)
c
c  ...work space for routine IJacobi
      dimension u2(1:MAXP+1,1:MAXP),du2dx(1:MAXP+1,1:MAXP),
     .          du2dt(1:MAXP+1,1:MAXP)
c     P. Gatto, Nov 24, 2014
ccc      dimension u2(1:MAXP+1,1:MAXP-1),du2dx(1:MAXP+1,1:MAXP-1),
ccc     .          du2dt(1:MAXP+1,1:MAXP-1)
c
c  ...work space
      dimension gradt(N),gradu2(N),cross1(2*N-3)
c
      iprint=0
c
      if (Nord.lt.2) then
        NrdofE=0
        return
      endif
c
c  ...dimension of curl
      m = 2*N-3
      t = Vshap(0)+Vshap(1)+Vshap(2)
      select case(N)
      case(2)
        if (abs(t-1.d0).gt.1.d-14) then
          write(*,*) 'shape_triaE: t = ',t
          stop 1
        endif
        idec=0
      case(3)
        idec=1
      end select
      gradt(1:N) = Vgrad(1:N,0)+Vgrad(1:N,1)+Vgrad(1:N,2)
c
c  ...loop through the two families of functions
      do iter=0,1
c
c  .....set up local affine coordinates
        do l=0,2
          l1 = mod(l+iter,3)
          vshapl(l) = Vshap(l1); vgradl(1:N,l) = Vgrad(1:N,l1)
        enddo
c
c  .....determine edge H(curl) shape functions for the first edge
        call shape_edgeE(vshapl,vgradl,N,Nord, nrdofE1,shapE1,curlE1)
c
c  .....use integrated Jacobi polynomials in y starting with p=1
        y = vshapl(2)
        call IJacobi(y,t,Nord-1,Nord,idec, u2,du2dx,du2dt)
c
c  .....shape function counter
        k=iter-1
c
c  .....watch for the structure of the loop to enforce hierarchy
c       in polynomial order, np = total order
        do np=1,Nord-1
c
c  .......order in y
          do j=1,np
c
c  .........order in x
            i=np-j
            k=k+2
            ShapE(1:N,k) = shapE1(1:N,i)*u2(i+1,j)
            gradu2(1:N) = du2dx(i+1,j)*vgradl(1:N,2)
            if (N.eq.3)
     .      gradu2(1:N) = gradu2(1:N) + du2dt(i+1,j)*gradt(1:N)
            call cross(N,shapE1(1:N,i),gradu2(1:N), cross1)
            CurlE(1:m,k) = curlE1(1:m,i)*u2(i+1,j) - cross1(1:m)
          enddo
c
c  .....end of loop through order of polynomials
        enddo
c
c  ...end of loop through families
      enddo
      NrdofE=k
c
c  ...check number of dof
      ndof = Nord*(Nord-1)
      if (ndof.ne.NrdofE) then
        write(*,*) 'shape_triaE: ndof,NrdofE = ',ndof,NrdofE
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7001) Vshap(0:2),Nord
 7001   format('shap_triaE: Vshap = ',3f8.3,' Nord = ',i2)
        do k=1,NrdofE
          write(*,7002) k,ShapE(1:N,k),CurlE(1:m,k)
 7002     format('k = ',i3,' ShapE, CurlE = ',6e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_triaE
c
c--------------------------------------------------------------------
c
c     routine name      - shape_triaE_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of face E1 bubble
c                         shape functions and their derivatives
c                         accounting for orientation
c
c     arguments:
c
c     in:
c             Vshap     - triangle vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             N         - space dimension
c             Nord      - polynomial order
c             Norient   - face orientation
c
c     out:
c             NrdofE    - number of dof
c             ShapE     - values of the shape functions at the point
c             CurlE     - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape_triaE_orient(Vshap,Vgrad,N,Nord,Norient,
     .                              NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:2),Vgrad(N,0:2),
     .          ShapE(N,MAXmdltE),CurlE(2*N-3,MAXmdltE)
c
c  ...affine coordinates and their gradient reordered according to
c     orientation
      dimension vshap_loc(0:2),vgrad_loc(N,0:2)
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
      vshap_loc(0) = Vshap(n0); vgrad_loc(1:N,0) = Vgrad(1:N,n0)
      vshap_loc(1) = Vshap(n1); vgrad_loc(1:N,1) = Vgrad(1:N,n1)
      vshap_loc(2) = Vshap(n2); vgrad_loc(1:N,2) = Vgrad(1:N,n2)
c
c  ...determine shape functions in the face coordinates
      call shape_triaE(vshap_loc,vgrad_loc,N,Nord, NrdofE,ShapE,CurlE)
c
      if ((iprint.eq.1).and.(NrdofE.gt.0)) then
        write(*,7001) Vshap(0:2),Nord, Norient
 7001   format('shap_triaE_orient: Vshap = ',3f8.3,
     .         ' Nord = ',i2,' Norient = ',i1)
        do k=1,NrdofE
          write(*,7002) k,ShapE(1:N,k),CurlE(1:2*N-3,k)
 7002     format('k = ',i3,' ShapE, CurlE = ',6e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_triaE_orient


