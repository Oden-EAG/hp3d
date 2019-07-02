c----------------------------------------------------------------------
c
c     routine name      - shape_triaH
c
c----------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of triangle H1 bubble
c                         shape functions and their derivatives
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
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             GradH     - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape_triaH(Vshap,Vgrad,N,Nord, NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:2),Vgrad(N,0:2),
     .          ShapH(MAXmdltH),GradH(N,MAXmdltH)
c
c  ...work space for routine shape_edgeH
      dimension shapH1(2:MAXP-1),gradH1(N,2:MAXP-1)
c
c  ...work space for routine IJacobi
      dimension u2(1:MAXP+1,1:MAXP),du2dx(1:MAXP+1,1:MAXP),
     .          du2dt(1:MAXP+1,1:MAXP)
c
c  ...workspace
      dimension gradt(N)
c
c
      iprint=0
c
      if (Nord.lt.3) then
        NrdofH=0
        return
      endif
c
      x = Vshap(1); t = Vshap(0)+Vshap(1)
      call shape_edgeH(Vshap,Vgrad,N,Nord-1, nrdofH1,shapH1,gradH1)
      y = Vshap(2)
      t = Vshap(0)+Vshap(1)+Vshap(2)
      gradt(1:N) = Vgrad(1:N,0)+Vgrad(1:N,1)+Vgrad(1:N,2)
      select case(N)
      case(2)
        if (abs(t-1.d0).gt.1.d-15) then
          write(*,*) 'shape_triaH: t = ',t
          stop 1
        endif
        idec=0
      case(3)
        idec=1
      end select
      call IJacobi(y,t,Nord,Nord,idec, u2,du2dx,du2dt)
c
      k=0
c
c  ...watch for the structure of the loop to enforce hierarchy
c     in polynomial order, np = total order
      do np=3,Nord
c
c  .....order in x
        do i=2,np-1
c
c  .......order in t
          j=np-i
          k=k+1
          ShapH(k) = shapH1(i)*u2(i,j)
          GradH(1:N,k) = gradH1(1:N,i)*u2(i,j)
     .                 + shapH1(i)*du2dx(i,j)*Vgrad(1:N,2)
          if (N.eq.3)
     .    GradH(1:N,k) = GradH(1:N,k)
     .                 + shapH1(i)*du2dt(i,j)*gradt(1:N)
        enddo
      enddo
      NrdofH=k
c
      if (iprint.eq.1) then
        write(*,7001) Vshap(0:2),Nord
 7001   format('shap_triaH: Vshap = ',3f8.3,' Nord = ',i2)
        do k=1,NrdofH
          write(*,7002) k,ShapH(k),GradH(1:N,k)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_triaH
c
c--------------------------------------------------------------------
c
c     routine name      - shape_triaH_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of face H1 bubble
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
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             GradH     - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape_triaH_orient(Vshap,Vgrad,N,Nord,Norient,
     .                              NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:2),Vgrad(N,0:2),
     .          ShapH(MAXmdltH),GradH(N,MAXmdltH)
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
      call shape_triaH(vshap_loc,vgrad_loc,N,Nord, NrdofH,ShapH,GradH)
c
c
      if (iprint.eq.1) then
        write(*,7001) Vshap(0:2),Nord, Norient
 7001   format('shap_triaH_orient: Vshap = ',3f8.3,
     .         ' Nord = ',i2,' Norient = ',i1)
        do k=1,NrdofH
          write(*,7002) k,ShapH(k),GradH(1:N,k)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_triaH_orient


