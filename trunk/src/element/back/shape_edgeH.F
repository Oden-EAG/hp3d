c-----------------------------------------------------------------------
c
c     routine name      - shape_edgeH
c
c-----------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of edge H1 shape
c                         functions (edge bubbles) and their gradients
c
c     arguments:
c
c     in:
c             Vshap     - edge vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             N         - space dimension (=2,3)
c             Nord      - polynomial order
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             GradH     - the corresponding gradients
c
c-----------------------------------------------------------------------
c
      subroutine shape_edgeH(Vshap,Vgrad,N,Nord, NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:1),Vgrad(N,0:1),
     .          ShapH(MAXP-1),GradH(N,MAXP-1)
c
c  ...work space for routine ILegendre
      dimension rL(2:MAXP),dLdx(2:MAXP),dLdt(2:MAXP)
c
c  ...work space
      dimension gradt(N)
c
      iprint=0
c
      if (Nord.le.1) return
      if ((N.eq.2).and.(N.eq.3)) then
        write(*,*) 'shape_edgeH: N = ',N
        stop 1
      endif
c
      x = Vshap(1); t = Vshap(0)+Vshap(1)
      gradt(1:N) = Vgrad(1:N,0) +  Vgrad(1:N,1)
c
      call ILegendre(x,t,Nord,1, rL,dLdx,dLdt)
      k=0
      do np=2,Nord
        k=k+1
        ShapH(k) = rl(np)
        GradH(1:N,k) = dLdx(np)*Vgrad(1:N,1)
        GradH(1:N,k) = GradH(1:N,k) + dLdt(np)*gradt(1:N)
      enddo
      NrdofH=k
c
c
      if (iprint.eq.1) then
        write(*,7001) Nord
 7001   format('shape_edgeH: Nord = ',i2)
        do i=0,1
          write(*,7002) i,Vshap(i),Vgrad(1:N,i)
 7002     format('i,Vshap,Vgrad = ',i2,2x,e12.5,2x,3e12.5)
        enddo
        do k=1,NrdofH
          write(*,7003) k,ShapH(k),GradH(1:2,k)
 7003     format('k = ',i2,' ShapH,GradH = ',e12.5,2x,3e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape_edgeH
c
c--------------------------------------------------------------------
c
c     routine name      - shape_edgeH_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of edge H1 shape
c                         functions and their derivatives accounting
c                         for orientations
c
c     arguments:
c
c     in:
c             Vshap     - edge vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             N         - space dimension
c             Nord      - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             GradH     - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape_edgeH_orient(Vshap,Vgrad,N,Nord,Norient,
     .                              NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:1),Vgrad(N,0:1),
     .          ShapH(MAXP-1),GradH(N,MAXP-1)
c
c  ...edge affine coordinates
      dimension vshap_loc(0:1), vgrad_loc(1:N,0:1)
c
      iprint=0
      if (Nord.le.1) then
        NrdofH=0; return
      endif
      if (iprint.eq.1) then
        write(*,7001) Vshap(0:1),Nord
        do ivar=1,N
          write(*,7003) ivar,Vgrad(ivar,0:1)
 7003     format('shape_edgeH_orient: Vgrad(',i1,'0:1) = ',2f8.3)
        enddo
      endif
c
c  ...rotate affine coordinates acoording to orientation
      select case(Norient)
      case(0)
        n0 = 0; n1 = 1
      case(1)
        n0 = 1; n1 = 0
      end select
      vshap_loc(0) = Vshap(n0); vgrad_loc(1:N,0) = Vgrad(1:N,n0)
      vshap_loc(1) = Vshap(n1); vgrad_loc(1:N,1) = Vgrad(1:N,n1)
c
      call shape_edgeH(vshap_loc,vgrad_loc,N,Nord, NrdofH,ShapH,GradH)
c
      if (iprint.eq.1) then
        write(*,7001) Vshap(0:1),Nord
 7001   format('shape_edgeH_orient: Vshap,Nord = ',2f8.3,2x,i2)
        do k=1,NrdofH
          write(*,7002) k,ShapH(k),GradH(1:N,k)
 7002       format('k = ',i2,' ShapH,GradH = ',e12.5,2x,3e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_edgeH_orient

