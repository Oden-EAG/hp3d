c----------------------------------------------------------------------
c
c   routine name       - build_u
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 10
c
c   purpose            - evaluates block u for contruction of H(curl)
c                        inner modes
c   arguments :
c     in:
c                ShapH - Barycentric Coordinates
c               DshapH - and their derivatives wrt master coor
c                 Nord - order of polynomials needed
c
c     out:
c                   Ub - evaluation of building block functions u
c                  Dub - evaluation of building block gradients;
c                        gradients wrt master coordinates
c
c----------------------------------------------------------------------
c
      subroutine build_u(ShapH,DshapH,Nord, Ub,Dub)
c
      use parameters
#include "syscom.blk"
c
c  ...exceed maximum order of approximation
      integer, parameter :: ndelta = 1
c
c  ...input output parameters
      dimension ShapH(2),DshapH(3,2)
      dimension Ub(1:MAXP+ndelta),Dub(3,1:MAXP+ndelta)
c
c  ...local variables
      dimension dx(3),ds(1:3)
      dimension poly(1:MAXP+ndelta),dpoly(1:MAXP+ndelta)
c
c----------------------------------------------------------------------
c
c ....test prints
      iprint = 0
c
      if (Nord.gt.MAXP+ndelta) then
        write(*,*)'build_u: Nord = ',Nord
        stop
      endif
c
c  ...temporary trick not to divide by zero
      if ( (ShapH(1).eq.0.d0).and.(ShapH(2).eq.0.d0) ) then
        ShapH(1) = 1e-15; ShapH(2) = 1e-15;
      endif
c
c  ...compute local coordinate
      x = (ShapH(2) - ShapH(1))/(ShapH(2) + ShapH(1))
      dx(1:3) = 2.d0*(DshapH(1:3,2)*ShapH(1) - DshapH(1:3,1)*ShapH(2))/
     .          (ShapH(2) + ShapH(1))**2
c
c  ...compute Weighted Integrated Jacobi polynomials
c     for n = 1:Nord
      a = 0.d0
      call build_Jacobi_W_I(a,x,Nord, poly,dpoly)
c
c  ...compute block Ub
      do i=1,Nord
        s = (ShapH(2) + ShapH(1))**i
        ds(1:3) = (ShapH(2) + ShapH(1))**(i-1)
     .           *(DshapH(1:3,2) + DshapH(1:3,1))*i
        Ub(i) = poly(i)*s
        Dub(1:3,i) = dpoly(i)*dx(1:3)*s + poly(i)*ds(1:3)
      enddo
c
c
      end subroutine
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - build_v
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 10
c
c   purpose            - evaluates block v for contruction of H(curl)
c                        inner modes
c   arguments :
c     in:
c                ShapH - Barycentric Coordinates
c               DshapH - and their derivatives wrt master coor
c                 Nord - order of polynomials needed
c
c     out:
c                   Vb - evaluation of building block functions v
c                  Dvb - evaluation of building block gradients;
c                        gradients wrt master coordinates
c
c----------------------------------------------------------------------
c
      subroutine build_v(ShapH,DshapH,Nord, Vb,Dvb)
c
      use parameters
#include "syscom.blk"
c
c  ...exceed maximum order of approximation
      integer, parameter :: ndelta = 1
c
c  ...input output parameters
      dimension ShapH(2),DshapH(3,2)
      dimension Vb(1:MAXP+ndelta,1:MAXP+ndelta),
     .          Dvb(3,1:MAXP+ndelta,1:MAXP+ndelta)
c
c  ...local variables
      dimension dx(3)
      dimension poly(1:MAXP+ndelta),dpoly(1:MAXP+ndelta)
c
c----------------------------------------------------------------------
c
c ....test prints
      iprint = 0
c
      if (Nord.gt.MAXP+ndelta) then
        write(*,*)'build_v: Nord = ',Nord
        stop
      endif
c
c  ...temporary trick not to divide by zero
      if ( ShapH(2).eq.1.d0 ) then
        ShapH(2) = 1.d0 - 1e-14
      endif

c
c  ...compute local coordinate
      x = (2.d0*ShapH(1) + ShapH(2) - 1.d0)/(1.d0 - ShapH(2))
      dx(1:3) = 2.d0*(DshapH(1:3,2)*ShapH(1) - DshapH(1:3,1)*ShapH(2) +
     .                DshapH(1:3,1))/(1.d0 - ShapH(2))**2
c
      do i=1,Nord
c
c  .....compute Weighted Integrated Jacobi polynomials
        a = 2*i-1
        call build_Jacobi_W_I(a,x,Nord, poly,dpoly)
c
c  .....compute block Vb
        do j=1,Nord
          Vb(i,j) = poly(j)*(1.d0 - ShapH(2))**j
          Dvb(1:3,i,j) = dpoly(j)*dx(1:3)*(1.d0 - ShapH(2))**j +
     .               poly(j)*j*(1.d0 - ShapH(2))**(j-1)*
     .               ( - DshapH(1:3,2))
        enddo
      enddo
c
c
      end subroutine
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - build_w
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 10
c
c   purpose            - evaluates block w for contruction of H(curl)
c                        inner modes
c   arguments :
c     in:
c                ShapH - Barycentric Coordinate
c               DshapH - and derivative wrt master coor
c                 Nord - order of polynomials needed
c
c     out:
c                   Wb - evaluation of building block functions w
c                  Dwb - evaluation of building block gradients;
c                        gradients wrt master coordinates
c
c----------------------------------------------------------------------
c
      subroutine build_w(ShapH,DshapH,Nord, Wb,Dwb)
c
      use parameters
#include "syscom.blk"
c
c     ...exceed maximum order of approximation
      integer, parameter :: ndelta = 1

c
c  ...input output parameters
      dimension DshapH(3)
      dimension Wb(1:MAXP+ndelta,1:MAXP+ndelta,1:MAXP+ndelta),
     .          Dwb(3,1:MAXP+ndelta,1:MAXP+ndelta,1:MAXP+ndelta)
c
c  ...local variables
      dimension dx(3)
      dimension poly(1:MAXP+ndelta),dpoly(1:MAXP+ndelta)
c
c----------------------------------------------------------------------
c
c ....test prints
      iprint = 0
c
      if (Nord.gt.MAXP+ndelta) then
        write(*,*)'build_w: Nord = ',Nord
        stop
      endif
c
c  ...compute local coordinate
      x = 2.d0*ShapH - 1.d0
      dx(1:3) = 2.d0*DshapH(1:3)
c
      do i=1,Nord
        do j=1,Nord
c
c  .......compute weighted integrated Jacobi polynomials
          a = 2*i+2*j-2
          call build_Jacobi_W_I(a,x,Nord, poly,dpoly)
c
c  .......compute block Wb
          do k=1,Nord
            Wb(i,j,k) = poly(k)
            Dwb(1:3,i,j,k) = dpoly(k)*dx(1:3)
          enddo
        enddo
      enddo
c
c
      end subroutine
