c----------------------------------------------------------------------
c
c     routine name       - shapEbe_0
c
c----------------------------------------------------------------------
c
c     latest revision    - Feb 11
c
c     purpose            - simple wrapper for legendre polynomial
c
c     arguments :
c     in:
c       Xi        - master element coordinates
c       Nord      - order of approximation for the nodes
c     out:
c       NrdofE    - number of element dof
c       ShapE     - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shapEbe_0(T,Nord,Norient,
     .                   NrdofE,ShapE)
c
      use parameters
#include "syscom.blk"
c
c  ...allow for exceeding the order of approximation
      integer, parameter :: ndelta = 2
c
c  ...input output parameters
      dimension ShapE(MAXP)
      dimension poly(0:MAXP+ndelta), dpoly(0:MAXP+ndelta)
c
c  ...compare the 3D edge
      dimension xin(3), norder(19), ne_orient(12), nf_orient(6),
     .  shapnE(3,MAXbrickE), curlnE(3,MAXbrickE)
c
c----------------------------------------------------------------------
c
      iprint=0
c
      ShapE(1:MAXP) = 0.d0
      if ( (T.lt.0.d0).or.(T.gt.1.d0) ) then
        write(*,*) ' shape1E: T OUT OF RANGE ', T
        call pause
      endif
c
      xi = 2.d0*T - 1.d0
      select case(Norient)
      case(0); s =  xi; dsdt =  2.d0
      case(1); s = -xi; dsdt = -2.d0
      end select
c
      NrdofE = Nord
      call legendre_I(s, Nord+ndelta, poly, dpoly)
      do i=1, NrdofE
        ShapE(i) = (dpoly(i)*dsdt)*(-1)**(i)
      enddo
      ShapE(1) = ShapE(1)/2.d0
c
      if (iprint.eq.1) then
        write(*,7001) T,Nord,Norient
 7001   format('shapEbe: T,Nord,Norient = ',f8.3,2i2)
c
        xin = 0.d0; norder = 0; ne_orient = 0; nf_orient = 0;
        nrdofnE = 0; shapnE = 0.d0; curlnE = 0.d0;
c
        xin(1) = T; norder(1) = Nord; ne_orient(1) = Norient;
        call shape3nE(xin,norder,ne_orient,nf_orient,
     .    nrdofnE,shapnE,curlnE)
        do k=1,NrdofE
          write(*,7002) k,ShapE(k), shapnE(1,k),
     .      abs(ShapE(k)-shapnE(1,k))
 7002     format('k = ',i2,' ShapE(k) = ',e12.5,
     .      ' Shape3nE(k) = ',e12.5,
     .      ' Diff = ',e12.5)
        enddo
      endif
      end subroutine
