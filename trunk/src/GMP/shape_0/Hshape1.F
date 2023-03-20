c----------------------------------------------------------------------
c
c   routine name       - Hshape1
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine evaluates 1D hierarchical C1-conforming
c                        shape functions, the first four are the standard
c                        Hermite polynomials, the next ones are
c                        double integrals of Legendre polynomials
c
c  arguments
c               Nord   - order of approximation
c               Xi     - master element coordinate
c     out:
c               Val    - values of the shape functions
c               Dval   - the corresponding first derivatives
c               Ddval  - the corresponding second derivatives
c
c----------------------------------------------------------------------
c
      subroutine Hshape1(Nord,Xi, Val,Dval,Ddval)
c
      use parameters , only : MAXP
      implicit none
c
      integer :: Nord
      real(8) :: Xi
      real(8) :: Val(MAXP+1), Dval(MAXP+1), Ddval(MAXP+1)
c
      integer :: j
c
c  ...check data
      if ((Nord.lt.3).or.(Nord.gt.MAXP)) then
        write(*,7001) Nord
 7001   format(' Hshape1: Nord = ',i3)
        stop 1
      endif
      if ((Xi.lt.0.d0).or.(Xi.gt.1.d0)) then
        write(*,7002) Xi
 7002   format(' Hshape1: WARNING! Xi = ',f8.3)
        call pause
      endif
c
      Val(1)   = 2.d0*Xi**3 - 3.d0*Xi**2 + 1.d0
      Dval(1)  = 6.d0*Xi**2 - 6.d0*Xi
      Ddval(1) = 12.d0*Xi - 6.d0
c
      Val(2)   = -2.d0*Xi**3 + 3.d0*Xi**2
      Dval(2)  = -6.d0*Xi**2 + 6.d0*Xi
      Ddval(2) = -12.d0*Xi + 6.d0
c
      Val(3)   =  Xi**3 - 2.d0*Xi**2 + Xi
      Dval(3)  =  3.d0*Xi**2 - 4.d0*Xi + 1.d0
      Ddval(3) =  6.d0*Xi - 4.d0
c
      Val(4)   = Xi**3 - Xi**2
      Dval(4)  = 3.d0*Xi**2 - 2.d0*Xi
      Ddval(4) = 6.d0*Xi - 2.d0
      if (Nord.eq.3) return
c
      Val(5)   = Xi**4 - 2.d0*Xi**3 + Xi**2
      Dval(5)  = 4.d0*Xi**3 - 6.d0*Xi**2 + 2.d0*Xi
      Ddval(5) = 12.d0*Xi**2 - 12.d0*Xi + 2.d0
c
      do j=6,Nord+1
        Val(j)   =   Val(j-1)*(2.d0*Xi-1.d0)
        Dval(j)  =  Dval(j-1)*(2.d0*Xi-1.d0) +  Val(j-1)*2.d0
        Ddval(j) = Ddval(j-1)*(2.d0*Xi-1.d0) + Dval(j-1)*4.d0
      enddo
c
c
      end subroutine Hshape1
