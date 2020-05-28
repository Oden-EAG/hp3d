c----------------------------------------------------------------------
c
c   routine name       - leg_coeffs
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Mar 02 (Pawel Solin)
c
c   purpose            - routine initiates constants for normalized
c                        integrated legendre polynomials divided
c                        by (1 - x*x)
c
c   arguments          - none
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine leg_coeffs
c
#include "syscom.blk"
c
#include "cint.blk"
c
c  ...to define 1D integrated Legendre polynomials divided
c     by (1 - x*x) -> this is a base for 1D shape functions
      common /legend0/ xleg0(9,9)
c
ccc      write(*,*) 'USING PRIMITIVES OF LEGENDRE POLYNOMIALS.'
c
c  ...to generate the coefficients in maple
c     > divide(int(orthopoly[P](N, y), y=-1..x), 1 - x^2, 'fN');
c
c
      xleg0(1:9,1:9) = 0.d0
c
      xleg0(1,1)=   -1.d0/2.d0
c
      xleg0(2,2)=   -1.d0/2.d0
c
      xleg0(3,1)=    1.d0/8.d0
      xleg0(3,3)=   -5.d0/8.d0
c
      xleg0(4,2)=    3.d0/8.d0
      xleg0(4,4)=   -7.d0/8.d0
c
      xleg0(5,1)=   -1.d0/16.d0
      xleg0(5,3)=    7.d0/8.d0
      xleg0(5,5)=   -21.d0/16.d0
c
      xleg0(6,2)=   -5.d0/16.d0
      xleg0(6,4)=    15.d0/8.d0
      xleg0(6,6)=   -33.d0/16.d0
c
      xleg0(7,1)=    5.d0/128.d0
      xleg0(7,3)=   -135.d0/128.d0
      xleg0(7,5)=    495.d0/128.d0
      xleg0(7,7)=   -429.d0/128.d0
c
      xleg0(8,2)=    35.d0/128.d0
      xleg0(8,4)=   -385.d0/128.d0
      xleg0(8,6)=    1001.d0/128.d0
      xleg0(8,8)=   -715.d0/128.d0
c
      xleg0(9,1)=   -7.d0/256.d0
      xleg0(9,3)=    77.d0/64.d0
      xleg0(9,5)=   -1001.d0/128.d0
      xleg0(9,7)=    1001.d0/64.d0
      xleg0(9,9)=   -2431.d0/256.d0
c
      return
      end
