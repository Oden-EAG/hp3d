c----------------------------------------------------------------------
c
c     routine name       - shapEb2_0
c
c----------------------------------------------------------------------
c
c     latest revision    - Feb 11
c
c     purpose            - routine evaluates H(curl) bubbles
c                          for a triangle or a quad
c
c     arguments :
c     in:
c       Type      - face type
c       Xi        - master element coordinates
c       Nord      - order of approximation for 'mdlt'
c       Nort      - orienttaion of 'mdlt'
c     out:
c       NrdofE    - number of element dof
c       ShapE     - values of shape functions
c       CurlE     - values of derivatives of the shape functions
c                   wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapEb2_0(Type,Xi,Nord,Nort, NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
c  ...input output parameters
      character(len=4) :: Type
      dimension Xi(2),vlambda(3),dlambda(2,3),
     .          ShapE(2,MAXmdlqE),CurlE(MAXmdlqE)
c-----------------------------------------------------------------------
c
      select case (Type)
c
      case('mdlt','tria')
c  .....affine coordinate
        vlambda(1)     =  1.d0 - Xi(1) - Xi(2)
        dlambda(1:2,1) = -1.d0;
        do i=1,2
          vlambda(    i+1) = Xi(i)
          dlambda(1:2,i+1) = 0.d0
          dlambda(i,  i+1) = 1.d0
        enddo
        call shapEbt(Xi, Nord,Nort,vlambda,dlambda,NrdofE,ShapE,CurlE)
ccc
ccc   kyungjoo :: will impliment the quad case later...not urgent
      case('mdlq','rect')
        call shapEbq(Xi,Nord,Nort, NrdofE,ShapE,CurlE)
c
      case default
        write(*,*) 'shapEb2: Type = ', Type
        call logic_error(ERR_INVALID_VALUE,
     .    __FILE__,__LINE__)
      end select
c
c
      end subroutine
