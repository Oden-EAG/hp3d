c----------------------------------------------------------------------
c
c   routine name       - collect_tangent
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 07
c
c   purpose            - routine collects tanget vector for a curve
c                        at one of its endpoints and its derivative
c                        wrt to the curve parameter
c
c   arguments
c     in:
c              Nc      - a Hermite curve
c              Np      - a point
c     out:
c              Tangent - the tangent vector
c
c----------------------------------------------------------------------
c
      subroutine collect_tangent(Nc,Np, Tangent,Dtangent)
c
      use GMP
c
#include "syscom.blk"
c
      dimension Tangent(1:3), Dtangent(1:3)
c
c  ...locate the point
      call locate(Np,CURVES(Nc)%EndPoNo(1:2),2, i)
      if (i.eq.0) then
        write(*,7001) Nc,Np
 7001   format('collect_tangent: Nc,Np = ',2i5)
        stop 1
      endif
c
      select case(i)
      case(1)
        Tangent(1:3)  = CURVES(Nc)%Rdata(1:3)
        Dtangent(1:3) = CURVES(Nc)%Rdata(7:9)
      case(2)
        Tangent(1:3) = -CURVES(Nc)%Rdata(4:6)
c
c  .....second order derivatives DO NOT change the sign
c       with orientation !!!
        Dtangent(1:3) = CURVES(Nc)%Rdata(10:12)
      end select
c
c
      end
