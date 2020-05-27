c----------------------------------------------------------------------
c
c   routine name       - Hbasis
c
c----------------------------------------------------------------------
c
c   latest revision    - August 03
c
c   purpose            - routine generate ratio surface
c
c  arguments
c               Id     - the number of the basis
c               Xi     - a value of the parameter
c     out:
c               Val    - value of the basis function
c               Dval   - Derivatives of the basis function
c
c----------------------------------------------------------------------
c
      subroutine Hbasis(Id,Xi,Val,Dval)
c
#include "syscom.blk"
c
c
      select case(Id)
      case(1)
         Val = 2*Xi**3 - 3*Xi**2 +1
         Dval = 6*Xi**2 - 6*Xi
      case(2)
         Val = -2*Xi**3 + 3*Xi**2
         Dval = -6*Xi**2 + 6*Xi
      case(3)
         Val =  Xi**3 - 2*Xi**2 + Xi
         Dval =  3*Xi**2 - 4*Xi + 1
      case(4)
         Val = Xi**3 - Xi**2
         Dval = 3*Xi**2 - 2*Xi
      case(5)
         Val = Xi**2*(1 - Xi)**2
         Dval = Xi*2*(1 - Xi)**2 - Xi**2*(1 - Xi)*2
      end select
c
c
      end
