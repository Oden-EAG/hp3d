c----------------------------------------------------------------------
c
c   routine name       - Dcartesian2polar
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Aug 08
c
c   purpose            - routine switches from cartesian to polar
c                        coordinates; it also returns derivatives of
c                        polar coordinates wrt catesian coordinates
c
c
c   arguments :
c     in:
c              r       - Cartesian coordinates of a point
c     out:
c        polar, Dpolar - polar coordinates of the point, derivatives
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine Dcartesian2polar(r, polar, Dpolar)
c
#include "syscom.blk"
c
      dimension r(2), polar(2), Dpolar(2,2)
c
      data eps /1.d-13/
c
      iprint=0
c
      x = r(1) 
      y = r(2)
c
c  ...determine polar coordinates for the point
      rho = sqrt(x**2+y**2)
      if (rho.lt.eps) then
        polar(2) = 0.d0
      else
        pi = acos(-1.d0)
        if (y.eq.0.d0) then
          if (x.lt.0.d0) then
            polar(2) = pi
          elseif (x.gt.0.d0) then
            polar(2) = 0.d0
          endif
        elseif (y.gt.0.d0) then
          polar(2) = acos(x/rho)
        else
          polar(2) = -acos(x/rho)
        endif
      endif
c      
      polar(1) = rho
c
      Dpolar(1,1) = cos(polar(2))                            !*** compute drho/d(x,y)
      Dpolar(1,2) = sin(polar(2))
c
      if (rho .gt. 0.d0) then                                !*** compute dtheta/d(x,y)
        Dpolar(2,1) = - sin(polar(2))/rho
        Dpolar(2,2) = cos(polar(2))/rho
      else
        write(*,*)'Dcartesian2polar: derivatives of theta
     . component are undefined.'
        stop
      end if
c
      end
