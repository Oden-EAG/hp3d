      subroutine arctan(Xc,Yc,Ang)
c
c----------------------------------------------------------------------
c
c   routine name       - arctan
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - jan, 1992
c
c   purpose            - routine computes the angle whose tangent
c                        is Yc/Xc
c
c   usage              - call arctan(Xc,Yc,Ang)
c
c   arguments
c             in       - Xc,Yc - usually the cartesian coordinates
c                                of given point
c             out      - Ang - the angle between radius of given
c                              point and x-axis ( 0 < Ang < 2*Pi )
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
#include "syscom.blk"
c
      data small,pi /1.d-10,3.14159265358979/
c
      if (dabs(Xc).le.small) then
         if (dabs(Yc).le.small) then
            Ang=0.0d0
         else
            if (Yc.ge.0.0d0) then
               Ang=0.5d0*pi
            else
               Ang = 1.5d0*pi
            endif
         endif
      else
         if (Xc.le.0.0d0) then
            Ang=atan(Yc/Xc)+pi
         else
            Ang=atan(Yc/Xc)
         endif
      endif
c
      return
      end

