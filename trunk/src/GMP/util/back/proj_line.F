c----------------------------------------------------------------------
c
c   routine name       - proj_line
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 04
c
c   purpose            - routine projects a vector along a line
c                        given by two points
c
c
c   arguments
c     in:
c                Vect  - vector to project
c                X1,X2 - points defining the line
c     out:
c                Vectl - projected vector
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine proj_line(Vect,X1,X2, Vectl)
c
#include "syscom.blk"
c
      dimension X1(3),X2(3),Vect(3),Vectl(3)
      dimension diff(3)
c
      diff(1:3) = X2(1:3) - X1(1:3)
      rl2 = diff(1)**2 + diff(2)**2 + diff(3)**2
      if (rl2.lt.1.d-20) then
        write(*,7001) rl2
 7001   format('proj_line: rl2 = ',e12.5)
        stop 1
      endif
      aux = Vect(1)*diff(1) + Vect(2)*diff(2) + Vect(3)*diff(3)
      Vectl(1:3) = diff(1:3)*aux/rl2
c
      end


