c----------------------------------------------------------------------
c
c   routine name       - 
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - 
c
c   purpose            - map j into the set {1,2,...,m} s.t. -1 -> m-1,
c                        0 -> m, 1 -> 1, ..., m -> m, m+1 -> 1, etc.
c
c
c   usage              - call 
c
c   arguments :
c     in:
c     out:
c
c   required  routines - 
c
c----------------------------------------------------------------------
c
      integer function iwrap(j,m)
#include "syscom.blk"
c
      if (j > 0) then
        iwrap = j-((j-1)/m)*m
      else
        iwrap = j+((m-j)/m)*m
      endif
c
      end function iwrap
