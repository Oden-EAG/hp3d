c-----------------------------------------------------------------------
c
c   routine name       - setzi
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - jan, 1990
c
c   purpose            - routine initiates an integer vector
c
c   usage              - call setzi(Ia,N)
c
c   arguments :
c     in:
c                Ia    - an integer vector
c                N     - number of elements to be zeroed
c
c   required  routines - none
c
c-----------------------------------------------------------------------
c
      subroutine setzi(Ia,N)
c
#include "syscom.blk"
c
      dimension Ia(*)
c
      do 10 i=1,N
        Ia(i)=0
   10 continue
c
      return
      end
