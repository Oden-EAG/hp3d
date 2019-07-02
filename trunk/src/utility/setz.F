c-----------------------------------------------------------------------c
c   routine name       - setz
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - jan 1990
c
c   purpose            - initiates a real array with zeros
c
c   usage              - call setz(a,n)
c
c   arguments :
c     in:
c                 a    - a real array
c                 n    - total length of the array      
c
c   required  routines - none
c
c-----------------------------------------------------------------------c
      subroutine setz(a,n)
c      
#include "syscom.blk"
c      
      dimension a(*)
c      
      do 10 i=1,n
        a(i)=0.d0
   10 continue
c      
      return
      end
