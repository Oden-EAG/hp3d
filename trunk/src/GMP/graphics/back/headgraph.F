c----------------------------------------------------------------------c----------------------------------------------------------------------
c
c   routine name       - headgraph
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - June,2002
c
c   purpose            - display the head model
c
c   arguments :
c
c   required  routines - 
c
c----------------------------------------------------------------------
c
      subroutine headgraph
c
c
#include "syscom.blk"
c
 10   continue
c      
      write(*,*)
      write(*,*)'headgraph: SELECT:'
      write(*,*)'           RETURN.............................0'
      write(*,*)'           OUTSIDE THE HEAD...................1'
      write(*,*)'           INDE THE HEAD......................2'
      write(*,*)'           INSIDE AND OUTSIDE THE HEAD........3'
      read(*,*) idec
c
c     
      if(idec.eq.0)then         
         return
      elseif(idec.eq.1)then
         call input_head(1)
         call graphg
c
      elseif(idec.eq.2)then
         call input_head(2)
         call graphg
c
      elseif(idec.eq.3)then
         call input_head(3)
         call graphg
c
      else
         write(*,*)'WRONG SELECTION',idec
         goto 10
      endif
c
      end
      
