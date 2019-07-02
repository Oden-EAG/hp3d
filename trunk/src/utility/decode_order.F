c-----------------------------------------------------------------------
c
c     routine name      - decode_orderb
c
c-----------------------------------------------------------------------
c
c     latest revision:  - Sep 05
c
c     purpose:          - routine decodes order of approximation for
c                         a brick middle node
c
c     arguments:
c
c     in:
c             Norderb   - the encoded order
c     out:
c             Nord1,2,3 - orders in directions 1,2,3
c
c-----------------------------------------------------------------------
c
      subroutine decode_orderb(Norderb, Nord1,Nord2,Nord3)
c
      use parameters , only : MAXP      
c      
#include "syscom.blk"
c
c     RECALL: Norderb = (Nord1*(MAXP+1)+Nord2)*(MAXP+1)+Nord3
c
      naux  = Norderb/(MAXP+1)
      Nord3 = Norderb-naux*(MAXP+1)
      Nord1 = naux/(MAXP+1)
      Nord2 = naux-Nord1*(MAXP+1)
c
      if (Nord1<1.or.Nord1>MAXP) then
        write(*,*)'decode_orderb: OUT OF RANGE Nord1 = ',Nord1
        stop
      endif
      if (Nord2<1.or.Nord2>MAXP) then
        write(*,*)'decode_orderb: OUT OF RANGE Nord2 = ',Nord2
        stop
      endif
      if (Nord3<1.or.Nord3>MAXP) then
        write(*,*)'decode_orderb: OUT OF RANGE Nord3 = ',Nord3
        stop
      endif
c
      end subroutine decode_orderb
c
c-----------------------------------------------------------------------
c
c     routine name      - decode_orderq
c
c-----------------------------------------------------------------------
c
c     latest revision:  - Sep 05
c
c     purpose:          - routine decodes order of approximation for
c                         a quad middle node
c
c     arguments:
c
c     in:
c             Norderq   - the encoded order
c     out:
c             Nord1,2   - orders in directions 1,2
c
c-----------------------------------------------------------------------
c
      subroutine decode_orderq(Norderq, Nord1,Nord2)
c
      use parameters , only : MAXP      
c      
#include "syscom.blk"
c
c     RECALL: Norderq = Nord1*(MAXP+1)+Nord2
c
      Nord1 = Norderq/(MAXP+1)
      Nord2 = Norderq-Nord1*(MAXP+1)
c
      if (Nord1<1.or.Nord1>MAXP) then
        write(*,*)'decode_orderq: OUT OF RANGE Nord1 = ',Nord1
        stop
      endif
      if (Nord2<1.or.Nord2>MAXP) then
        write(*,*)'decode_orderq: OUT OF RANGE Nord2 = ',Nord2
        stop
      endif
c
      end subroutine decode_orderq
