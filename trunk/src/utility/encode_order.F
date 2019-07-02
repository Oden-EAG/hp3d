c-----------------------------------------------------------------------
c
c     routine name      - encode_orderb
c
c-----------------------------------------------------------------------
c
c     latest revision:  - Sep 05
c
c     purpose:          - routine encodes order of approximation for
c                         a brick middle node
c
c     arguments:
c
c     in:
c             Nord1,2,3 - orders in directions 1,2,3
c     out:
c             Norderb   - the encoded order
c
c-----------------------------------------------------------------------
c
      subroutine encode_orderb(Nord1,Nord2,Nord3, Norderb)
c
      use parameters , only : MAXP
c      
#include "syscom.blk"
c
      if (Nord1<1.or.Nord1>MAXP) then
        write(*,*)'encode_orderb: OUT OF RANGE Nord1 = ',Nord1
        stop
      endif
      if (Nord2<1.or.Nord2>MAXP) then
        write(*,*)'encode_orderb: OUT OF RANGE Nord2 = ',Nord2
        stop
      endif
      if (Nord3<1.or.Nord3>MAXP) then
        write(*,*)'encode_orderb: OUT OF RANGE Nord3 = ',Nord3
        stop
      endif
c
      Norderb = (Nord1*(MAXP+1)+Nord2)*(MAXP+1)+Nord3
c
      end subroutine encode_orderb
c
c-----------------------------------------------------------------------
c
c     routine name      - encode_orderq
c
c-----------------------------------------------------------------------
c
c     latest revision:  - Sep 05
c
c     purpose:          - routine encodes order of approximation for
c                         a quad middle node
c
c     arguments:
c
c     in:
c             Nord1,2   - orders in directions 1,2
c     out:
c             Norderq   - the encoded order
c
c-----------------------------------------------------------------------
c
      subroutine encode_orderq(Nord1,Nord2, Norderq)
c
      use parameters , only : MAXP
c      
#include "syscom.blk"
c
      if (Nord1<1.or.Nord1>MAXP) then
        write(*,*)'encode_orderq: OUT OF RANGE Nord1 = ',Nord1
        stop
      endif
      if (Nord2<1.or.Nord2>MAXP) then
        write(*,*)'encode_orderq: OUT OF RANGE Nord2 = ',Nord2
        stop
      endif
c
      Norderq = Nord1*(MAXP+1)+Nord2
c
      end subroutine encode_orderq
