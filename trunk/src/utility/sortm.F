      subroutine sortm(Nrvistr,Igtr,Rgtrz)
c
c----------------------------------------------------------------------
c
c   routine name       - sort
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - jan, 1992
c
c   purpose            - routine sorts all visible triangles in order
c                        back-to-front according to z-coordinate
c                        of their mid-points (array RGTRZ) and stores
c                        sorted numbers in the array IGTR
c
c   usage              - call sortz
c
c   arguments          - none
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
c
#include "syscom.blk"


      dimension Igtr(*),Rgtrz(*)


c
c  ...initialize index array
      do 10 i=1,NRVISTR
        IGTR(i) = i
   10 continue
      if (NRVISTR.eq.1) return

c
c  ...start sorting (heapsort algorithm)
      l=NRVISTR/2+1
      ir=NRVISTR
   20 continue
      if (l.gt.1) then
        l=l-1
        index=IGTR(l)
        q=RGTRZ(index)
      else
        index=IGTR(ir)
        q=RGTRZ(index)
        IGTR(ir)=IGTR(1)
        ir = ir - 1
        if (ir.eq.1) then
          IGTR(1) = index

ccccc     return
          go to 700

        endif
      endif
      i=l
      j=l*2
   30 continue
      if (j.gt.ir) go to 40



cwr05.16.00
ccccc if ((j.lt.ir).and.(RGTRZ(IGTR(j)).lt.RGTRZ(IGTR(j+1)))) j=j+1
      if(j.lt.ir) then
        if(RGTRZ(IGTR(j)).lt.RGTRZ(IGTR(j+1))) j=j+1
      endif



      if (q.lt.RGTRZ(IGTR(j))) then
        IGTR(i) = IGTR(j)
        i=j
        j=j+j
      else
        j=ir+1
      endif
      go to 30
   40 continue
      IGTR(i) = index
      go to 20
c

  700 continue
c      write(*,*)'SORT: IGTR=',(igtr(i),i=1,nrvistr)
c      pause

      end








