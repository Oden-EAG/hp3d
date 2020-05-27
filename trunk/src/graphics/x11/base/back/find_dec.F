c----------------------------------------------------------------------
c
c   routine name       - find_dec
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Dec 98
c
c   purpose            - determine value of a refinement flag
c
c   usage              - call find_dec(Xmod, Nflag)
c
c   arguments :
c     in:
c              Xmod    - coordinates of an input point
c     out:
c              Nflag   = 0 for h-refinements
c                      = 1,...,8 for p-refinements
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine find_dec(Xmod, Nflag)
c
      use graphmod
#include "syscom.blk"
c
      dimension Xmod(2)
c
      iprint=0
c
      if (iprint.eq.1) then
        write(*,*) 'find_dec: Xmod = ',Xmod
        write(*,*) 'find_dec: xy_box = '
        do i=0,8
          write(*,7001) ((xy_box(k,j,i),k=1,2),j=1,4)
 7001     format(4(2f6.1,1x))
        enddo
        call pause
      endif
c
      Nflag=-1
      do i=0,8
        if ((Xmod(1).ge.xy_box(1,1,i)).and.
     .      (Xmod(1).le.xy_box(1,2,i)).and.
     .      (Xmod(2).ge.xy_box(2,1,i)).and.
     .      (Xmod(2).le.xy_box(2,4,i)) ) then
           Nflag=i
           exit
        endif
ccc        d = (xy_box(1,i) - Xmod(1))**2 + (xy_box(2,i) - Xmod(2))**2
ccc        d = sqrt(d)
ccc        if (d.lt.dmin) then
ccc          dmin = d
ccc          Nflag=i
ccc        endif
      enddo
      if (iprint.ge.1) then
        write(*,*) 'find_dec: Nflag = ',Nflag
      endif
c
c
      return
      end



