      subroutine start_clock(Icount)
#include "syscom.blk"
c
      call system_clock(Icount,ir,im)
c
      end subroutine start_clock
c
c
c
      subroutine stop_clock(Dtime,Icount)
#include "syscom.blk"
c
      call system_clock(ic,ir,im)
c
      if (ic > Icount) then
        Dtime = dble(ic-Icount)/dble(ir)
      elseif (ic < Icount) then
        Dtime = dble(im+ic-Icount)/dble(ir)
      else
        Dtime = 0.d0
      endif
c
      end subroutine stop_clock
c
c
c
      subroutine reset_clock(Dtime,Icount)
#include "syscom.blk"
c
      call stop_clock(Dtime,Icount)
      call start_clock(Icount)
c
      end subroutine reset_clock
