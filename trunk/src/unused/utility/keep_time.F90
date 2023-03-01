   subroutine start_clock(Icount)
      implicit none
      integer :: Icount
      integer :: ir,im
!
      call system_clock(Icount,ir,im)
!
   end subroutine start_clock
!
!
!
   subroutine stop_clock(Dtime,Icount)
      implicit none
      real(8) :: Dtime
      integer :: Icount
      integer :: ic,ir,im
!
      call system_clock(ic,ir,im)
!
      if (ic > Icount) then
        Dtime = dble(ic-Icount)/dble(ir)
      elseif (ic < Icount) then
        Dtime = dble(im+ic-Icount)/dble(ir)
      else
        Dtime = 0.d0
      endif
!
   end subroutine stop_clock
!
!
!
   subroutine reset_clock(Dtime,Icount)
      implicit none
      real(8) :: Dtime
      integer :: Icount
!
      call stop_clock(Dtime,Icount)
      call start_clock(Icount)
!
   end subroutine reset_clock
