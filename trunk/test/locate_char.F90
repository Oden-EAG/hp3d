!> @date Mar 2023
program test_locate_char
!
   implicit none
!
   integer :: NPASS
!
   character(4) :: a(10)
   integer      :: loc
!
   NPASS = 1
!
   a = 'none'; a(7) = 'curl'
!
   call locate_char('curl',a,10, loc)
   if (loc .ne. 7) NPASS = 0
!
   call locate_char('cont',a,10, loc)
   if (loc .ne. 0) NPASS = 0
!
!
   if (NPASS.ne.1) stop 1
!
end program test_locate_char
