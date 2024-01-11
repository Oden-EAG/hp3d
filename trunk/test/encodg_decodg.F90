!> @date Mar 2023
program test_encodg_decodg
!
   implicit none
!
   integer :: NPASS
!
   integer :: i,len,mod
!
   integer :: a(5),b(5)
   integer :: nick
!
   NPASS = 1
!
   mod = 10; len = 5
!
   a = (/ 5,1,8,9,2 /)
   call encodg(a,mod,len, nick)
   call decodg(nick,mod,len, b)
!
   do i=1,len
      if (a(i) .ne. b(i)) NPASS = 0
   enddo
!
!
   if (NPASS.ne.1) stop 1
!
end program test_encodg_decodg
