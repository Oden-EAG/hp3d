!> @date Mar 2023
program test_encod_decod
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
   integer    :: c(10),d(10)
   integer(8) :: nickLong
!
   NPASS = 1
!
   mod = 10; len = 5
!
   a = (/ 5,1,8,9,2 /)
   call encod(a,mod,len, nick)
   call decod(nick,mod,len, b)
!
   do i=1,len
      if (a(i) .ne. b(i)) NPASS = 0
   enddo
!
   mod = 10; len = 10
!
   c = (/ 5,1,8,9,2, 6,2,4,3,9 /)
   call encodLong(c,mod,len, nickLong)
   call decodLong(nickLong,mod,len, d)
!
   do i=1,len
      if (c(i) .ne. d(i)) NPASS = 0
   enddo
!
   if (NPASS .ne. 1) then
      write(*,*) 'test_encod_decod FAILED.'
   else
      write(*,*) 'test_encod_decod PASSED.'
   endif
!
end program test_encod_decod
