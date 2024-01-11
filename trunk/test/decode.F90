!> @date Mar 2023
program test_decode
!
   implicit none
!
   integer :: NPASS
!
   integer :: nx,ny,nz
   integer :: nick
!
   NPASS = 1
!
!  Test: ddecode
   nick = 538;
!
   call ddecode(nick, nx,ny,nz)
!
   if (nx .ne. 5) NPASS = 0
   if (ny .ne. 3) NPASS = 0
   if (nz .ne. 8) NPASS = 0
!
!  Test: decode
   nick = 19
!
   call decode(nick, nx,ny)
!
   if (nx .ne. 1) NPASS = 0
   if (ny .ne. 9) NPASS = 0
!
!  Test: decode2
   nick = 81762
!
   call decode2(nick, nx,ny)
!
   if (nx .ne. 817) NPASS = 0
   if (ny .ne.  62) NPASS = 0
!
!
   if (NPASS.ne.1) stop 1
!
end program test_decode
