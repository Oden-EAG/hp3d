!> @date Mar 2023
program test_near
!
   implicit none
!
   integer :: NPASS
!
   real(8) :: a,b
!
#if C_MODE
   complex(8) :: c,d
#else
   real(8) :: c,d
#endif
!
   logical, external :: dnear,znear
!
   NPASS = 1
!
!  Test: dnear
   a = 1.d0; b = 1.d0
   if (.not. dnear(a,b)) NPASS = 0
!
   b = a - 1.d-16
   if (.not. dnear(a,b)) NPASS = 0
!
   b = 1.01d0
   if (dnear(a,b)) NPASS = 0
!
!
!  Test: znear
#if C_MODE
   c = (1.d0,-1.d0)
   d = (1.d0,-1.d0)
   if (.not. znear(c,d)) NPASS = 0
!
   d = c - (1.d-16,-1.d-16)
   if (.not. znear(c,d)) NPASS = 0
!
   d = (1.d0,-1.01d0)
   if (znear(c,d)) NPASS = 0
#else
   c = 1.d0; d = 1.d0
   if (.not. znear(c,d)) NPASS = 0
!
   d = c - 1.0d-16
   if (.not. znear(c,d)) NPASS = 0
!
   d = 1.01d0
   if (znear(c,d)) NPASS = 0
#endif
!
   write(*,*) NPASS
!
   if (NPASS .ne. 1) then
      write(*,*) 'test_near FAILED.'
   else
      write(*,*) 'test_near PASSED.'
   endif
!
end program test_near
