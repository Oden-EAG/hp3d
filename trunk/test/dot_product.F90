!> @date Mar 2023
program test_dot_product
!
   implicit none
!
   integer :: NPASS
!
   real(8) :: A(3),B(3)
   real(8) :: ab
!
#if C_MODE
   complex(8) :: C(3)
   complex(8) :: ac,val
#else
   real(8) :: C(3)
   real(8) :: ac,val
#endif
!
   logical, external :: dnear,znear
!
   NPASS = 1
!
!  Test: dot_product
   A = (/ 1.d0, 2.d0, 3.d0 /)
   B = (/ 4.d0,-3.d0, 1.d0 /)
!
   call dot_product(A,B, ab)
!
   if (.not. dnear(ab,1.d0)) NPASS = 0
!
!  Test: zdot_product
#if C_MODE
   C = (/ (4.d0,-4.d0), (-3.d0,3.0d0), (1.d0,-1.d0) /)
   val = (1.d0,-1.d0)
#else
   C = (/ 4.d0,-3.d0, 1.d0 /)
   val = 1.d0
#endif
!
   call zdot_product(A,C, ac)
!
   if (.not. znear(ac,val)) NPASS = 0
!
!
   if (NPASS.ne.1) stop 1
!
end program test_dot_product
