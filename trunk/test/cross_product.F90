!> @date Mar 2023
program test_cross_product
!
   implicit none
!
   integer :: NPASS
!
   real(8) :: A(3),B(3)
   real(8) :: AB(3)
!
   real(8) :: F(2),G(2)
   real(8) :: fg
!
   real(8) :: VEC(3)
   real(8) :: sVEC(1)
   real(8) :: val
!
#if C_MODE
   complex(8) :: C(3),D(3)
   complex(8) :: AC(3),CD(3)
   complex(8) :: zVEC(3)
#else
   real(8)    :: C(3),D(3)
   real(8)    :: AC(3),CD(3)
   real(8)    :: zVEC(3)
#endif
!
   integer :: i
   logical, external :: dnear,znear
!
   NPASS = 1
!
!  Test: cross_product
   A  = (/  1.d0,  2.d0, 2.d0 /)
   B  = (/  2.d0, -3.d0,-4.d0 /)
   AB = (/ -2.d0,  8.d0,-7.d0 /)
!
   call cross_product(A,B, VEC)
!
   do i=1,3
      if (.not. dnear(VEC(i),AB(i))) NPASS = 0
   enddo
!
!  Test: zcross_product
#if C_MODE
   C  = (/ ( 2.d0,1.d0), ( -3.d0,-1.d0), (-4.d0,  1.d0) /)
   D  = (/ ( 1.d0,2.d0), (  2.d0, 1.d0), ( 2.d0, -2.d0) /)
   AC = (/ (-2.d0,4.d0), (  8.d0, 1.d0), (-7.d0, -3.d0) /)
   CD = (/ ( 1.d0,6.d0), (-12.d0,-5.d0), ( 4.d0, 11.d0) /)
#else
   C  = (/  2.d0, -3.d0,-4.d0 /)
   D  = (/  1.d0,  2.d0, 2.d0 /)
   AC = (/ -2.d0,  8.d0,-7.d0 /)
   CD = (/  2.d0, -8.d0, 7.d0 /)
#endif
!
   call zcross_product(A,C, zVEC)
   do i=1,3
      if (.not. znear(zVEC(i),AC(i))) NPASS = 0
   enddo
!
!  Test: zz_cross_product
   call zz_cross_product(C,D, zVEC)
   do i=1,3
      if (.not. znear(zVEC(i),CD(i))) NPASS = 0
   enddo
!
!  Test: cross_product2D
   F  = (/  1.d0, -4.d0 /)
   G  = (/  2.d0, -3.d0 /)
   fg = 5.d0
!
   call cross_product2D(F,G, val)
!
   if (.not. dnear(val,fg)) NPASS = 0
!
!  Test: cross
   call cross(3,A,B, VEC)
!
   do i=1,3
      if (.not. dnear(VEC(i),AB(i))) NPASS = 0
   enddo
!
   call cross(2,F,G, sVec)
!
   if (.not. dnear(sVec(1),fg)) NPASS = 0
!
!
   if (NPASS.ne.1) stop 1
!
end program test_cross_product
