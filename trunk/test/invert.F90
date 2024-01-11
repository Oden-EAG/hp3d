!> @date Mar 2023
program test_invert
!
   implicit none
!
   integer :: NPASS
!
   real(8)    :: A(3,3),invA(3,3),detA
   real(8)    :: mat(3,3),det
#if C_MODE
   complex(8) :: B(3,3),invB(3,3),detB
   complex(8) :: zmat(3,3),zdet
#else
   real(8)    :: B(3,3),invB(3,3),detB
   real(8)    :: zmat(3,3),zdet
#endif
!
   integer :: i,j
!
   logical, external :: dnear,znear
!
   NPASS = 1
!
   A    = reshape( (/  1.d0, -1.d0,  3.d0,   &
                      -2.d0,  2.d0,  1.d0,   &
                      -1.d0,  0.d0,  3.d0 /) &
                  ,(/ 3,3 /))
   invA = reshape( (/ 8.571428571428572d-01, &
                      4.285714285714285d-01, &
                     -1.0d0                , &
                      7.142857142857144d-01, &
                      8.571428571428571d-01, &
                     -1.0d0                , &
                      2.857142857142857d-01, &
                      1.428571428571428d-01, &
                      0.0d0               /) &
                  ,(/ 3,3 /))
   detA = 7.d0
!
#if C_MODE
   B    = reshape( (/ ( 1.d0, 1.d0), (-1.d0, 1.d0), (3.d0, 0.d0),   &
                      (-2.d0,-1.d0), ( 2.d0, 2.d0), (1.d0, 0.d0),   &
                      (-1.d0,-1.d0), ( 0.d0, 3.d0), (3.d0,-2.d0) /) &
                  ,(/ 3,3 /))
   invB = reshape( (/ ( 3.739837398373985d-01,-3.658536585365854d-01),   &
                      ( 1.707317073170732d-01, 1.300813008130082d-01),   &
                      (-4.471544715447155d-01, 2.439024390243901d-02),   &
                      ( 2.195121951219513d-01,-3.089430894308944d-01),   &
                      ( 4.552845528455285d-01,-9.756097560975605d-02),   &
                      (-4.146341463414635d-01, 6.504065040650407d-02),   &
                      ( 5.691056910569101d-02,-1.788617886178863d-01),   &
                      ( 1.056910569105691d-01,-2.845528455284553d-01),   &
                      ( 4.065040650406504d-02, 3.008130081300814d-01) /) &
                  ,(/ 3,3 /))
   detB = (15.d0,12.d0)
#else
   B = A; invB = invA; detB = detA
#endif
!
!  Test: get_ddet
   call get_ddet(A, det)
   if (.not. dnear(det,detA)) NPASS = 0
!
!  Test: get_zdet
   call get_zdet(B, zdet)
   if (.not. znear(zdet,detB)) NPASS = 0
!
!  Test: dinvert
   call dinvert(A, mat,det)
   if (.not. dnear(det,detA)) NPASS = 0
   do j=1,3
      do i=1,3
         if (.not. dnear(mat(i,j),invA(i,j))) NPASS = 0
      enddo
   enddo
!
!  Test: zinvert
   call zinvert(B, zmat,zdet)
   if (.not. znear(zdet,detB)) NPASS = 0
   do j=1,3
      do i=1,3
         if (.not. znear(zmat(i,j),invB(i,j))) NPASS = 0
      enddo
   enddo
!
!
   if (NPASS.ne.1) stop 1
!
end program test_invert
