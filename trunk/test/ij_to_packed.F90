!> @date Mar 2023
program test_ij_to_packed
!
   implicit none
!
   integer :: NPASS
!
   integer      :: i,j,n
   character(1) :: uplo
!
   integer, external :: ij_to_packed
!
   NPASS = 1
!
   uplo = 'U'
   n = 3; i = 2; j = 2
   if (ij_to_packed(i,j,n,uplo) .ne. 3) NPASS = 0
   n = 3; i = 2; j = 3
   if (ij_to_packed(i,j,n,uplo) .ne. 5) NPASS = 0
!
   uplo = 'L'
   n = 3; i = 2; j = 2
   if (ij_to_packed(i,j,n,uplo) .ne. 4) NPASS = 0
   n = 3; i = 3; j = 2
   if (ij_to_packed(i,j,n,uplo) .ne. 5) NPASS = 0
!
!
   if (NPASS.ne.1) stop 1
!
end program test_ij_to_packed
