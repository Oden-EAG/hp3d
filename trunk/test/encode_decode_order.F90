!> @date Mar 2023
program test_encode_decode_order
!
   use parameters, only : MAXP
!
   implicit none
!
   integer :: NPASS
!
   integer :: i,len,mod
!
   integer :: nordx,nordy,nordz,norder
!
   integer, parameter :: px = min(2,MAXP)
   integer, parameter :: py = min(5,MAXP)
   integer, parameter :: pz = min(3,MAXP)
!
   NPASS = 1
!
   call encode_orderb(px,py,pz, norder)
   call decode_orderb(norder, nordx,nordy,nordz)
!
   if (nordx .ne. px) NPASS = 0
   if (nordy .ne. py) NPASS = 0
   if (nordz .ne. pz) NPASS = 0
!
   call encode_orderq(px,py, norder)
   call decode_orderq(norder, nordx,nordy)
!
   if (nordx .ne. px) NPASS = 0
   if (nordy .ne. py) NPASS = 0
!
   if (NPASS .ne. 1) then
      write(*,*) 'test_encode_decode_order FAILED.'
   else
      write(*,*) 'test_encode_decode_order PASSED.'
   endif
!
end program test_encode_decode_order
