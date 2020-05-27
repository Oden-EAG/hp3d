integer function factorial(n)
! factorial of n
  integer :: i,n
!
  if (n .eq. 0) then
    factorial = 1
    return
  endif
  factorial = 1
  do i = 2, n
    factorial = factorial*i
  enddo
end function factorial
