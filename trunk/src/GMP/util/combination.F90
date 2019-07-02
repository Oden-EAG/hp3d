integer function combination(n, m)
! how many ways can we choose m objects from n objects?
  integer :: factorial
!
  if (n .lt. m ) then
    write(*,*)'combinations: n should be greater than m; n, m = ',n,m
    stop
  endif        
  combination = factorial(n)/(factorial(m)*factorial(n - m))
end function combination
