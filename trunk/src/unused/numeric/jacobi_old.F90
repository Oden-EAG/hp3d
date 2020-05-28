!>  Purpose :
!!    build_block of weighted integrated jacobi polynomials and
!!    their derivatives
!!
!! @param Nord order of polynomial
!! @param A     weights
!! @param X     variable
!! @param Poly  polynomials
!! @param Dpoly derivatives
!! @revision Dec 10
subroutine build_jacobi_w_i (A,X,Nord, Poly,Dpoly)

  use error

  implicit none

  ! functions in use
  real, external ::  jacobi_w_i
  real, external :: djacobi_w_i

  ! arguments in out
  integer, intent(in) &
       :: Nord
  real*8, intent(in) &
       :: A, X
  real*8, intent(out),dimension(1:) &
       :: Poly, Dpoly

  ! local variables
  integer :: n

  ! == Pre-condition
  if (Nord.gt.size(Poly)) then
     write(*,*)'Nord = ',Nord, ' Size = ', size(Poly)
     call logic_error(ERR_OUT_OF_RANGE,__FILE__,__LINE__)
  endif

  ! == Body
  do n=1,Nord
     Poly(n)  =  jacobi_w_i(A,X,n)
     Dpoly(n) = djacobi_w_i(A,X,n)
  enddo
  ! == Post-condition

end subroutine build_jacobi_w_i



!----------------------------------------------------------------------
!
!  2-terms recurrence relation from Wolfram
!
!  c1 P_n+1^(a,b) = c2 P_n^(a,b) + c3 P_n-1^(a,b)
!       P_0^(a,b) = 1
!       P_1^(a,b) = [2(a + 1) + (a + b + 2) (x - 1)]/2
!
!  c1 = 2(n + 1) (n + a + b + 1) (2n + a + b)
!  c2 = (2n + a + b + 1) (a^2 - b^2) (2n + a + b + 2) *
!       (2n + a + b + 1) (2n + a + b)
!  c3 = -2(n + a) (n + b) (2n + a + b + 2)
!
!----------------------------------------------------------------------
!>  Purpose :
!!    evaluate weighted jacobi polynomial
!!
!! @param Nord order of polynomial
!! @param A     weights
!! @param B     weights
!! @param X     variable
!! @revision Dec 10
real*8 function jacobi_w (A,B,X,Nord)
  use error

  implicit none

  ! arguements in out
  integer, intent(in) &
       :: Nord
  real*8, intent(in) &
       :: A, B, X

  ! local variables
  real*8, dimension(0:Nord) &
       :: poly
  real*8 &
       :: c1, c2, c3, c4
  integer &
       :: iprint=0, n

  ! P_0^(a,b), P_1^(a,b)
  poly(0) = 1.d0
  poly(1) = (2.d0*(A+1.d0) + (A+B+2.d0)*(X+1.d0))/2.d0

  ! P_n^(a,b) , n >=2
  do n=1,Nord-1
     !  .....coefficients
     c1 = 2.d0*(n+1.d0)*(n+A+B+1.d0)*(2.d0*n+A+B)
     c2 = (2.d0*n+A+B+1.d0)*(A**2-B**2)
     c3 = (2.d0*n+A+B+2.d0)*(2.d0*n+A+B+1.d0)*(2.d0*n+A+B)
     c4 = -2.d0*(n+A)*(n+B)*(2.d0*n+A+B+2.d0)
     !  .....recurrence relation
     poly(n+1) = ((c2 + c3*X)*poly(n) + c4*poly(n-1))/c1
  enddo

  !
  jacobi_w = poly(Nord)
  !
end function jacobi_w

!>  Purpose :
!!    evaluate derivative of weighted jacobi polynomial
!!
!! @param Nord order of polynomial
!! @param A     weights
!! @param B     weights
!! @param X     variable
!! @revision Dec 10
real*8 function djacobi_w(A,B,X,Nord)
  !----------------------------------------------------------------------
  !
  !   for Wolfram (and many other sources)
  !
  !   P_n^(a,b)' = (n + a + b + 1)/2 P_n-1^(a+1,b+1)
  !
  !----------------------------------------------------------------------
  use error

  implicit none

  ! functions in use
  real, external :: jacobi_w

  ! arguements in out
  integer, intent(in) &
       :: Nord
  real*8, intent(in) &
       :: A, B, X

  ! local variables
  integer &
       :: iprint = 0

  if (Nord.eq.0) then
     djacobi_w = 0.d0
  else
     djacobi_w = jacobi_w(A+1,B+1,X,Nord-1)*(Nord+A+B+1.d0)/2.d0
  endif

end function djacobi_w

!>  Purpose :
!!    evaluate weighted integrated jacobi polynomial
!!
!! @param Nord order of polynomial
!! @param A     weights
!! @param X     variable
!! @revision Dec 10
real*8 function jacobi_w_i(A,X,Nord)
  !----------------------------------------------------------------------
  !   p_n^a := P_n^(a,0)
  !
  !   \hat{p}_0^a = 1
  !   \hat{p}_n^a = \int_{-1}^x p_{n-1}^a , n >= 1
  !   Recurrence relations from Beuchler's papers:
  !
  !   n >= 2
  !   \hat{p}_n^a = c1 p_n^a + c2 p_{n-1}^a + c3 p_{n-2}^a
  !
  !   c1 = 2(n + a)/[(2n + a - 1) (2n + a)]
  !   c2 = 2a/[(2n + a - 2) (2n + a)]
  !   c3 = -2(n - 1)/[(2n + a - 1) (2n + a - 2)]
  !
  !   n = 1
  !   \hat{p}_n^a = c1 ( p_n^{a-1} + p_{n-1}^{a-1} )
  !   c1 = 1/(2n + a - 1)
  !----------------------------------------------------------------------
  use error

  implicit none
  ! functions in use
  real, external :: jacobi_w
  ! arguements in out
  integer, intent(in) &
       :: Nord
  real*8, intent(in) &
       :: A, X

  ! local variables
  real*8 &
       :: c1,c2,c3
  integer &
       :: iprint = 0

  select case (Nord)
  case(0)
     Jacobi_W_I = 1.d0
     !  ...use a 2 terms recurrence formula for \hat{p}_1^a
  case(1)
     c1 = 1.d0/(2.d0*Nord + A - 1.d0)

     jacobi_w_i = c1*jacobi_w(A-1.d0,0.d0,X,Nord) &
          + c1*jacobi_w(A-1.d0,0.d0,X,Nord-1)

     !  ...use a 3 terms recurrence formula for \hat{p}_n^a , n >=2
  case default
     c1 = 2.d0*(A+Nord)/((2.d0*Nord+A-1.d0)*(2.d0*Nord+A))
     c2 = 2.d0*A/((2.d0*Nord+A-2.d0)*(2.d0*Nord+A))
     c3 = -2.d0*(Nord-1.d0)/((2.d0*Nord+A-1.d0)*(2.d0*Nord+A-2.d0))

     Jacobi_W_I = c1*Jacobi_W(A,0.d0,X,Nord) + &
                 c2*Jacobi_W(A,0.d0,X,Nord-1) + &
                 c3*Jacobi_W(A,0.d0,X,Nord-2)
  end select

end function jacobi_w_i

!>  Purpose :
!!    evaluate derivative of weighted integrated jacobi polynomial
!!
!! @param Nord order of polynomial
!! @param A     weights
!! @param X     variable
!! @revision Dec 10
real*8 function djacobi_w_i(A,X,Nord)
  use error

  implicit none

  ! functions in use
  real, external :: djacobi_w

  ! arguements in out
  integer, intent(in) &
       :: Nord
  real*8, intent(in) &
       :: A, X

  ! local variables
  real*8  &
       :: c1, c2, c3
  integer :: &
       iprint = 0

  select case (Nord)
  case(0)
     dJacobi_W_I = 0.d0

     !  ...use a 2 terms recurrence formula for \hat{p}_1^a
  case(1)
     c1 = 1.d0/(2.d0*Nord+A-1.d0)

     dJacobi_W_I = c1*dJacobi_W(A-1.d0,0.d0,X,Nord) &
          + c1*dJacobi_W(A-1.d0,0.d0,X,Nord-1)

     !  ...use a 3 terms recurrence formula for \hat{p}_n^a , n >=2
  case default
     c1 = 2.d0*(A+Nord)/((2.d0*Nord+A-1.d0)*(2.d0*Nord+A))
     c2 = 2.d0*A/((2.d0*Nord+A-2.d0)*(2.d0*Nord+A))
     c3 = -2.d0*(Nord-1.d0)/((2.d0*Nord+A-1.d0)*(2.d0*Nord+A-2.d0))

     dJacobi_W_I = c1*dJacobi_W(A,0.d0,X,Nord) + &
          c2*dJacobi_W(A,0.d0,X,Nord-1) + &
          c3*dJacobi_W(A,0.d0,X,Nord-2)
  end select

end function djacobi_w_i
