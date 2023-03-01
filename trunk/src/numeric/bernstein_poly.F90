!----------------------------------------------------------------------
!> @date Feb 2023
!----------------------------------------------------------------------
subroutine bernstein_poly(I,N,T, Poly,dPoly)
!----------------------------------------------------------------------
   implicit none
!----------------------------------------------------------------------
!..DUMMY ARGUMENTS
   integer, intent(in)  :: I,N
   real(8), intent(in)  :: T
   real(8), intent(out) :: Poly,dPoly
!----------------------------------------------------------------------
!..EXTERNAL FUNCTION
   integer,external :: fact
!----------------------------------------------------------------------
!
   if (I.lt.0) then
      Poly = 0.d0; dPoly = 0.d0
   elseif (I.gt.N) then
      Poly = 0.d0; dPoly = 0.d0
   else
      Poly = fact(N)/(fact(I)*fact(N-I)) * T**I * (1.d0 - T)**(N-I)
!.....set derivative equal zero for now...
      dPoly = 0.d0
   endif
!
end subroutine bernstein_poly
