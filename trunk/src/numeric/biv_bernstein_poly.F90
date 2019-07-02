!----------------------------------------------------------------------      
      subroutine biv_bernstein_poly(I,J,N,U,V, Poly,dPoly)
!----------------------------------------------------------------------
      implicit none
!----------------------------------------------------------------------
!     DUMMY ARGUMENTS
      integer,              intent(in)  :: I,J,N      
      real*8,               intent(in)  :: U,V
      real*8,               intent(out) :: Poly
      real*8, dimension(2), intent(out) :: dPoly
!----------------------------------------------------------------------
!     EXTERNAL FUNCTION
      integer :: fact
!----------------------------------------------------------------------
!     
      if ((I + J).gt.N) then
        write(*,*)'biv_Bernstein_poly: invalid values for I,J,N!'
        write(*,7000)I,J,N
 7000   format(' I = ',i2,'; J = ',i2,'; N = ',i2) 
        stop      
      endif
!
      Poly = fact(N)/(fact(I)*fact(J)*fact(N-I-J))*U**I*V**J*  &
             (1.d0-U-V)**(N-I-J)
      dPoly = 0.d0
!
      end subroutine biv_bernstein_poly
