!---------------------------------------------------------------------------------
!
      subroutine radial_derivative(s,j,k, c)
!
!---------------------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------------------
!     DUMMY ARGUMENTS
      real(8), intent(in)  :: s
      integer, intent(in)  :: j,k
      real(8), intent(out) :: c
!     LOCAL VARIABLES
      real(8)              :: poly1,poly2,poly3,dpoly
      integer              :: iprint
      integer,parameter    :: deg = 6
!---------------------------------------------------------------------------------
!
      iprint=0
!
      call Bernstein_poly(j-2,deg-k-2,s, poly1,dpoly)
      call Bernstein_poly(j-1,deg-k-2,s, poly2,dpoly)
      call Bernstein_poly(j  ,deg-k-2,s, poly3,dpoly)
      c = deg*(deg-k)*(deg-k-1)*(poly1 - 2.d0*poly2 + poly3)*(-1.d0)**(k+1)
      if (iprint.eq.1) then
        write(*,1000)k,j,c
 1000   format(' radial_derivative: k = ',i1,'; j = ',i1,'; --> c = ',e12.5)
      endif
!
      end subroutine radial_derivative
