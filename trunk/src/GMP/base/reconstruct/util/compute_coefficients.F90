!---------------------------------------------------------------------------------
!
      subroutine compute_coefficients(s,i,j, c) 
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------------------
!     DUMMY ARGUMENTS
      real(8), intent(in)  :: s
      integer, intent(in)  :: i,j
      real(8), intent(out) :: c
!     LOCAL VARIABLES      
      real(8)              :: poly1,poly2,poly3,dpoly
      integer              :: iprint
      integer,parameter    :: deg = 7
!---------------------------------------------------------------------------------
!
      iprint=0
!
      call Bernstein_poly(i-2,deg-j-2,s, poly1,dpoly)
      call Bernstein_poly(i-1,deg-j-2,s, poly2,dpoly)
      call Bernstein_poly(i  ,deg-j-2,s, poly3,dpoly)
      c = deg*(deg-j)*(deg-j-1)*(poly1 - 2.d0*poly2 + poly3)*(-1.d0)**(j+1)
      if (iprint.eq.1) then
        write(*,1000)i,j,c
 1000   format(' compute_coefficients: i = ',i1,'; j = ',i1,'; --> c = ',e12.5)          
      endif
!
      end subroutine


!---------------------------------------------------------------------------------
!
      subroutine rad_der_coefficient(s,i,j,d, c) 
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------------------
!     DUMMY ARGUMENTS
      real(8), intent(in)  :: s
      integer, intent(in)  :: i,j,d
      real(8), intent(out) :: c
!     LOCAL VARIABLES      
      real(8)              :: poly1,poly2,poly3,poly4,dpoly
      integer              :: iprint
      integer,parameter    :: deg = 7
!---------------------------------------------------------------------------------
!
      iprint=0

      if (j.ge.2) then
        c = 0.d0
        return
      endif

      select case(d)
!  ...radial derivative      
      case(0)
        call Bernstein_poly(i,deg-j,s, poly1,dpoly)
        c = deg*poly1*(-1.d0)**(j+1)
!  ...1st derivative
      case(1)
        call Bernstein_poly(i-1,deg-j-1,s, poly1,dpoly)
        call Bernstein_poly(i  ,deg-j-1,s, poly2,dpoly)
        c = deg*(deg-j)*(poly1 - poly2)*(-1.d0)**(j+1)
!  ...2nd derivative              
      case(2)
        call Bernstein_poly(i-2,deg-j-2,s, poly1,dpoly)
        call Bernstein_poly(i-1,deg-j-2,s, poly2,dpoly)
        call Bernstein_poly(i  ,deg-j-2,s, poly3,dpoly)
         c = deg*(deg-j)*(deg-j-1)*(poly1 - 2.d0*poly2 + poly3)*(-1.d0)**(j+1)
!  ...3rd derivative         
      case(3)
        call Bernstein_poly(i-3,deg-j-3,s, poly1,dpoly)
        call Bernstein_poly(i-2,deg-j-3,s, poly2,dpoly)
        call Bernstein_poly(i-1,deg-j-3,s, poly3,dpoly)
        call Bernstein_poly(i  ,deg-j-3,s, poly4,dpoly)
        c = deg*(deg-j)*(deg-j-1)*(deg-j-2)*  &
            (poly1 - 3.d0*poly2 + 3.d0*poly3 - poly4)*(-1.d0)**(j+1)
      endselect

      if (iprint.eq.1) then
        write(*,1000)d,i,j,c
 1000   format(' rad_der_coefficient: d = ',i1,'; i = ',i1,'; j = ',i1,'; --> c = ',e12.5)          
      endif
!
      end subroutine
