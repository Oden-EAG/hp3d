c----------------------------------------------------------------------
c      
      subroutine biv_Bernstein_poly(U,V,W, BB)
c
c----------------------------------------------------------------------
      use kinds
c----------------------------------------------------------------------
      implicit none
c----------------------------------------------------------------------
c  ...set degree      
      integer, parameter :: deg = 4
c----------------------------------------------------------------------
c     DUMMY ARGUMENTS
      real(DP), intent(in)                                :: U,V,W
      real(DP), dimension(0:deg,0:deg,0:deg), intent(out) :: BB
c----------------------------------------------------------------------
c     EXTERNAL FUNCTION
      integer :: fact
c----------------------------------------------------------------------
c     LOCAL VARIABLES
      real(DP) :: c
      integer  :: i,j,k
c----------------------------------------------------------------------
c
      do i=0,deg
        do j=0,deg
          do k=0,deg
            c = fact(i+j+k)/(fact(i)*fact(j)*fact(k))
            BB(i,j,k) = c*U**i*V**j*W**k
          enddo
        enddo
      enddo
c
c
      end
