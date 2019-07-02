      subroutine normal_derivative(s, coeff)
      implicit none
      real*8,external :: Bern_poly
      real*8 :: c,s,poly1,poly2,poly3,dpoly
      integer :: deg
      integer :: j,k
      real*8,dimension(0:1,0:5) :: coeff

      deg = 6
!      write(*,*)'bern_poly',Bern_poly(6,5,0.2d0)

      do k = 0,1
        do j = 0,(deg-k)
          call Bernstein_poly(j-2,deg-k-2,s, poly1,dpoly)
          call Bernstein_poly(j-1,deg-k-2,s, poly2,dpoly)
          call Bernstein_poly(j  ,deg-k-2,s, poly3,dpoly)

          c = deg*(deg-k)*(deg-k-1)*(     Bern_poly(j-2,deg-k-2,s) -                &
                                     2.d0*Bern_poly(j-1,deg-k-2,s) +                &
                                          Bern_poly(j  ,deg-k-2,s))*(-1.d0)**(k+1)
!          c = deg*(deg-k)*(deg-k-1)*(poly1 - 2.d0*poly2 + poly3)*(-1.d0)**(k+1)
          coeff(k,j) = c

!          write(*,1000)k,j,c
 1000     format(' k = ',i1,'; j = ',i1,'; --> c = ',e12.5)
        enddo
      enddo



      end
