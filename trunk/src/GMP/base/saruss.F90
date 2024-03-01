!-----------------------------------------------------------------------
!
!   routine name       - saruss
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - solve a linear 3 times 3 system of equations
!                        using determinants (Cramer's rule)
!
!   arguments :
!     in:
!                A     - lhs matrix
!                B     - rhs vector
!
!     out:
!                C     - the solution
!                Nfl   - return flag
!                        = 0 if non-singular system
!                        = 1 if numerically singular system
!
!-----------------------------------------------------------------------
!
   subroutine saruss(A,B, C,Nfl)
!
      implicit none
!
      integer :: Nfl
      real(8) :: A(3,3),B(3),C(3)
!
      real(8) :: detA,detX,detY,detZ
      integer :: i,j
!
      real(8), parameter :: eps = 1.0d-18
!
      Nfl=0
!
!  ...compute the main determinant
      detA = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
           + A(1,3)*A(2,1)*A(3,2) - (A(1,3)*A(2,2)*A(3,1) &
           + A(2,3)*A(3,2)*A(1,1 )+A(3,3)*A(1,2)*A(2,1))
!
!
      if (abs(B(1)).le.eps.and.abs(B(2)).le.eps.and.abs(B(3)).le.eps) &
        then
        C(1:3) = 0.d0
        return
      endif
!
      if (abs(detA).le.eps) then
        Nfl=1
        write(*,7001) detA
 7001   format('saruss: SINGULAR MATRIX detA = ',e12.5)
        do i=1,3
          write(*,7002) (A(i,j),j=1,3)
 7002     format(3e12.5)
        enddo
        call pause
        return
      endif
!
!
      detX = B(1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*B(3) &
            +A(1,3)*B(2)*A(3,2)-(A(1,3)*A(2,2)*B(3) &
            +A(2,3)*A(3,2)*B(1)+A(3,3)*A(1,2)*B(2))
!
      detY = A(1,1)*B(2)*A(3,3)+B(1)*A(2,3)*A(3,1) &
            +A(1,3)*A(2,1)*B(3)-(A(1,3)*B(2)*A(3,1) &
            +A(2,3)*B(3)*A(1,1)+A(3,3)*B(1)*A(2,1))
!
      detZ = A(1,1)*A(2,2)*B(3)+A(1,2)*B(2)*A(3,1) &
            +B(1)*A(2,1)*A(3,2)-(B(1)*A(2,2)*A(3,1) &
            +B(2)*A(3,2)*A(1,1)+B(3)*A(1,2)*A(2,1))
!
      C(1) = detX/detA
      C(2) = detY/detA
      C(3) = detZ/detA
!
!
   end subroutine saruss
