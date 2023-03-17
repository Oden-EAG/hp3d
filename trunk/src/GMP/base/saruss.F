c-----------------------------------------------------------------------
c
c   routine name       - saruss
c
c-----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - solve a linear 3 times 3 system of equations
c                        using determinants (Cramer's rule)
c
c   arguments :
c     in:
c                A     - lhs matrix
c                B     - rhs vector
c
c     out:
c                C     - the solution
c                Nfl   - return flag
c                        = 0 if non-singular system
c                        = 1 if numerically singular system
c
c-----------------------------------------------------------------------
c
      subroutine saruss(A,B, C,Nfl)
c
      implicit none
c
      integer :: Nfl
      real(8) :: A(3,3),B(3),C(3)
c
      real(8) :: detA,detX,detY,detZ
      integer :: i,j
c
      real(8), parameter :: eps = 1.0d-18
c
      Nfl=0
c
c  ...compute the main determinant
      detA = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     .     + A(1,3)*A(2,1)*A(3,2) - (A(1,3)*A(2,2)*A(3,1)
     .     + A(2,3)*A(3,2)*A(1,1 )+A(3,3)*A(1,2)*A(2,1))
c
c
      if (abs(B(1)).le.eps.and.abs(B(2)).le.eps.and.abs(B(3)).le.eps)
     .  then
        C(1:3) = 0.d0
        return
      endif
c
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
c
c
      detX = B(1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*B(3)
     .      +A(1,3)*B(2)*A(3,2)-(A(1,3)*A(2,2)*B(3)
     .      +A(2,3)*A(3,2)*B(1)+A(3,3)*A(1,2)*B(2))
c
      detY = A(1,1)*B(2)*A(3,3)+B(1)*A(2,3)*A(3,1)
     .      +A(1,3)*A(2,1)*B(3)-(A(1,3)*B(2)*A(3,1)
     .      +A(2,3)*B(3)*A(1,1)+A(3,3)*B(1)*A(2,1))
c
      detZ = A(1,1)*A(2,2)*B(3)+A(1,2)*B(2)*A(3,1)
     .      +B(1)*A(2,1)*A(3,2)-(B(1)*A(2,2)*A(3,1)
     .      +B(2)*A(3,2)*A(1,1)+B(3)*A(1,2)*A(2,1))
c
      C(1) = detX/detA
      C(2) = detY/detA
      C(3) = detZ/detA
c
c
      end subroutine saruss
