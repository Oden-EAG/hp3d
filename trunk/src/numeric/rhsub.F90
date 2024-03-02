!-----------------------------------------------------------------------
!
!   routine name       - rhsub
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - backward gauss substitution
!
!   arguments in: a    - m by m matrix
!                 b    - vector of length m containing the right
!                        hand side
!                 n    - row dimension of matrix a exactly as
!                        specified in the dimension statement in the
!                        calling program
!                 m    - number of equations
!             out:x    - vector of length m containing the solution
!
!-----------------------------------------------------------------------
!
   subroutine rhsub(a,x,b,n,m)
!
      implicit none
!
      integer, intent(in) :: n,m
      real(8) :: a(n,n),x(*),b(*)
!
      integer :: i,j,j1,m1,ib
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      m1=m-1
!
!.....begin forward reduction of right hand side
      do i=1,m1
        j1=i+1
        do j=j1,m
          b(j) = b(j)-b(i)*a(j,i)/a(i,i)
        enddo
      enddo
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001)
 7001   format('rhsub: RHS AFTER FORWARD ELIMINATION')
        write(*,7002) b(1:m)
 7002   format(12e10.3)
      endif
#endif
!
!.....begin back substitution
      x(m)=b(m)/a(m,m)
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7003) m,b(m),a(m,m),x(m)
 7003   format('rhsub: i,b(i),a(i,i),x(i) = ',i3,3e12.5)
      endif
#endif
      do i=1,m1
        ib=m-i
        j1=ib+1
        do j=j1,m
          b(ib) = b(ib) - a(ib,j)*x(j)
#if HP3D_DEBUG
          if (iprint.eq.1.and.ib.eq.9) then
            write(*,*) 'ib,j,b(ib) = ',ib,j,b(ib)
          endif
#endif
        enddo
        x(ib)=b(ib)/a(ib,ib)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7003) ib,b(ib),a(ib,ib),x(ib)
        endif
#endif
      enddo
!
   end subroutine rhsub
