!-----------------------------------------------------------------------
!
!   routine name       - tri
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - forward Gauss elimination
!
!   arguments
!      in:        n    - row dimension of matrix "a" exactly as
!                        specified in the dimension statement in the
!                        calling program
!                 m    - number of equations
!     inout:      a    - Matrix which is Gauss factorized in-place
!
!-----------------------------------------------------------------------
!
      subroutine tri(a,n,m)
!
      implicit none
!
      integer, intent(in)    :: n,m
      real(8), intent(inout) :: a(n,n)
!
      integer :: i,j,j1,k,m1
      real(8) :: fac,tiny
!
      real(8), parameter :: eps = 1.0d-15
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      m1=m-1
      tiny=1.e-30
!
!.....eliminate degree of freedom i
      do i=1,m1
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7001) i,a(i,i)
 7001     format('tri: i,a(i,i) = ',i3,e12.5)
        endif
#endif
!
!.......check for excessively small pivot
        if (abs(a(i,i)).lt.tiny) then
          write(*,9999) i,a(i,i)
 9999     format('tri: reduction failed due to small pivot', &
            /,' equation no.',i5,' pivot ',e12.5)
          stop 1
        endif
        j1=i+1
!
!.......modify rows j
        do j=j1,m
!
!  .......skip if the leading coefficient is zero
          if (dabs(a(j,i)) < eps) cycle
          fac = a(j,i)/a(i,i)
          do k=j1,m
            a(j,k) = a(j,k) - a(i,k)*fac
          enddo
        enddo
      enddo
!
!  ...write out the last pivot
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) m,a(m,m)
      endif
#endif
!
      end subroutine tri
