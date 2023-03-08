c-----------------------------------------------------------------------
c
c   routine name       - tri
c
c-----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - forward Gauss elimination
c
c   arguments
c      in:        n    - row dimension of matrix "a" exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 m    - number of equations
c     inout:      a    - Matrix which is Gauss factorized in-place
c
c-----------------------------------------------------------------------
c
      subroutine tri(a,n,m)
c
      implicit none
c
      integer, intent(in)    :: n,m
      real(8), intent(inout) :: a(n,n)
c
      integer :: i,j,j1,k,m1
      real(8) :: fac,tiny
c
#if DEBUG_MODE
      integer :: iprint=0
#endif
c
      m1=m-1
      tiny=1.e-30
c
c.....eliminate degree of freedom i
      do i=1,m1
#if DEBUG_MODE
        if (iprint.eq.1) then
          write(*,7001) i,a(i,i)
 7001     format('tri: i,a(i,i) = ',i3,e12.5)
        endif
#endif
c
c.......check for excessively small pivot
        if (abs(a(i,i)).lt.tiny) then
          write(*,9999) i,a(i,i)
 9999     format('tri: reduction failed due to small pivot',
     .      /,' equation no.',i5,' pivot ',e12.5)
          stop 1
        endif
        j1=i+1
c
c.......modify rows j
        do j=j1,m
c
c  .......skip if the leading coefficient is zero
          if (a(j,i).eq.0.d0) cycle
          fac = a(j,i)/a(i,i)
          do k=j1,m
            a(j,k) = a(j,k) - a(i,k)*fac
          enddo
        enddo
      enddo
c
c  ...write out the last pivot
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7001) m,a(m,m)
      endif
#endif
c
      end subroutine tri
