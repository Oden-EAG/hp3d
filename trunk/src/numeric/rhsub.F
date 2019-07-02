c-----------------------------------------------------------------------
c
c   routine name       - rhsub
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - april 2, 1988
c
c   purpose            - backward gauss substitution
c
c   usage              - call rhsub(a,x,b,n,m)
c
c   arguments in: a    - m by m matrix
c                 b    - vector of length m containing the right
c                        hand side
c                 n    - row dimension of matrix a exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 m    - number of equations
c             out:x    - vector of length m containing the solution
c
c   required routines  - none
c
c-----------------------------------------------------------------------
c
      subroutine rhsub(a,x,b,n,m)
c
#include "syscom.blk"
      common /crhsub/ iprint
c
      dimension a(n,n),x(*),b(*)
c
      m1=m-1
c
c.....begin forward reduction of right hand side
      do i=1,m1
        j1=i+1
        do j=j1,m
          b(j) = b(j)-b(i)*a(j,i)/a(i,i)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,7001)
 7001   format('rhsub: RHS AFTER FORWARD ELIMINATION')
        write(*,7002) b(1:m)
 7002   format(12e10.3)
      endif
c
c.....begin back substitution
      x(m)=b(m)/a(m,m)
      if (iprint.eq.1) then
        write(*,7003) m,b(m),a(m,m),x(m)
 7003   format('rhsub: i,b(i),a(i,i),x(i) = ',i3,3e12.5)
      endif
      do i=1,m1
        ib=m-i
        j1=ib+1
        do j=j1,m
          b(ib) = b(ib) - a(ib,j)*x(j)
          if (iprint.eq.1.and.ib.eq.9) then
            write(*,*) 'ib,j,b(ib) = ',ib,j,b(ib)
          endif
        enddo
        x(ib)=b(ib)/a(ib,ib)
        if (iprint.eq.1) then
          write(*,7003) ib,b(ib),a(ib,ib),x(ib)
        endif
      enddo
c
c
      end
