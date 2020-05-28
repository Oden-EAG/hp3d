c-----------------------------------------------------------------------
c
c   routine name       - zrhsub
c
c-----------------------------------------------------------------------
c
c   latest revision    - april 2, 1988
c
c   purpose            - backward gauss substitution, complex rhs case
c
c   arguments in: a    - m by m matrix
c                 zb    - vector of length m containing the right
c                        hand side
c                 n    - row dimension of matrix a exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 m    - number of equations
c             out:zx   - vector of length m containing the solution
c
c-----------------------------------------------------------------------
c
      subroutine zrhsub(a,zx,zb,n,m)
c
#include "syscom.blk"
c
      dimension a(n,n),zx(*),zb(*)
c
      if (m.eq.0) return
c
      m1=m-1
c.....begin forward reduction of right hand side
      do 20  i=1,m1
             j1=i+1
             do 10 j=j1,m
                   zb(j)=zb(j)-zb(i)*a(j,i)/a(i,i)
   10              continue
c
 20          continue
c
c.....begin back substitution
      zx(m)=zb(m)/a(m,m)
c
      do 40 i=1,m1
         ib=m-i
         j1=ib+1
c
         do 30 j=j1,m
               zb(ib)=zb(ib)-a(ib,j)*zx(j)
   30          continue
c
               zx(ib)=zb(ib)/a(ib,ib)
40       continue
c
      end
