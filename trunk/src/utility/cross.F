c
c----------------------------------------------------------------------
c
c   routine name       - cross
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 13
c
c   purpose            - routine evaluates cross product of two
c                        vectors in R^2 or R^3
c
c   arguments :
c     in:    
c           N          - space dimension N=2,3
c           A,B        - vectors in R^N
c     out: 
c           C          - scalar (R^2) or vector (R^3)
c
c----------------------------------------------------------------------
c
      subroutine cross(N,A,B, C)
c
#include "syscom.blk"
c
      dimension A(N),B(N), C(2*N-3)
c
      if (N.eq.3) then    
        C(1) =   A(2)*B(3) - A(3)*B(2)
        C(2) = - A(1)*B(3) + A(3)*B(1)
      endif
      C(2*N-3) =   A(1)*B(2) - A(2)*B(1)
c
      end subroutine cross
