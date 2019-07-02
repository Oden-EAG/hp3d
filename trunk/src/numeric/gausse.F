c-----------------------------------------------------------------------
c
c   routine name       - gausse
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - april 2, 1988
c
c   purpose            - solution of a full matrix linear system
c                        of equations using gauss elimination
c                        without pivoting
c
c   usage              - call gausse(gk,igk,gf,u,n)
c
c   arguments in: gk   - n by n matrix
c                 igk  - row dimension of matrix gk exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 gf   - vector of length n (the right-hand side)
c                 n    - number of equations
c             out:u    - vector of length n containing the solution
c
c   required routines  - rhsub,tri
c
c-----------------------------------------------------------------------
c
      subroutine gausse(gk,igk,gf,u,n)
c
#include "syscom.blk"
c
c
      dimension gk(igk,*),gf(*),u(*)
c
      call tri(gk,igk,n)
      call rhsub(gk,u,gf,igk,n)
c
      return
      end
