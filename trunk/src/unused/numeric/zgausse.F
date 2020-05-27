c-----------------------------------------------------------------------
c
c   routine name       - zgausse
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Dec 97
c
c   purpose            - solution of a full matrix linear system
c                        of equations using gauss elimination
c                        without pivoting - complex rhs case
c
c   usage              - call gausse(gk,igk,zgf,zu,n)
c
c   arguments in: gk   - n by n matrix
c                 igk  - row dimension of matrix gk exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 zgf   - vector of length n (the right-hand side)
c                 n    - number of equations
c             out:zu    - vector of length n containing the solution
c
c   required routines  - rhsub,tri
c
c-----------------------------------------------------------------------
c
      subroutine zgausse(gk,igk,zgf,zu,n)
c
#include "syscom.blk"
c
c
      dimension gk(igk,*),zgf(*),zu(*)
c
      call tri(gk,igk,n)

      call zrhsub(gk,zu,zgf,igk,n)
c
      return
      end

