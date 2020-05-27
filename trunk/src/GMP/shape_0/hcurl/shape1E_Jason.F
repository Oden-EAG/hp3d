
c----------------------------------------------------------------------
c
c   routine name       - shape1E
c
c----------------------------------------------------------------------
c
c   latest revision    - Jul 06
c
c   purpose            - routine evaluates 1D hierarchical shape
c                        functions with their derivatives
c
c   arguments :
c     in:
c              Xi      - coordinate between 0 and 1
c              Norder  - order of approximation
c     out:
c              Nrdof   - number of dof (shape functions)
c              Vshap   - values of the shape functions at the point
c              Dvshap  - their derivatives
c
c----------------------------------------------------------------------
c
      subroutine shape1E_Jason(Xi,Norder,  Nrdof,Vshap)
c
      use parameters, only:  NPMAX => MAXP

#include "syscom.blk"
c
      dimension Vshap(*)
c
      dimension v0(NPMAX+2),d0(NPMAX+2)
c
      Nrdof = Norder
      call shape1_Jason(Xi,Norder, nvoid,v0,d0)
      do i=1,Nrdof
        Vshap(i) = d0(i+1)
      enddo
c
      end subroutine shape1E_Jason
