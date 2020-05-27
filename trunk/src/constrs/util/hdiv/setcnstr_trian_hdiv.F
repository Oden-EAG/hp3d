c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian_div
c
c----------------------------------------------------------------------
c
c   latest revision    - May 14
c
c   purpose            - routine evaluates constraint coefficients for
c                        the triangular master element and H(div) (L2)
c                        constrained appoximation
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine setcnstr_trian_hdiv
c
      use parameters
      use constraints
#include "syscom.blk"
c
      RRTQ = 0.d0
c
      call setcnstr_trian_iso_hdiv
ccc      do i=2,4
ccc        call setcnstr_trian_aniso_hdiv(i)
ccc      enddo
c
      end subroutine setcnstr_trian_hdiv
