c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian
c
c----------------------------------------------------------------------
c
c   latest revision    - May 12
c
c   purpose            - routine evaluates constraint coefficients for
c                        the triangular master element
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine setcnstr_trian_hcurl
c      
      use parameters
      use constraints
#include "syscom.blk"
c
      RRTE = 0.d0; RRQE = 0.d0
c
      call setcnstr_trian_iso_hcurl
ccc      do i=2,4
ccc        call setcnstr_trian_aniso_hcurl(i)
ccc      enddo
c      
      end subroutine setcnstr_trian_hcurl
