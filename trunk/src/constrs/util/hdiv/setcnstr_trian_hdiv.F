c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian_hdiv
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
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
      implicit none
c
      RRTQ = 0.d0
c
      call setcnstr_trian_iso_hdiv
c
      end subroutine setcnstr_trian_hdiv
