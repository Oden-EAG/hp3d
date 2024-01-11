c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian_hcurl
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
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
      implicit none
c
      RRTE = 0.d0
c
      call setcnstr_trian_iso_hcurl
c
      end subroutine setcnstr_trian_hcurl
