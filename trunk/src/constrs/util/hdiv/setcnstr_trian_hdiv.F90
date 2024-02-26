!----------------------------------------------------------------------
!
!   routine name       - setcnstr_trian_hdiv
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine evaluates constraint coefficients for
!                        the triangular master element and H(div) (L2)
!                        constrained appoximation
!
!   arguments          - none
!
!----------------------------------------------------------------------
!
      subroutine setcnstr_trian_hdiv
!
      use parameters
      use constraints
      implicit none
!
      RRTQ = 0.d0
!
      call setcnstr_trian_iso_hdiv
!
      end subroutine setcnstr_trian_hdiv
