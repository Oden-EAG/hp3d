!----------------------------------------------------------------------
!
!   routine name       - setcnstr_trian_hcurl
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine evaluates constraint coefficients for
!                        the triangular master element
!
!   arguments          - none
!
!----------------------------------------------------------------------
!
   subroutine setcnstr_trian_hcurl
!
      use parameters
      use constraints
      implicit none
!
      RRTE = 0.d0
!
      call setcnstr_trian_iso_hcurl
!
   end subroutine setcnstr_trian_hcurl
