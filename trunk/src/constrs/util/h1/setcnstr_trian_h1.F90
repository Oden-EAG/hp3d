!----------------------------------------------------------------------
!
!   routine name       - setcnstr_trian_h1
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine evaluates constraint coefficients for
!                        the triangular master element
!
!   arguments          - none
!
!----------------------------------------------------------------------
!
   subroutine setcnstr_trian_h1
!
      use parameters
      use constraints
!
      implicit none
      integer :: i
!
      RRTH = 0.d0; RRQH = 0.d0
!
      call setcnstr_trian_iso_h1
      do i=2,4
         call setcnstr_trian_aniso_h1(i)
      enddo
!
   end subroutine setcnstr_trian_h1
