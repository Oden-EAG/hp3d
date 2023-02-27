c----------------------------------------------------------------------
c
c   routine name       - setcnstr_trian
c
c----------------------------------------------------------------------
c
c   latest revision    - Jun 10
c
c   purpose            - routine evaluates constraint coefficients for
c                        the triangular master element
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine setcnstr_trian_h1
c
      use parameters
      use constraints
c
      implicit none
      integer :: i
c
      RRTH = 0.d0; RRQH = 0.d0
c
      call setcnstr_trian_iso_h1
      do i=2,4
        call setcnstr_trian_aniso_h1(i)
      enddo
c
      end subroutine setcnstr_trian_h1
