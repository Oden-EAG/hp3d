      subroutine init_cnstr
      use constraints , only : INITIALIZED_CONSTR
      implicit none
c
      if (INITIALIZED_CONSTR) return
c
c  ...H1       
      call setcnstr_edge_h1
      call setcnstr_trian_h1
c
c  ...H(curl) and H(div)
      call setcnstr_edge_hcurl
      call setcnstr_trian_hcurl
      call setcnstr_trian_hdiv
c
      INITIALIZED_CONSTR = .TRUE.
c      
c
      end subroutine init_cnstr
