subroutine init_cnstr
!
   use constraints, only: INITIALIZED_CONSTR
!
   implicit none
!
   if (INITIALIZED_CONSTR) return
!
!..H1
   call setcnstr_edge_h1
   call setcnstr_trian_h1
!
!..H(curl) and H(div)
   call setcnstr_edge_hcurl
   call setcnstr_trian_hcurl
   call setcnstr_trian_hdiv
!
   INITIALIZED_CONSTR = .TRUE.
!
end subroutine init_cnstr
