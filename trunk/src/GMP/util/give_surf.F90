!----------------------------------------------------------------------
!   routine name       - give_surf
!----------------------------------------------------------------------
!   latest revision    - Dec 08
!
!   purpose            - Given a point and a list surfaces, routine
!                        identifies the surfaces on which the point
!                        is located
!
!   arguments         
!     in:
!            Xp        - point coordinates
!            Nr_confm  - number of surfaces
!            Ns_confm  - list of surfaces
!     out:  
!            Nrs       - number of surfaces to which the point belongs       
!            Nos       - surface numbers
!---------------------------------------------------------------------
!
      subroutine give_surf(Xp,Nr_confm,Ns_confm, Nrs,Nos)
!
      use control
      use GMP
#include "syscom.blk"     
!
      dimension Xp(3),Ns_confm(10),Nos(10)
      dimension dfdx(3)
!
      Nrs = 0
      do is = 1, Nr_confm
        ns = Ns_confm(is)
        call surf(ns,Xp, fval,dfdx)
        if (abs(fval) .lt. GEOM_TOL) then
          Nrs = Nrs + 1
          Nos(Nrs) = ns
        endif
      enddo
!
      end subroutine give_surf
