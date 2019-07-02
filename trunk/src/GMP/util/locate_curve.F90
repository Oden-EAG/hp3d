!----------------------------------------------------------------------
!   routine name       - locate_curve
!----------------------------------------------------------------------
!   latest revision    - Aug 08
!
!   purpose            - Routine locates a curve on a list
!
!   arguments         
!     in:
!            Np1,Np2   - endpoints of a curve
!            List      - list of endpoints for a collection of
!                        curves
!            Nlist     - length of the list
!     out:
!            Ifound    = 0 if the curve is not on the list
!                      = 1 if the curve has been found
!---------------------------------------------------------------------
!
      subroutine locate_curve(Np1,Np2,List,Nlist, Ifound)
!
#include "syscom.blk"     
!
      dimension List(2,Nlist)
!
      Ifound = 0
      do l = 1, Nlist
! ......recall to account for both possible orietations      
        if ((Np1 .eq. List(1,l)) .and. (Np2 .eq. List(2,l)))  Ifound = 1
        if ((Np2 .eq. List(1,l)) .and. (Np1 .eq. List(2,l)))  Ifound = 1
      enddo
!
      end subroutine locate_curve
!
