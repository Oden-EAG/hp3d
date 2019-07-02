c
c----------------------------------------------------------------------
c
c   routine name       - dot_product
c
c----------------------------------------------------------------------
c
c   latest revision    - Apr 14
c
c   purpose            - routine evaluates dot product  of two
c                        vectors in R^3
c
c   arguments :
c     in:    
c           A,B        - vectors in R^3
c     out: 
c           Prod       - value of the product
c
c----------------------------------------------------------------------
c
#include "implicit_none.h"

      subroutine dot_product(A,B, Prod)
c
#include "syscom.blk"
c
      dimension A(3),B(3)
c
      Prod = 0.d0
      do i=1,3
        Prod = Prod + A(i)*B(i)
      enddo
c
      end subroutine dot_product



!                                                                    
!                                                                    
!--------------------------------------------------------------------
!                                                                     
!     routine name      - zdot_product
!                                                                     
!-------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - Feb 2018
!                                                                     
!     purpose:          - compute dot product of real and
!                         (possibly) complex valued vector
!                                                                    
!---------------------------------------------------------------------
!
      subroutine zdot_product(Rn, Za, Zprod)
!
      IMPLICIT NONE

      VTYPE  :: Za(3), Zprod
      real*8 :: Rn(3)
!
!---------------------------------------------------------------------
!
      Zprod = Rn(1)*Za(1)+Rn(2)*Za(2)+Rn(3)*Za(3)
!
!
      end subroutine zdot_product