!
!----------------------------------------------------------------------
!                                                                     
!     routine name      - exact
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - return exact solution value at a point
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:              
!             X         - a point in physical space
!             Icase     - node case (specifies what variables are supported)            
!     out:              
!             ValH      - value of the H1 solution
!             DvalH     - corresponding first derivatives
!             D2valH    - corresponding second derivatives
!             DvalE     - value of the H(curl) solution
!             DdvalE    - corresponding first derivatives
!             Dd2valE   - corresponding second derivatives
!             DvalV     - value of the H(div) solution
!             DdvalV    - corresponding first derivatives
!             Dd2valV   - corresponding second derivatives
!             DvalQ     - value of the H(div) solution
!             DdvalQ    - corresponding first derivatives
!             Dd2valQ   - corresponding second derivatives
!
!----------------------------------------------------------------------
!
#include "implicit_none.h" 
!
   subroutine exact(X,Icase, ValH,DvalH,D2valH, &
                             ValE,DvalE,D2valE, &
                             ValV,DvalV,D2valV, &
                             ValQ,DvalQ,D2valQ)
   use data_structure3D
   use common_prob_data, ONLY: OMEGA
   implicit none
!   
!------------------------------------------------------------------------------
!
   real*8,dimension(3),             intent(in)  :: X
   integer,                         intent(in)  :: Icase
   VTYPE,dimension(  MAXEQNH    ), intent(out) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3), intent(out) :: D2valH
   VTYPE,dimension(3,MAXEQNE    ), intent(out) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
   VTYPE,dimension(3,MAXEQNV    ), intent(out) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
   VTYPE,dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
   VTYPE,dimension(  MAXEQNQ    ), intent(out) ::   ValQ
   VTYPE,dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
   VTYPE,dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!  
   VTYPE                   :: p
   VTYPE, dimension(3)     :: gradp
   VTYPE, dimension(3,3)   :: grad2p
!   
!------------------------------------------------------------------------------
!
!..initialize exact solution
   ValH = ZERO ; DvalH = ZERO ; D2valH = ZERO
   ValE = ZERO ; DvalE = ZERO ; D2valE = ZERO
   ValV = ZERO ; DvalV = ZERO ; D2valV = ZERO
   ValQ = ZERO ; DvalQ = ZERO ; D2valQ = ZERO

!..get solution for the pressure
   call acoustics_solution(X, p, gradp, grad2p)
!

   ValQ(1)           = p                             ! pressure
   DvalQ(1,  1:3)    = gradp                         ! gradient of the pressure
   D2valQ(1,1:3,1:3) = grad2p                        ! Hessian of the pressure
   ValQ(2:4)         = - gradp  / (ZIMG * OMEGA)     ! velocity  
   DvalQ(2:4,1:3)    = - grad2p / (ZIMG * OMEGA)     ! Hessian of velocity

!
   ValH (1)          = p        !  p \hat      (H1)
   DvalH(1,1:3)      = gradp    !  p \hat      1st der
   D2valH(1,1:3,1:3) = grad2p   !  p \hat      2nd der
!
   ValV(1:3,1)       = -gradp  / (ZIMG * OMEGA)    ! u \hat     (Hdiv)
   DvalV(1:3,1,1:3)  = -grad2p / (ZIMG * OMEGA)    ! u \hat    1st der
!
! 
   end subroutine exact