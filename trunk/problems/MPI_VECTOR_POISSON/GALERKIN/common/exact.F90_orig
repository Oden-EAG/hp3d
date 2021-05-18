!----------------------------------------------------------------------
!
!     routine name      - exact
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!     purpose:          - return exact solution value at a point X
!
!     arguments:
!        in:
!             X         - a point in physical space
!             Icase     - node case (specifies what variables are supported)
!        out:
!             ValH      - value of the H1 solution
!             DvalH     - corresponding first derivatives
!             D2valH    - corresponding second derivatives
!             DvalE     - value of the H(curl) solution
!             DdvalE    - corresponding first derivatives
!             Dd2valE   - corresponding second derivatives
!             DvalV     - value of the H(div) solution
!             DdvalV    - corresponding first derivatives
!             Dd2valV   - corresponding second derivatives
!             DvalQ     - value of the L2 solution
!             DdvalQ    - corresponding first derivatives
!             Dd2valQ   - corresponding second derivatives
!
!----------------------------------------------------------------------
subroutine exact(X,Icase, ValH,DvalH,D2valH, &
                          ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, &
                          ValQ,DvalQ,D2valQ)
!
   use data_structure3D
!
   implicit none
!
!------------------------------------------------------------------------------
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase
!
   real(8),dimension(  MAXEQNH    ), intent(out) ::   ValH
   real(8),dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
   real(8),dimension(  MAXEQNH,3,3), intent(out) :: D2valH
   real(8),dimension(3,MAXEQNE    ), intent(out) ::   ValE
   real(8),dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
   real(8),dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
   real(8),dimension(3,MAXEQNV    ), intent(out) ::   ValV
   real(8),dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
   real(8),dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
   real(8),dimension(  MAXEQNQ    ), intent(out) ::   ValQ
   real(8),dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
   real(8),dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!
   real(8)                   :: u
   real(8), dimension(3)     :: gradu
   real(8), dimension(3,3)   :: gradgradu
!
!------------------------------------------------------------------------------
!
!..initialize exact solution
   ValH = ZERO ; DvalH = ZERO ; D2valH = ZERO
   ValE = ZERO ; DvalE = ZERO ; D2valE = ZERO
   ValV = ZERO ; DvalV = ZERO ; D2valV = ZERO
   ValQ = ZERO ; DvalQ = ZERO ; D2valQ = ZERO
!
!..get exact solution
   call solution(X, u,gradu,gradgradu)
!
   ValH (1)          = u         ! H1
   DvalH(1,1:3)      = gradu     ! 1st derivative [gradient]
   D2valH(1,1:3,1:3) = gradgradu ! 2nd derivative [Hessian]
!
!
end subroutine exact
