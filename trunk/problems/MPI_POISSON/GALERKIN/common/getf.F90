!----------------------------------------------------------------------
!
!     routine name      - getf
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!     purpose:          - return source term value at a point X
!
!     arguments:
!
!     in:
!             Mdle      - element (middle node) number
!             X         - a point in physical space
!     out:
!             Fval      - Value of source term at the point X
!
!----------------------------------------------------------------------
subroutine getf(Mdle,X, Fval)
!
   use data_structure3D
   use control         , only : NEXACT
   use parameters      , only : ZERO
!
   implicit none
!-------------------------------------------------------------------
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3)
   real(8), intent(out) :: Fval
!-------------------------------------------------------------------
   real(8),dimension(  MAXEQNH    ) ::   ValH
   real(8),dimension(  MAXEQNH,3  ) ::  DvalH
   real(8),dimension(  MAXEQNH,3,3) :: D2valH
   real(8),dimension(3,MAXEQNE    ) ::   ValE
   real(8),dimension(3,MAXEQNE,3  ) ::  DvalE
   real(8),dimension(3,MAXEQNE,3,3) :: D2valE
   real(8),dimension(3,MAXEQNV    ) ::   ValV
   real(8),dimension(3,MAXEQNV,3  ) ::  DvalV
   real(8),dimension(3,MAXEQNV,3,3) :: D2valV
   real(8),dimension(  MAXEQNQ    ) ::   ValQ
   real(8),dimension(  MAXEQNQ,3  ) ::  DvalQ
   real(8),dimension(  MAXEQNQ,3,3) :: D2valQ
!
   integer :: icase
!-------------------------------------------------------------------
!
   Fval = ZERO
   icase = 1
!
!..make sure exact solution is available
   if (NEXACT==0) then
      write(*,*) 'getf: source term cannot be computed;', &
                     '  exact solution is unknown. stop.'
      stop
   endif
!
!..compute exact solution
   call exact(X,icase, ValH,DvalH,d2valH,ValE,DvalE,d2valE, &
                       ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
!..f = -L(u)
   Fval = -(D2valH(1,1,1)+D2valH(1,2,2)+D2valH(1,3,3))
!
end subroutine getf
