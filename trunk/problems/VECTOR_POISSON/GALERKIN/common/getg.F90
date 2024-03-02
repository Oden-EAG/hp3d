!----------------------------------------------------------------------
!
!     routine name      - getg
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2021
!
!> @brief         - return Neumann BC value at a point X
!
!     arguments:
!
!     in:
!             Mdle      - element (middle node) number
!             X         - a point in physical space
!             Rn        - normal to the boundary
!     out:
!             Gval      - Value of source term at the point X
!
!----------------------------------------------------------------------
!
subroutine getg(Mdle,X,Rn, Gval)
!
   use data_structure3D
   use control         , only : NEXACT
   use parameters      , only : ZERO
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3),Rn(3)
   real(8), intent(out) :: Gval
!
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
!
!-------------------------------------------------------------------
!
   Gval = ZERO
   icase = 1
!
   select case(NEXACT)
!
!  ...exact solution unknown
      case(0)
         write(*,*) 'getg: DATA UNDEFINED'
         stop 1
!
!  ...use the exact solution to compute the data
      case(1,2)
!
!     ...compute exact solution
         call exact(X,icase, ValH,DvalH,d2valH,ValE,DvalE,d2valE, &
                             ValV,DvalV,d2valV,valQ,dvalQ,d2valQ)
!
!     ...g = normal derivative
         Gval = DvalH(1,1)*Rn(1) + DvalH(1,2)*Rn(2) + DvalH(1,3)*Rn(3)
!
   end select
!
end subroutine getg
