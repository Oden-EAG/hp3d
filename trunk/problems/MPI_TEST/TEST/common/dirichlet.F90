!
!----------------------------------------------------------------------
!
!     routine name      - dirichlet
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          - return dirichlet data at a point
!
!     arguments:
!
!     in:
!             X         - a point in physical space
!             Icase     - node case (specifies what variables are supported)
!     out:
!             ValH      - value of the H1 solution
!             DvalH     - corresponding first derivatives
!             DvalE     - value of the H(curl) solution
!             DdvalE    - corresponding first derivatives
!             DvalV     - value of the H(div) solution
!             DdvalV    - corresponding first derivatives
!
!----------------------------------------------------------------------
!
   subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
   use control    , only : NEXACT, GEOM_TOL
   use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
   use common_prob_data
!
   implicit none
!
   real(8), intent(in)  :: X(3)
   integer, intent(in)  :: Icase,Mdle
!..exact solution
   complex(8),dimension(  MAXEQNH    ) ::   ValH
   complex(8),dimension(  MAXEQNH,3  ) ::  DvalH
   complex(8),dimension(  MAXEQNH,3,3) :: d2valH
   complex(8),dimension(3,MAXEQNE    ) ::   ValE
   complex(8),dimension(3,MAXEQNE,3  ) ::  DvalE
   complex(8),dimension(3,MAXEQNE,3,3) :: d2valE
   complex(8),dimension(3,MAXEQNV    ) ::   ValV
   complex(8),dimension(3,MAXEQNV,3  ) ::  DvalV
   complex(8),dimension(3,MAXEQNV,3,3) :: d2valV
   complex(8),dimension(  MAXEQNQ    ) ::   valQ
   complex(8),dimension(  MAXEQNQ,3  ) ::  dvalQ
   complex(8),dimension(  MAXEQNQ,3,3) :: d2valQ
!
!..printing flag
   integer :: iprint
!
!--------------------------------------------------------------------
!
   iprint = 0
!
!..initialize
   ValH = ZERO; DvalH = ZERO
   ValE = ZERO; DvalE = ZERO
   ValV = ZERO; DvalV = ZERO
!
   end subroutine dirichlet
