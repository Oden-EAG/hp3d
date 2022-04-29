!----------------------------------------------------------------------
!
!     routine name      - dirichlet
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2022
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
!             ValE      - value of the H(curl) solution
!             DvalE     - corresponding first derivatives
!             ValV      - value of the H(div) solution
!             DvalV     - corresponding first derivatives
!
!----------------------------------------------------------------------
#include "typedefs.h"
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
   use control    , only : NEXACT
   use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
   use common_prob_data
!   
   implicit none
!   
   real*8, dimension(3),          intent(in)  :: X
   integer,                       intent(in)  :: Icase,Mdle
!..exact solution
   complex*16,dimension(  MAXEQNH    ) ::   ValH
   complex*16,dimension(  MAXEQNH,3  ) ::  DvalH
   complex*16,dimension(  MAXEQNH,3,3) :: d2valH
   complex*16,dimension(3,MAXEQNE    ) ::   ValE
   complex*16,dimension(3,MAXEQNE,3  ) ::  DvalE
   complex*16,dimension(3,MAXEQNE,3,3) :: d2valE
   complex*16,dimension(3,MAXEQNV    ) ::   ValV
   complex*16,dimension(3,MAXEQNV,3  ) ::  DvalV
   complex*16,dimension(3,MAXEQNV,3,3) :: d2valV
   complex*16,dimension(  MAXEQNQ    ) ::   valQ
   complex*16,dimension(  MAXEQNQ,3  ) ::  dvalQ
   complex*16,dimension(  MAXEQNQ,3,3) :: d2valQ
!   

   ! real*8, parameter :: c1 = 0.2857142857142857d0 !0.1428571428571428d0
   ! real*8, parameter :: c2 = 0.7142857142857143d0 ! 0.8571428571428571d0
   real*8, parameter :: c1 = 0.1428571428571428d0
   real*8, parameter :: c2 = 0.8571428571428571d0

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
   select case(NEXACT)
!      
!..unknown exact solution
   case(0)
!
      select case(PROB_KIND)
      case(PROB_SCAT_CUBE_PML, PROB_SCAT_SPHR_PML)
!   
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)

!     ...change sign (scattered wave = - incident wave)
!
         ValH = -ValH;  DvalH = -DvalH; d2valH = -d2valH
         ValV = -ValV;  DvalV = -DvalV; d2valV = -d2valV
!            
         if (x(1) .ge. c2) then 
            ValH = ZERO;  DvalH = ZERO; d2valH = ZERO
            ValV = ZERO;  DvalV = ZERO; d2valV = ZERO
         endif
         if (x(1) .lt. c1) then 
            ValH = ZERO;  DvalH = ZERO; d2valH = ZERO
            ValV = ZERO;  DvalV = ZERO; d2valV = ZERO
         endif
         if (x(2) .ge. c2) then 
            ValH = ZERO;  DvalH = ZERO; d2valH = ZERO
            ValV = ZERO;  DvalV = ZERO; d2valV = ZERO
         endif
         if (x(2) .lt. c1) then 
            ValH = ZERO;  DvalH = ZERO; d2valH = ZERO
            ValV = ZERO;  DvalV = ZERO; d2valV = ZERO
         endif            
         if (x(3) .ge. c2) then 
            ValH = ZERO;  DvalH = ZERO; d2valH = ZERO
            ValV = ZERO;  DvalV = ZERO; d2valV = ZERO
         endif
         if (x(3) .lt. c1) then 
            ValH = ZERO;  DvalH = ZERO; d2valH = ZERO
            ValV = ZERO;  DvalV = ZERO; d2valV = ZERO
         endif

      case(PROB_SCAT_CAVITY,PROB_SCAT_SPHERE,PROB_SCAT_PLANE_PML)
!
!     ...scattering problem
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                          ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
!     ...change sign (scattered wave = - incident wave)
!
         ValH = -ValH;  DvalH = -DvalH; d2valH = -d2valH
         ValV = -ValV;  DvalV = -DvalV; d2valV = -d2valV

      end select   
!      
!..known exact solution
   case(1,2)
!  ...use the exact solution to determine Dirichlet data
      call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
!
   case default
      write(*,*)'dirichlet: UNKNOWN EXACT SOLUTION FLAG', NEXACT
      stop 1
   end select

!
!
   end subroutine dirichlet
