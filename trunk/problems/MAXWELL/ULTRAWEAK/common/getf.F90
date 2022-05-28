!----------------------------------------------------------------------
!
!     routine name      - getf
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2022
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
#include "typedefs.h"
subroutine getf(Mdle,X, Fval)
!
   use data_structure3D
   use control         , only : NEXACT
   use parameters      , only : ZERO, ZIMG
   use common_prob_data
!
   implicit none
!-------------------------------------------------------------------
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3)
   complex*16, dimension(6),  intent(out) :: Fval
!-------------------------------------------------------------------
   VTYPE,dimension(  MAXEQNH    ) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: D2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: D2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  DvalV
   VTYPE,dimension(3,MAXEQNV,3,3) :: D2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   ValQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  DvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: D2valQ
!
   complex*16 :: coeff

   integer :: icase
!-------------------------------------------------------------------
!
   Fval = ZERO
   icase = 1
!

   select case(NEXACT)

!..unknown exact solution
   case(0)
!  ...compute exact solution
!
      Fval = ZERO
!..manufactured solution
   case(1)
!  ...compute exact solution
      call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                          ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)

!
!  curl(Ε) + iωμΗ =  0
   coeff = ZIMG*OMEGA*MU
   Fval(1) = DvalE(3,1,2) - DvalE(2,1,3) + coeff*ValE(1,2) 
   Fval(2) = DvalE(1,1,3) - DvalE(3,1,1) + coeff*ValE(2,2) 
   Fval(3) = DvalE(2,1,1) - DvalE(1,1,2) + coeff*ValE(3,2) 
!
!  curl(H) -iωεE = J 
   coeff = ZIMG*OMEGA*EPS
   Fval(4) = DvalE(3,2,2) - DvalE(2,2,3) - coeff*ValE(1,1)  
   Fval(5) = DvalE(1,2,3) - DvalE(3,2,1) - coeff*ValE(2,1)  
   Fval(6) = DvalE(2,2,1) - DvalE(1,2,2) - coeff*ValE(3,1) 
!   
! 
!..known exact solution. Homogeneous RHS
   case(2)
!
       Fval = ZERO
   end select

end subroutine getf


!
!----------------------------------------------------------------------
!                                                                     
!     routine name      - getg
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - May 2022
!                                                                     
!     purpose:          - return boundary source term value at a point
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:   
!             Mdle      - element (middle node) number           
!             X         - a point in physical space on the boundary
!             Rn        - normal vector on the face
!             NBCflag   - boundary condition flag
!     out:              
!             Gval      - Value of the boundary source term at the point
!
!-----------------------------------------------------------------------

subroutine getg(Mdle,X,Rn,NBCflag,Gval)
!      
   use data_structure3D
   use control   , only: NEXACT
   use parameters, only: ZIMG, ZERO
   use common_prob_data
   
   implicit none
!   
!-----------------------------------------------------------------------
!
   integer,     intent(in)  :: Mdle, NBCflag
   real*8,      intent(in)  :: X(3), Rn(3)
   complex*16,  intent(out) :: Gval(3)
   integer :: icase

   VTYPE,dimension(  MAXEQNH    ) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: D2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: D2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  DvalV
   VTYPE,dimension(3,MAXEQNV,3,3) :: D2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   ValQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  DvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: D2valQ

   complex*16 :: rnE(3), rn2E(3), rnH(3)

!
!-----------------------------------------------------------------------

   icase=0
!
   Gval = ZERO
!
   select case(NEXACT)
!
!..unknown exact OR known exact solution but homogeneous RHS
   case(0)
!
   select case(PROB_KIND)
   case(PROB_FREESPACE)
      select case(NBCflag)
      case(9)
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         call zcross_product(Rn,valE(1:3,1),rnE)
         call zcross_product(Rn,rnE, rn2E)
         call zcross_product(Rn,valE(1:3,2),rnH)
!
         Gval = rnH - rn2E
!         
         Gval = Gval*exp((-x(1)**6-x(2)**6-(x(3))**6)*1000)

      end select 

   case(PROB_SCAT_CAVITY)

!..Homogeneous impedance data (absorbing bc)
      Gval = ZERO

   end select   



!..manufactured solution
   case(1,2)

      select case(NBCflag)
!
!  ...impedance BC
      case(9)

         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         call zcross_product(Rn,valE(1:3,1),rnE)
         call zcross_product(Rn,rnE, rn2E)
         call zcross_product(Rn,valE(1:3,2),rnH)
!
         Gval = rnH - rn2E
!
      end select
!
   end select
!
!
end subroutine getg


