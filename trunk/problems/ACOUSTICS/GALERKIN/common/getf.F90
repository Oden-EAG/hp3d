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
#include "typedefs.h"
subroutine getf(Mdle,X, Fval)
!
   use data_structure3D
   use control         , only : NEXACT
   use parameters      , only : ZERO
   use common_prob_data, only : OMEGA
!
   implicit none
!-------------------------------------------------------------------
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3)
   complex(8), intent(out) :: Fval
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
!..f = -L(p) - \omega^2 p
   Fval = - D2valH(1,1,1) - d2valH(1,2,2) - d2valH(1,3,3) - OMEGA**2*ValH(1)

end subroutine getf


!
!----------------------------------------------------------------------
!                                                                     
!     routine name      - getg
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
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
   use common_prob_data, only: OMEGA
   
   implicit none
!   
!-----------------------------------------------------------------------
!
   integer,     intent(in)  :: Mdle, NBCflag
   real*8,      intent(in)  :: X(3), Rn(3)
   complex*16,  intent(out) :: Gval
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
!
!-----------------------------------------------------------------------

   icase=0
!
   select case(NEXACT)
!
!..unknown exact
   case(0)
      Gval = ZERO
!
!
!..manufactured solution
   case(1)
      select case(NBCflag)
      case(9)
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                              ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         Gval =  dvalH(1,1)*Rn(1) + dvalH(1,2)*Rn(2) + dvalH(1,3)*Rn(3)  &
               +  ZIMG*OMEGA*valH(1)
!         
      end select  
!     
!..known exact solution. Homogeneous RHS
   case(2)
!      
      select case(NBCflag)
      case(9)
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                              ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         Gval =  dvalH(1,1)*Rn(1) + dvalH(1,2)*Rn(2) + dvalH(1,3)*Rn(3)  &
               +  ZIMG*OMEGA*valH(1)
!         
      end select  
!
   end select
!
!
end subroutine getg