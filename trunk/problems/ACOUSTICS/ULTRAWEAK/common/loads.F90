!
!----------------------------------------------------------------------
!                                                                     
!     routine name      - getf
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - return source term value at a point
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:   
!             Mdle      - element (middle node) number           
!             X         - a point in physical space
!     out:              
!             Fval      - Value of source term at the point
!
!----------------------------------------------------------------------
!
#include "typedefs.h"

   subroutine getf(Mdle,X, Fval)
!      
   use data_structure3D
   use control   , only: NEXACT
   use parameters, only: ZIMG, ZERO
   use common_prob_data_UW, only: OMEGA
  
   implicit none
!------------------------------------------------------------------------------
   integer,                   intent(in)  :: Mdle
   real*8, dimension(3),      intent(in)  :: X
   complex*16, dimension(4),  intent(out) :: Fval
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

   icase=0

   select case(NEXACT)

!..unknown exact solution
   case(0)
! !  ...compute exact solution
!       call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
!                           ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
! !
!       Fval(1)   = ZIMG*OMEGA*ValQ(1)   + DvalQ(2,1) + DvalQ(3,2) + DvalQ(4,3)
!       Fval(2:4) = ZIMG*OMEGA*ValQ(2:4) + DvalQ(1,1:3)      
!
      Fval = ZERO
!..manufactured solution
   case(1)
!  ...compute exact solution
      call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                          ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
      Fval(1)   = ZIMG*OMEGA*ValQ(1)   + DvalQ(2,1) + DvalQ(3,2) + DvalQ(4,3)
      Fval(2:4) = ZIMG*OMEGA*ValQ(2:4) + DvalQ(1,1:3)
! 
!..known exact solution. Homogeneous RHS
   case(2)
!
       Fval = ZERO
   end select
!
! 
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
   use common_prob_data_UW
  
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
 
   select case(PROB_KIND)
   case(PROB_CAVITY,PROB_FREESPACE,PROB_SPHERE)
      select case(NBCflag)
      case(3)
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         Gval = -valV(1,1)*Rn(1)-valV(2,1)*Rn(2)-valV(3,1)*Rn(3) + valH(1)
!         

         ! Gval = Gval*exp((-x(1)**6-x(2)**6-x(3)**6)/64.d0)
         Gval = Gval*exp((-x(1)**6-x(2)**6-x(3)**6)*1000)

      end select 

   case(PROB_SCAT_CAVITY,PROB_SCAT_SPHERE, PROB_SCAT_PLANE_PML)

!..Homogeneous impedance data (absorbing bc)
   
     Gval = ZERO

   case default
      
      write(*,1000) PROB_KIND
 1000 format('loads: Load not set up for this case of problem. PROB_KIND = ', i1)
   end select   
!
!
!..manufactured solution
   case(1)
      select case(NBCflag)
      case(3)
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         Gval = -valV(1,1)*Rn(1)-valV(2,1)*Rn(2)-valV(3,1)*Rn(3) + valH(1)
!         
      end select  
!     
!..known exact solution. Homogeneous RHS
   case(2)
!      
      select case(NBCflag)
      case(3)
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         Gval = -valV(1,1)*Rn(1)-valV(2,1)*Rn(2)-valV(3,1)*Rn(3) + valH(1)
!         
      end select  
!
   end select
!
!
   end subroutine getg

!----------------------------------------------------------------------
!                                                                     
!     routine name      - getpml
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - Nov 18
!                                                                     
!     purpose:          - return pml derivatives at a physical point
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:   
!             X         - a point in physical space
!     out:              
!             Dfpml     - Value of the derivatives at the point
!
!----------------------------------------------------------------------
!
#include "typedefs.h"

   subroutine getpml(X, Dfpml)
!      
   use data_structure3D
   use parameters, only: ZIMG, ZERO
   use common_prob_data_UW, only: OMEGA, PROB_KIND, PROB_SCAT_CUBE_PML
  
   implicit none
!------------------------------------------------------------------------------
   real*8, dimension(3),      intent(in)  :: X
   complex*16, dimension(3),  intent(out) :: Dfpml
! 
   real*8  :: n 
   real*8  :: lbeg, lend
   real*8, parameter :: c  = 50.0d0
   ! real*8, parameter :: c1 = 0.2857142857142857d0 !0.1428571428571428d0
   ! real*8, parameter :: c2 = 0.7142857142857143d0 ! 0.8571428571428571d0
   real*8, parameter :: c1 = 0.1428571428571428d0
   real*8, parameter :: c2 = 0.8571428571428571d0
   
   ! real*8, parameter :: c1 = 0.25d0
   ! real*8, parameter :: c2 = 0.75d0
   ! real*8, parameter :: c2 = 0.7142857142857143d0

   complex*16 :: zcoeff

!
!-----------------------------------------------------------------------------
!
!
   n = 2.0d0         
   Dfpml = ZONE
   if (PROB_KIND .ne. PROB_SCAT_CUBE_PML) return
!
!
!..stretch in x direction 
   if (X(1) .ge. c2 ) then
      lbeg = c2
      lend = 1.0d0
      zcoeff = ZIMG*(n*c/OMEGA/(lend-lbeg)**int(n)) 
      Dfpml(1) = ZONE + zcoeff * ((x(1)-lbeg)**int((n-1.d0)))
   endif   
   if (X(1) .le. c1 ) then
      lbeg = c1
      lend = 0.0d0
      zcoeff = ZIMG*(n*c/OMEGA/(lend-lbeg)**int(n)) 

      Dfpml(1) = ZONE - zcoeff * ((x(1)-lbeg)**int((n-1.d0)))
   endif

!..stretch in y direction 
   if (X(2) .ge. c2 ) then
      lbeg = c2
      lend = 1.0d0
      zcoeff = ZIMG*(n*c/OMEGA/(lend-lbeg)**int(n)) 
      Dfpml(2) = ZONE + zcoeff * ((x(2)-lbeg)**int((n-1.d0)))
   endif   
   if (X(2) .le. c1 ) then
      lbeg = c1
      lend = 0.0d0
      zcoeff = ZIMG*(n*c/OMEGA/(lend-lbeg)**int(n)) 

      Dfpml(2) = ZONE - zcoeff * ((x(2)-lbeg)**int((n-1.d0)))
   endif

!..stretch in z direction 
   if (X(3) .ge. c2 ) then
      lbeg = c2
      lend = 1.0d0
      zcoeff = ZIMG*(n*c/OMEGA/(lend-lbeg)**int(n)) 
      Dfpml(3) = ZONE + zcoeff * ((x(3)-lbeg)**int((n-1.d0)))
   endif   
   if (X(3) .le. c1 ) then
      lbeg = c1
      lend = 0.0d0
      zcoeff = ZIMG*(n*c/OMEGA/(lend-lbeg)**int(n)) 

      Dfpml(3) = ZONE - zcoeff * ((x(3)-lbeg)**int((n-1.d0)))
   endif
!
! 
   end subroutine getpml
