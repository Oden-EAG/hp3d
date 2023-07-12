!
#include "typedefs.h"
!
!---------------------------------------------------------------------
!> @brief      Returns source term value at a point
!!
!> @param[in]  Mdle     - middle node number
!> @param[in]  X        - a point in physical space
!> @param[out] Fval     - value of source term at the point
!!
!> @date       July 2023
!----------------------------------------------------------------------
   subroutine getf(Mdle,X, Fval)
!      
      use data_structure3D
      use control   , only: NEXACT
      use parameters, only: ZIMG, ZERO
      use common_prob_data_UW, only: OMEGA
!
      implicit none
!
      integer,    intent(in)  :: Mdle
      real(8),    intent(in)  :: X(3)
      complex(8), intent(out) :: Fval(4)
!
      complex(8) :: ValH  (  MAXEQNH    )
      complex(8) :: DvalH (  MAXEQNH,3  )
      complex(8) :: D2valH(  MAXEQNH,3,3)
      complex(8) :: ValE  (3,MAXEQNE    )
      complex(8) :: DvalE (3,MAXEQNE,3  )
      complex(8) :: D2valE(3,MAXEQNE,3,3)
      complex(8) :: ValV  (3,MAXEQNV    )
      complex(8) :: DvalV (3,MAXEQNV,3  )
      complex(8) :: D2valV(3,MAXEQNV,3,3)
      complex(8) :: ValQ  (  MAXEQNQ    )
      complex(8) :: DvalQ (  MAXEQNQ,3  )
      complex(8) :: D2valQ(  MAXEQNQ,3,3)
!
      integer :: icase = 0
!
!----------------------------------------------------------------------
!
      select case(NEXACT)

!  ...unknown exact solution
      case(0)

!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         Fval(1)   = ZIMG*OMEGA*ValQ(1)   + DvalQ(2,1) + DvalQ(3,2) + DvalQ(4,3)
         Fval(2:4) = ZIMG*OMEGA*ValQ(2:4) + DvalQ(1,1:3)
!
!  ...manufactured solution
      case(1)
!
!     ...compute exact solution
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
         Fval(1)   = ZIMG*OMEGA*ValQ(1)   + DvalQ(2,1) + DvalQ(3,2) + DvalQ(4,3)
         Fval(2:4) = ZIMG*OMEGA*ValQ(2:4) + DvalQ(1,1:3)
! 
!  ...known exact solution. Homogeneous RHS
      case(2)
!
         Fval(:) = ZERO
!
      end select
!
! 
   end subroutine getf


!---------------------------------------------------------------------
!> @brief      Returns impedance source term at a point
!!
!> @param[in]  Mdle     - middle node number
!> @param[in]  X        - a point in physical space
!> @param[in]  Rn       - surface normal
!> @param[in]  NBCflag  - boudnary condition flag
!> @param[out] Gval     - value of impedance source term at the point
!!
!> @date       July 2023
!----------------------------------------------------------------------
   subroutine getg(Mdle,X,Rn,NBCflag,Gval)
!      
      use data_structure3D
      use control   , only: NEXACT
      use parameters, only: ZIMG, ZERO
      use common_prob_data_UW
!
      implicit none
!
      integer,     intent(in)  :: Mdle, NBCflag
      real(8),     intent(in)  :: X(3), Rn(3)
      complex(8),  intent(out) :: Gval
!
      complex(8) :: ValH  (  MAXEQNH    )
      complex(8) :: DvalH (  MAXEQNH,3  )
      complex(8) :: D2valH(  MAXEQNH,3,3)
      complex(8) :: ValE  (3,MAXEQNE    )
      complex(8) :: DvalE (3,MAXEQNE,3  )
      complex(8) :: D2valE(3,MAXEQNE,3,3)
      complex(8) :: ValV  (3,MAXEQNV    )
      complex(8) :: DvalV (3,MAXEQNV,3  )
      complex(8) :: D2valV(3,MAXEQNV,3,3)
      complex(8) :: ValQ  (  MAXEQNQ    )
      complex(8) :: DvalQ (  MAXEQNQ,3  )
      complex(8) :: D2valQ(  MAXEQNQ,3,3)
!
      integer :: icase = 0
!
!-----------------------------------------------------------------------
!

      select case(NEXACT)
!  ...unknown exact
      case(0)
 
         select case(PROB_KIND)
!
         case(PROB_CAVITY,PROB_FREESPACE)
!
            select case(NBCflag)
            case(3)
!           ...compute exact solution
               call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                                   ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
               Gval = -valV(1,1)*Rn(1)-valV(2,1)*Rn(2)-valV(3,1)*Rn(3) + valH(1)
!
!           ...for gaussian beam problem, window near corner only
               if (IEXACT_PROB .eq. IEXACT_GAUSS) then
                  Gval = Gval*exp((-x(1)**6-x(2)**6-x(3)**6)*1000)
               endif
!
            end select
!
         case default
            write(*,1000) PROB_KIND
 1000 format('loads: Load not set up for this case of problem. PROB_KIND = ', i1)
            stop 1
         end select
!
!  ...manufactured solution
      case(1)
         select case(NBCflag)
         case(3)
!        ...compute exact solution
            call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                                ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
            Gval = -valV(1,1)*Rn(1)-valV(2,1)*Rn(2)-valV(3,1)*Rn(3) + valH(1)
!         
         end select
!     
!  ...known exact solution. Homogeneous RHS
      case(2)
!      
         select case(NBCflag)
         case(3)
!        ...compute exact solution
            call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                                ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
            Gval = -valV(1,1)*Rn(1)-valV(2,1)*Rn(2)-valV(3,1)*Rn(3) + valH(1)
!         
         end select
!
      end select
!
   end subroutine getg
