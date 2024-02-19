!
#include "typedefs.h"
!
!------------------------------------------------------------------------------
!> @brief      Evaluate source terms at a point
!!
!> @param[in]  Mdle  - element (middle node) number
!> @param[in]  X     - physical coordinates
!> @param[out] Zfval - rhs f
!> @param[out] ZJval - rhs J
!!
!> @date       July 2023
!------------------------------------------------------------------------------
subroutine getf(Mdle,X, ZJval)
!
   use control     , only: NEXACT
   use parameters  , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
   use commonParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3)
   VTYPE  , intent(out) :: ZJval(3)
!
!..exact solution
   VTYPE :: valH  (  MAXEQNH    )
   VTYPE :: dvalH (  MAXEQNH,3  )
   VTYPE :: d2valH(  MAXEQNH,3,3)
   VTYPE :: valE  (3,MAXEQNE    )
   VTYPE :: dvalE (3,MAXEQNE,3  )
   VTYPE :: d2valE(3,MAXEQNE,3,3)
   VTYPE :: valV  (3,MAXEQNV    )
   VTYPE :: dvalV (3,MAXEQNV,3  )
   VTYPE :: d2valV(3,MAXEQNV,3,3)
   VTYPE :: valQ  (  MAXEQNQ    )
   VTYPE :: dvalQ (  MAXEQNQ,3  )
   VTYPE :: d2valQ(  MAXEQNQ,3,3)
!
!..miscellaneus
   integer :: iload,ivar,ibeg,icomp,jcomp,k,l,iprint
   VTYPE   :: zaux
!
!------------------------------------------------------------------------------
!
   iprint = 0
!
   if (iprint.eq.1) then
      write(*,7001) Mdle,X
 7001 format(' getf: Mdle,X = ',i8,2x,3(f8.3,2x))
   endif
!
!..initialize source terms
   ZJval = ZERO
!
   select case(NEXACT)
!  ............................
!  ...UNKNOWN EXACT SOLUTION...
      case(0)
!
         ZJval(1:3) = ZERO
!
!  ........................................
!  ...KNOWN EXACT SOLUTION, NON-ZERO RHS...
      case(1)
!     ...compute exact solution
         call exact(X,Mdle, valH,dvalH,d2valH,valE,dvalE,d2valE,  &
                            valV,dvalV,d2valV,valQ,dvalQ,d2valQ)
!
!        ...RHS = curl H - iωεE
            zaux = ZI*OMEGA*EPSILON
            ZJval(1) = dvalE(3,2,2) - dvalE(2,2,3) - zaux*valE(1,1)
            ZJval(2) = dvalE(1,2,3) - dvalE(3,2,1) - zaux*valE(2,1)
            ZJval(3) = dvalE(2,2,1) - dvalE(1,2,2) - zaux*valE(3,1)
!
!  ...........................................
!  ...KNOWN EXACT SOLUTION, HOMOGENEOUS RHS...
      case(2)
         continue
   end select
!
end subroutine getf





!------------------------------------------------------------------------------
!> @brief      Evaluate impedance load on boundary (from mfd. solution)
!!
!> @param[in]  Mdle     - element (middle node) number
!> @param[in]  X        - physical coordinates
!> @param[in]  Rn       - surface normal
!> @param[out] Imp_val  - impedance value at point
!!
!> @date       July 2023
!------------------------------------------------------------------------------
subroutine get_bdSource(Mdle,X,Rn, Imp_val)
!
   use control          , only : NEXACT
   use parameters       , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
   use commonParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3)
   real(8), intent(in)  :: Rn(3)
   VTYPE  , intent(out) :: Imp_val(3)
!
!------------------------------------------------------------------------------
!
!..exact solution
   VTYPE :: zvalH  (  MAXEQNH    )
   VTYPE :: zdvalH (  MAXEQNH,3  )
   VTYPE :: zd2valH(  MAXEQNH,3,3)
   VTYPE :: zvalE  (3,MAXEQNE    )
   VTYPE :: zdvalE (3,MAXEQNE,3  )
   VTYPE :: zd2valE(3,MAXEQNE,3,3)
   VTYPE :: zvalV  (3,MAXEQNV    )
   VTYPE :: zdvalV (3,MAXEQNV,3  )
   VTYPE :: zd2valV(3,MAXEQNV,3,3)
   VTYPE :: zvalQ  (  MAXEQNQ    )
   VTYPE :: zdvalQ (  MAXEQNQ,3  )
   VTYPE :: zd2valQ(  MAXEQNQ,3,3)
!
!..miscellaneus
   integer :: iload,ivar,ibeg,icomp,jcomp,k,l
   complex(8) :: zaux
   VTYPE, dimension(3) :: rntimesE,rn2timesE,rntimesH
   real(8) :: impedanceConstant
!
!..printing flag
   integer :: iprint = 0
!
!------------------------------------------------------------------------------
!
      if (iprint.eq.1) then
         write(*,7001) Mdle,X
 7001    format(' get_bdSource: Mdle,X = ',i8,2x,3(f8.3,2x))
      endif
!
!     initialize source terms
      Imp_val = ZERO
!
   select case(NEXACT)
      case(0)
         call exact(X,Mdle, zvalH,zdvalH,zd2valH, &
                            zvalE,zdvalE,zd2valE, &
                            zvalV,zdvalV,zd2valV, &
                            zvalQ,zdvalQ,zd2valQ)
!     ...n x E
         call zcross_product(Rn,zvalE(1:3,1), rntimesE)
!     ...n x (n x E)
         call zcross_product(Rn,rntimesE, rn2timesE)
!     ...n x H
         call zcross_product(Rn,zvalE(1:3,2), rntimesH)
!
         Imp_val = rntimesH - rn2timesE
!
!     ...for Gauss beam problem, window the solution near corner
         if (ISOL.eq.5) then
            Imp_val = Imp_val*exp((-x(1)**6-x(2)**6-x(3)**6)*1000)
         endif
!
!  ...exact solution known manufactured solution
      case(1)
         call exact(X,Mdle, zvalH,zdvalH,zd2valH, &
                            zvalE,zdvalE,zd2valE, &
                            zvalV,zdvalV,zd2valV, &
                            zvalQ,zdvalQ,zd2valQ)
!
!     ...n x E
         call zcross_product(Rn,zvalE(1:3,1), rntimesE)
!     ...n x (n x E)
         call zcross_product(Rn,rntimesE, rn2timesE)
!     ...n x H
         call zcross_product(Rn,zvalE(1:3,2), rntimesH)
!
!     ...g = n x H - gamma * n x (n x E)
         Imp_val = rntimesH - GAMMA*rn2timesE
!
      case(2)
         continue
!  ...exact solution unknown
      case default
         write(*,*) 'get_bdSource: UNSPECIFIED NEXACT'
         stop
   end select
   if (iprint.eq.1) then
      write(*,7002) Imp_val
 7002 format('get_bsource: Imp_val = ',2e12.5)
   endif
!
end subroutine get_bdSource
!
