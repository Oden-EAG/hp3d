!
#include "typedefs.h"
!
!------------------------------------------------------------------------------
!> @brief      Evaluate source terms at a point
!!
!> @param[in]  Mdle - element (middle node) number
!> @param[in]  X    - physical coordinates
!> @param[out] ZJ   - rhs J^imp
!> @param[out] ZL   - rhs L^fdy
!!
!> @date       October 2024
!------------------------------------------------------------------------------
subroutine getf(Mdle,Xp, ZJ,ZL)
!
   use control     , only: NEXACT
   use parameters  , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
   use commonParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: Xp(3)
   VTYPE  , intent(out) :: ZJ(3),ZL(3)
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
   VTYPE :: zaux,zKH(3),zKE(3)
!
#if HP3D_DEBUG
   integer :: iprint
   iprint = 0
#endif
!
!------------------------------------------------------------------------------
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,7001) Mdle,Xp
 7001 format(' getf: Mdle,Xp = ',i8,2x,3(f8.3,2x))
   endif
#endif
!
!..initialize source terms
   ZJ = ZERO; ZL = ZERO
!
   select case(NEXACT)
!  ............................
!  ...UNKNOWN EXACT SOLUTION...
      case(0)
!  ...They remain as ZERO
!
!  ........................................
!  ...KNOWN EXACT SOLUTION, NON-ZERO RHS...
      case(1)
!     ...compute exact solution
         call exact(Xp,Mdle,valH,dvalH,d2valH,valE,dvalE,d2valE,  &
                            valV,dvalV,d2valV,valQ,dvalQ,d2valQ)
!
!        ...ZJ = curl H - i K.H -iωεE
!        ...first term
            ZJ(1) = dvalE(3,2,2) - dvalE(2,2,3)
            ZJ(2) = dvalE(1,2,3) - dvalE(3,2,1)
            ZJ(3) = dvalE(2,2,1) - dvalE(1,2,2)
!        ...apply the transformation matrix K to H
            call apply_matrixK(Mdle,Xp,valE(1:3,2),zKH)
            zaux = ZI*OMEGA*EPSILON
!        ...add 2nd and 3rd terms            
            ZJ(1:3) = ZJ(1:3) - ZI*zKH - zaux*valE(1:3,1)
!
!        ...ZL = curl E - i K.E +iωμH
!        ...first term
            ZL(1) = dvalE(3,1,2) - dvalE(2,1,3)
            ZL(2) = dvalE(1,1,3) - dvalE(3,1,1)
            ZL(3) = dvalE(2,1,1) - dvalE(1,1,2)
!        ...apply the transformation matrix K to E
            call apply_matrixK(Mdle,Xp,valE(1:3,1),zKE)
            zaux = ZI*OMEGA*MU
!        ...add 2nd and 3rd terms            
            ZL(1:3) = ZL(1:3) - ZI*zKE + zaux*valE(1:3,2)
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
   VTYPE, dimension(3) :: rntimesE,rn2timesE,rntimesH
!
#if HP3D_DEBUG
!..printing flag
   integer :: iprint
   iprint = 0
#endif
!
!------------------------------------------------------------------------------
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
         write(*,7001) Mdle,X
 7001    format(' get_bdSource: Mdle,X = ',i8,2x,3(f8.3,2x))
      endif
#endif
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
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,7002) Imp_val
 7002 format('get_bsource: Imp_val = ',2e12.5)
   endif
#endif
!
end subroutine get_bdSource
!
