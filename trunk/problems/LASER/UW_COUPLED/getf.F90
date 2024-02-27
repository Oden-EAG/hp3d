!------------------------------------------------------------------------------
!                       subroutine getf
!------------------------------------------------------------------------------
!> Purpose : source term
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  X     - physical coordinates
!! @param[out] Zfval - rhs f
!! @param[out] ZJval - rhs J
!------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine getf(Mdle,X, Zfval,ZJval)
!
   use control     , only: NEXACT
   use parameters  , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
   use commonParam
   use laserParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(in)  :: X(3)
   VTYPE  , intent(out) :: Zfval
   VTYPE  , intent(out) :: ZJval(3)
!
!..auxiliary variable
   VTYPE :: zaux
!
!..exact solution
   VTYPE,dimension(  MAXEQNH    ) ::   valH
   VTYPE,dimension(  MAXEQNH,3  ) ::  dvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   valE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  dvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   valV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  dvalV
   VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   valQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
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
      write(*,7001) Mdle,X
 7001 format(' getf: Mdle,X = ',i8,2x,3(f8.3,2x))
   endif
#endif
!
!..initialize source terms
   Zfval = ZERO
   ZJval = ZERO
!
   select case(NEXACT)
!  ............................
!  ...UNKNOWN EXACT SOLUTION...
      case(0)
         select case(NO_PROBLEM)
!
!        ...single time step of heat
            case(1)
               Zfval = ZERO
               ZJval(1:3) = ZERO
!
!        ...transient heat equation
            case(2)
               Zfval = ZERO
               ZJval(1:3) = ZERO
!
!        ...time-harmonic Maxwell
            case(3,4)
               Zfval = ZERO
               ZJval(1:3) = ZERO
         end select
!
!  ........................................
!  ...KNOWN EXACT SOLUTION, NON-ZERO RHS...
      case(1)
!     ...compute exact solution
         call exact(X,Mdle, valH,dvalH,d2valH,valE,dvalE,d2valE,  &
                            valV,dvalV,d2valV,valQ,dvalQ,d2valQ)
!
         select case(NO_PROBLEM)
!
!        ...single time step
            case(1)
!           ...account for anisotropic short fiber operator yet
               if (ANISO_HEAT .eq. 1) then
                  Zfval = valH(1)-DELTA_T*ALPHA_0 *  &
                          ( d2valH(1,1,1) + &
                            d2valH(1,2,2) + &
                            ALPHA_Z*ALPHA_Z*d2valH(1,3,3) )
               else
                  Zfval = valH(1)-DELTA_T*ALPHA_0 *  &
                          ( d2valH(1,1,1) + &
                            d2valH(1,2,2) + &
                            d2valH(1,3,3) )
               endif
!
!        ...transient heat equation
            case(2)
               write(*,*) 'getf: INCONSISTENCY. stop.'
               stop
!
!        ...time harmonic Maxwell - signal
            case(3)
!           ...RHS = curl H - iωεE
               zaux = ZI*OMEGA*OMEGA_RATIO_SIGNAL*EPSILON + SIGMA
               ZJval(1) = dvalE(3,2,2) - dvalE(2,2,3) - zaux*valE(1,1)
               ZJval(2) = dvalE(1,2,3) - dvalE(3,2,1) - zaux*valE(2,1)
               ZJval(3) = dvalE(2,2,1) - dvalE(1,2,2) - zaux*valE(3,1)
!           ...account for extra term when solving envelope formulation
               if (ENVELOPE) then
!              ...-ik (e_z x H), where e_z x H = (-H_y, H_x, 0)
                  ZJval(1) = ZJval(1) + ZI*WAVENUM_SIGNAL*valE(2,2)
                  ZJval(2) = ZJval(2) - ZI*WAVENUM_SIGNAL*valE(1,2)
               endif
!
!        ...time harmonic Maxwell - pump
            case(4)
!           ...RHS = curl H - iωεE
               zaux = ZI*OMEGA*OMEGA_RATIO_PUMP*EPSILON + SIGMA
               ZJval(1) = dvalE(3,4,2) - dvalE(2,4,3) - zaux*valE(1,3)
               ZJval(2) = dvalE(1,4,3) - dvalE(3,4,1) - zaux*valE(2,3)
               ZJval(3) = dvalE(2,4,1) - dvalE(1,4,2) - zaux*valE(3,3)
!           ...account for extra term when solving envelope formulation
               if (ENVELOPE) then
!              ...-ik (e_z x H), where e_z x H = (-H_y, H_x, 0)
                  ZJval(1) = ZJval(1) + ZI*WAVENUM_PUMP*valE(2,4)
                  ZJval(2) = ZJval(2) - ZI*WAVENUM_PUMP*valE(1,4)
               endif
         end select
!
!  ...........................................
!  ...KNOWN EXACT SOLUTION, HOMOGENEOUS RHS...
      case(2)
!     ...do nothing
   end select
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,7010) Zfval
 7010 format(' getf: Zfval = ',2e12.5)
      call pause
   endif
#endif
!
end subroutine getf
!
!
!------------------------------------------------------------------------------
!                       subroutine get_bdSource
!
! > purpose: computes impedance BC from manufactured solution
!
!------------------------------------------------------------------------------
!
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
   VTYPE,dimension(  MAXEQNH    ) ::   zvalH
   VTYPE,dimension(  MAXEQNH,3  ) ::  zdvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: zd2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   zvalE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  zdvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: zd2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   zvalV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  zdvalV
   VTYPE,dimension(3,MAXEQNV,3,3) :: zd2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   zvalQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  zdvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: zd2valQ
!
   VTYPE,dimension(3) :: rntimesE,rn2timesE,rntimesH
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
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
!  ...known HOMOGENEOUS solution: do nothing
!     impedance constant = 0 (absorbing bc)
      case(0,2)
!
!  ...exact solution known manufactured solution
      case(1)
         call exact(X,Mdle, zvalH,zdvalH,zd2valH, &
                            zvalE,zdvalE,zd2valE, &
                            zvalV,zdvalV,zd2valV, &
                            zvalQ,zdvalQ,zd2valQ)
!     ...exact for signal
         if(NO_PROBLEM.eq.3) then
!        ...n x E
            call zcross_product(Rn,zvalE(1:3,1), rntimesE)
!        ...n x (n x E)
            call zcross_product(Rn,rntimesE, rn2timesE)
!        ...n x H
            call zcross_product(Rn,zvalE(1:3,2), rntimesH)
!     ...exact for pump
         elseif(NO_PROBLEM.eq.4) then
            call zcross_product(Rn,zvalE(1:3,3), rntimesE)
            call zcross_product(Rn,rntimesE, rn2timesE)
            call zcross_product(Rn,zvalE(1:3,4), rntimesH)
         else
            write(*,*) 'error in get_bdsource: NO_PROBLEM must be 3 or 4'
            stop
         endif
!
!     ...g = n x H - gamma * n x (n x E)
         Imp_val = rntimesH - GAMMA*rn2timesE
!
!     ...g should be zero for exact solution (absorbing BC TE10 mode)
         if (abs(sum(Imp_val(1:3))) .ge. 1.0D-14) then
            write(*,*) 'Imp_val is = ', Imp_val
            stop
         endif
!
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
