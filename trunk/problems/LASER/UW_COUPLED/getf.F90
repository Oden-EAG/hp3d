!------------------------------------------------------------------------------
!                       subroutine getf
!------------------------------------------------------------------------------
!> Purpose : source term
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  X     - physical coordinates
!! @param[out] Zfval,zJval - rhs
!------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine getf(Mdle,X, Zfval,zJval)
!
   use control     , only : NEXACT
   use parameters  , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
   use CommonParam
   use LaserParam
!
   implicit none
!
   integer,               intent(in)  :: Mdle
   real*8 , dimension(3), intent(in)  :: X
   VTYPE  ,               intent(out) :: Zfval
   VTYPE  , dimension(3), intent(out) :: ZJval
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
!..miscellaneus
   integer :: iload,ivar,ibeg,icomp,jcomp,k,l,iprint
!
!------------------------------------------------------------------------------
!
!..iprint: 0 - silent
!          1 - verbose
   iprint = 0
!
   if (iprint.eq.1) then
      write(*,7001) Mdle,X
 7001 format(' getf: Mdle,X = ',i8,2x,3(f8.3,2x))
   endif
!
!..initialize source terms
   Zfval = ZERO
   ZJval = ZERO
!
   select case(NEXACT)
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                      |
!==============================================================================
      case(0)
!
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
!        ...time harmonic Maxwell
            case(3,4)
               Zfval = ZERO
               ZJval(1:3) = ZERO
         end select

!
!==============================================================================
!  KNOWN EXACT SOLUTION, NON-ZERO RHS
!==============================================================================
      case(1)
!
!        compute exact solution
         call exact(X,Mdle, valH,dvalH,d2valH, valE,dvalE,d2valE, &
                           valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
!
         select case(NO_PROBLEM)
!
!        ...single time step
            case(1)
               Zfval = -KAPPA*DELTAT*(D2valH(1,1,1) + D2valH(1,2,2) + D2valH(1,3,3)) + valH(1)
!
!        ...transient heat equation
            case(2)
               write(*,*) 'getf: INCONSISTENCY'; stop 1
!
!        ...time harmonic Maxwell - signal
            case(3)
               zaux = ZI*OMEGA*OMEGA_RATIO_SIGNAL*EPSILON + SIGMA
               zJval(1) = DvalE(3,2,2) - DvalE(2,2,3) - zaux*ValE(1,1)
               zJval(2) = DvalE(1,2,3) - DvalE(3,2,1) - zaux*ValE(2,1)
               zJval(3) = DvalE(2,2,1) - DvalE(1,2,2) - zaux*ValE(3,1)
!
!        ...time harmonic Maxwell - pump
            case(4)
               zaux = ZI*OMEGA*OMEGA_RATIO_PUMP*EPSILON + SIGMA
               zJval(1) = DvalE(3,4,2) - DvalE(2,4,3) - zaux*ValE(1,3)
               zJval(2) = DvalE(1,4,3) - DvalE(3,4,1) - zaux*ValE(2,3)
               zJval(3) = DvalE(2,4,1) - DvalE(1,4,2) - zaux*ValE(3,3)
         end select
!
!
!==============================================================================
!  KNOWN EXACT SOLUTION , HOMOGENEOUS RHS                                     |
!==============================================================================
      case(2)
!     ...do nothing
!
   end select
!
   if (iprint.eq.1) then
      write(*,7010) Zfval
 7010 format(' getf: Zfval = ',2e12.5)
      pause
   endif
!
end subroutine getf
!
!
!------------------------------------------------------------------------------
!                       subroutine getf_Newton
!------------------------------------------------------------------------------
!> Purpose : source term
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  X     - physical coordinates
!! @param[out] Zfval,zJval - rhs
!------------------------------------------------------------------------------
!
!
subroutine getf_Newton(Mdle,X, Zfval,zJval)
!
      use control          , only : NEXACT
      use parameters       , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
      use CommonParam
      use LaserParam
!
      implicit none
!
      integer,                  intent(in)  :: Mdle
      real*8,dimension(3),      intent(in)  :: X
      VTYPE,                    intent(out) :: Zfval
      VTYPE,dimension(6),       intent(out) :: ZJval
      VTYPE                                 :: zaux1,zaux2
!------------------------------------------------------------------------------
!
!     exact solution
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
!     miscellaneus
      integer :: iload,ivar,ibeg,icomp,jcomp,k,l
!
!     printing flag
      integer :: iprint
!------------------------------------------------------------------------------
!
!     printing flag : 0 - silent ; 1 - verbose
      iprint=0
!
      if (iprint.eq.1) then
        write(*,7001) Mdle,X
 7001   format(' getf: Mdle,X = ',i8,2x,3(f8.3,2x))
      endif
!
!     initialize source terms
      Zfval = ZERO
      ZJval = ZERO
!
      select case(NEXACT)
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                      |
!==============================================================================
      case(0)
!
        select case(NO_PROBLEM)
!
!  .....single time step of heat
        case(1)
        Zfval = 1.d0
        ZJval(1:6) = ZERO
!
!  .....transient heat equation
        case(2)
        Zfval = 1.d0
        ZJval(1:6) = ZERO

!  .....time harmonic Maxwell
        case(3,4)
        Zfval = 1.d0
        ZJval(1:6) = ZERO
      end select
!
!==============================================================================
!  KNOWN EXACT SOLUTION, NON-ZERO RHS
!==============================================================================
      case(1)
!
!       compute exact solution
        call exact(X,Mdle, valH,dvalH,d2valH, valE,dvalE,d2valE, &
                           valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
!
        select case(NO_PROBLEM)
!
!  .....single time step
        case(1)
          Zfval = -KAPPA*DELTAT*(D2valH(1,1,1) + D2valH(1,2,2) + D2valH(1,3,3)) + valH(1)
!
!  .....transient heat equation
        case(2)
          write(*,*) 'getf: INCONSISTENCY'; stop 1

!  .....time harmonic Maxwell - signal
        case(3)
        zaux1 = ZI*OMEGA*OMEGA_RATIO_SIGNAL*EPSILON + SIGMA
        zaux2 = ZI*OMEGA*OMEGA_RATIO_SIGNAL*MU
!  ..... curl H - (i*omega*epsilon+sigma) E = J(4:6)
        zJval(1) = DvalE(3,2,2) - DvalE(2,2,3) - zaux1*ValE(1,1)
        zJval(2) = DvalE(1,2,3) - DvalE(3,2,1) - zaux1*ValE(2,1)
        zJval(3) = DvalE(2,2,1) - DvalE(1,2,2) - zaux1*ValE(3,1)
!  ..... curl E + (i*omega*mu) H = J(1:3)
        zJval(4) = DvalE(3,1,2) - DvalE(2,1,3) + zaux2*ValE(1,2)
        zJval(5) = DvalE(1,1,3) - DvalE(3,1,1) + zaux2*ValE(2,2)
        zJval(6) = DvalE(2,1,1) - DvalE(1,1,2) + zaux2*ValE(3,2)


!  .....time harmonic Maxwell - pump
        case(4)
        zaux1 = ZI*OMEGA*OMEGA_RATIO_PUMP*EPSILON + SIGMA
        zaux2 = ZI*OMEGA*OMEGA_RATIO_PUMP*MU
!  ..... curl H - (i*omega*epsilon+sigma) E = J(1:3)
        zJval(1) = DvalE(3,4,2) - DvalE(2,4,3) - zaux1*ValE(1,3)
        zJval(2) = DvalE(1,4,3) - DvalE(3,4,1) - zaux1*ValE(2,3)
        zJval(3) = DvalE(2,4,1) - DvalE(1,4,2) - zaux1*ValE(3,3)
!  ..... curl E + (i*omega*mu) H = J(4:6)
        zJval(4) = DvalE(3,3,2) - DvalE(2,3,3) + zaux2*ValE(1,4)
        zJval(5) = DvalE(1,3,3) - DvalE(3,3,1) + zaux2*ValE(2,4)
        zJval(6) = DvalE(2,3,1) - DvalE(1,3,2) + zaux2*ValE(3,4)
        end select
!
!
!==============================================================================
!  KNOWN EXACT SOLUTION , HOMOGENEOUS RHS                                     |
!==============================================================================
      case(2)
!
      endselect
!
      if (iprint.eq.1) then
        write(*,7010) Zfval
 7010   format(' getf_Newton: Zfval = ',2e12.5)
        write(*,*)
        write(*,*) 'getf_Newton: zJval'
        write(*,*)  zJval
        call pause
      endif
!
!
end subroutine getf_Newton
!
!
!------------------------------------------------------------------------------
!                       subroutine get_bdSource
!------------------------------------------------------------------------------
!
subroutine get_bdSource(Mdle,X,Rn, Imp_val)
!
      use control          , only : NEXACT
      use assembly         , only : NR_RHS
      use data_structure3D , only : NR_COMP,ADRES,NRINDEX
      use parameters       , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
      use CommonParam
!
      implicit none
      integer,                  intent(in)  :: Mdle
      real*8,dimension(3),      intent(in)  :: X
      real*8,dimension(3),      intent(in)  :: Rn
      VTYPE,dimension(3),       intent(out) :: Imp_val
!------------------------------------------------------------------------------
!
!     exact solution
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
!     miscellaneus
      integer :: iload,ivar,ibeg,icomp,jcomp,k,l
      complex*16 :: zaux
      VTYPE,dimension(3) ::   rntimesE,rn2timesE
      VTYPE,dimension(3) ::   rntimesH
      real*8                    :: impedanceConstant
      real*8                    :: E   ! vector field
      real*8, dimension(3)      :: dE  ! 1st derivative
      real*8, dimension(3,3)    :: d2E ! 2nd derivative
!
!     printing flag
      integer :: iprint
!------------------------------------------------------------------------------
!
!     printing flag : 0 - silent ; 1 - verbose
      iprint=0
!
      if (iprint.eq.1) then
        write(*,7001) Mdle,X
 7001   format(' get_bdSource: Mdle,X = ',i8,2x,3(f8.3,2x))
      endif
!
!     initialize source terms
      Imp_val = ZERO
!
      select case(NEXACT)
!  ... known HOMOGENEOUS solution: do nothing
      case(0,2)
!
!  ...exact solution known manufactured solution
      case(1)

          call exact(X,Mdle, zvalH,zdvalH,zd2valH, &
                        zvalE,zdvalE,zd2valE, &
                        zvalV,zdvalV,zd2valV, &
                        zvalQ,zdvalQ,zd2valQ)
          !... exact for signal
          if(NO_PROBLEM.eq.3) then
            call zcross_product(Rn,zvalE(1:3,1), rntimesE)
            call zcross_product(Rn,rntimesE, rn2timesE)
            call zcross_product(Rn,zvalE(1:3,2), rntimesH)
          !... exact for pump
          elseif(NO_PROBLEM.eq.4) then
            call zcross_product(Rn,zvalE(1:3,3), rntimesE)
            call zcross_product(Rn,rntimesE, rn2timesE)
            call zcross_product(Rn,zvalE(1:3,4), rntimesH)
          else
            write(*,*) 'error in get_bdsource: NO_PROBLEM must be 3 or 4'
            stop
          endif

!
          Imp_val = rntimesH - ((GAMMA))*rn2timesE
          !write(*,*) 'Imp_val is = ', Imp_val
!
!  ...exact solution unknown
      case default
        write(*,*) 'get_bdSource: UNSPECIFIED NEXACT';stop
      end select
      if (iprint.eq.1) then
        write(*,7002) Imp_val
 7002   format('get_bsource: Imp_val = ',2e12.5)
      endif
!
end subroutine get_bdSource
!
