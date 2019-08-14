!------------------------------------------------------------------------------
!> Purpose : exact (manufactured) solution
!
!> @param[in]  X      - a point in physical space
!> @param[in]  Mdle   - element (middle node) number
!> @param[out] ValH   - value of the H1 solution
!> @param[out] DvalH  - corresponding first derivatives
!> @param[out] D2valH - corresponding second derivatives
!> @param[out] ValE   - value of the H(curl) solution
!> @param[out] DvalE  - corresponding first derivatives
!> @param[out] D2valE - corresponding second derivatives
!> @param[out] ValV   - value of the H(div) solution
!> @param[out] DvalV  - corresponding first derivatives
!> @param[out] D2valV - corresponding second derivatives
!> @param[out] ValQ   - value of the H(div) solution
!> @param[out] DvalQ  - corresponding first derivatives
!> @param[out] D2valQ - corresponding second derivatives
!------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine exact(Xp,Mdle, ValH,DvalH,D2valH, ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, ValQ,DvalQ,D2valQ)
!
  use data_structure3D
  use commonParam
  use laserParam
!
  implicit none
  real*8,dimension(3),            intent(in)  :: Xp
  integer                       , intent(in)  :: Mdle
  VTYPE,dimension(  MAXEQNH    ), intent(out) ::   ValH
  VTYPE,dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
  VTYPE,dimension(  MAXEQNH,3,3), intent(out) :: D2valH
  VTYPE,dimension(3,MAXEQNE    ), intent(out) ::   ValE
  VTYPE,dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
  VTYPE,dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
  VTYPE,dimension(3,MAXEQNV    ), intent(out) ::   ValV
  VTYPE,dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
  VTYPE,dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
  VTYPE,dimension(  MAXEQNQ    ), intent(out) ::   ValQ
  VTYPE,dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
  VTYPE,dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!
!------------------------------------------------------------------------------
!     Space for temporary solutions
!
   VTYPE                    :: E
   VTYPE,dimension(3)       :: dE
   VTYPE,dimension(3,3)     :: d2E
   integer                  :: icomp, fld
!
!..index: signal, pump (E-trace)
   integer, parameter :: s = 1, p = 3
!
!..auxiliary variables
   real*8 :: k,r,n
!
!------------------------------------------------------------------------------
!
!..initialize exact solution
   ValH=ZERO ; DvalH=ZERO ; D2valH=ZERO
   ValE=ZERO ; DvalE=ZERO ; D2valE=ZERO
   ValV=ZERO ; DvalV=ZERO ; D2valV=ZERO
   ValQ=ZERO ; DvalQ=ZERO ; D2valQ=ZERO
!..set default field to signal
   fld = 0
!
!..heat variables
   select case (NO_PROBLEM)
   case(1,2)
      call mfd_solutions(Xp,fld, ValH(1),DvalH(1,1:3),D2valH(1,1:3,1:3))
      ValV(1:3,1) = DvalH(1,1:3)                ! Hdiv
      DvalV(1:3,1,1:3) = D2valH(1,1:3,1:3)      ! 1st derivative
!  ...account for anisotropic short fiber heat operator
      if (ANISO_HEAT .eq. 1) then
         ValV(3,1) = ALPHA_Z*ALPHA_Z*ValV(3,1)
         DvalV(3,1,1:3) = ALPHA_Z*ALPHA_Z*DvalV(3,1,1:3)
      endif
!
   case(3,4)
!..ISOL=12 has radial and azimuthal transverse fields (needs both Ex and Ey)
   if (ISOL .eq. 12) then
      ICOMP_EXACT = 1; icomp = ICOMP_EXACT
!  ...signal E-field trace
      call mfd_solutions(Xp,fld, ValE(icomp,s),DvalE(icomp,s,1:3),D2valE(icomp,s,1:3,1:3))
!  ...y component
      ICOMP_EXACT = 2; icomp = ICOMP_EXACT
!  ...signal E-field trace
      call mfd_solutions(Xp,fld, ValE(icomp,s),DvalE(icomp,s,1:3),D2valE(icomp,s,1:3,1:3))
!  ...signal H-field trace
      !ValE(1,2) = -1.74d0 * ValE(2,s) ! Hx ~ -Ey
      !ValE(2,2) =  1.74d0 * ValE(1,s) ! Hy ~  Ex
      !ValE(3,2) = 0                   ! Hz = 0
!
!..LP modes (polarized in x or y)
   elseif (ISOL .ge. 13 .and. ISOL .le. 15) then
      icomp = ICOMP_EXACT
!  ...signal E-field trace
      fld = 0
      call mfd_solutions(Xp,fld, ValE(icomp,s),DvalE(icomp,s,1:3),D2valE(icomp,s,1:3,1:3))
!  ...pump E-field trace
      fld = 1
      call mfd_solutions(Xp,fld, ValE(icomp,p),DvalE(icomp,p,1:3),D2valE(icomp,p,1:3,1:3))
!
!..LP modes birefringent fiber (using both, polarized in x and y)
   elseif (ISOL .eq. 20 .or. ISOL .eq. 21) then
      ICOMP_EXACT = 1; icomp = ICOMP_EXACT
!  ...signal E-field trace (E_x)
      fld = 0
      call mfd_solutions(Xp,fld, ValE(icomp,s),DvalE(icomp,s,1:3),D2valE(icomp,s,1:3,1:3))
!  ...pump E-field trace (E_x)
      fld = 1
      call mfd_solutions(Xp,fld, ValE(icomp,p),DvalE(icomp,p,1:3),D2valE(icomp,p,1:3,1:3))
!
      ICOMP_EXACT = 2; icomp = ICOMP_EXACT
!  ...signal E-field trace (E_y)
      fld = 0
      call mfd_solutions(Xp,fld, ValE(icomp,s),DvalE(icomp,s,1:3),D2valE(icomp,s,1:3,1:3))
!  ...pump E-field trace (E_y)
      fld = 1
      call mfd_solutions(Xp,fld, ValE(icomp,p),DvalE(icomp,p,1:3),D2valE(icomp,p,1:3,1:3))
!
!..normal procedure (computing one component only)
   else
!  ...set icomp (the only non-zero component of the Etrc solution)
      icomp = ICOMP_EXACT
!
!  ...initialize variables
      E = ZERO; dE = ZERO; d2E = ZERO
!
      call mfd_solutions(Xp,fld, E,dE,d2E)
!
!  ...E-field value
      ValE(icomp,s) = E ! signal E-field trace; (icomp,2) is signal H-field trace
      ValE(icomp,p) = E ! pump   E-field trace; (icomp,4) is pump   H-field trace
!
!  ...E-field 1st order derivatives
      DvalE(icomp,s,1:3) = dE(1:3) ! signal
      DvalE(icomp,p,1:3) = dE(1:3) ! pump
!
!  ...E-field 2nd order derivatives
      D2valE(icomp,s,1:3,1:3) = d2E(1:3,1:3) ! signal
      D2valE(icomp,p,1:3,1:3) = d2E(1:3,1:3) ! pump
!
!  ...Experiment: setting both transverse components non-zero (E_x, E_y)
!      ValE(icomp+1,1:MAXEQNE)           = ValE(icomp,1:MAXEQNE)
!      DvalE(icomp+1,1:MAXEQNE,1:3)      = DvalE(icomp,1:MAXEQNE,1:3)
!      D2valE(icomp+1,1:MAXEQNE,1:3,1:3) = D2valE(icomp,1:MAXEQNE,1:3,1:3)
!
!  ...2nd H(curl) ATTRIBUTE = curl of the first attribute/-i omega \mu
!     signal H-field value (H-field trace)
      ValE(1,2)   = DvalE(3,1,2) - DvalE(2,1,3)
      ValE(2,2)   = DvalE(1,1,3) - DvalE(3,1,1)
      ValE(3,2)   = DvalE(2,1,1) - DvalE(1,1,2)
      ValE(1:3,2) = ValE(1:3,2)/(-ZI*OMEGA*OMEGA_RATIO_SIGNAL*MU)
!
!  ...pump H-field value (H-field trace)
      ValE(1,4)   = DvalE(3,3,2) - DvalE(2,3,3)
      ValE(2,4)   = DvalE(1,3,3) - DvalE(3,3,1)
      ValE(3,4)   = DvalE(2,3,1) - DvalE(1,3,2)
      ValE(1:3,4) = ValE(1:3,4)/(-ZI*OMEGA*OMEGA_RATIO_PUMP*MU)
!
!  ...signal H-field 1st order derivatives
      DvalE(1,2,1) = D2valE(3,1,2,1) - D2valE(2,1,3,1)
      DvalE(1,2,2) = D2valE(3,1,2,2) - D2valE(2,1,3,2)
      DvalE(1,2,3) = D2valE(3,1,2,3) - D2valE(2,1,3,3)
!
      DvalE(2,2,1) = D2valE(1,1,3,1) - D2valE(3,1,1,1)
      DvalE(2,2,2) = D2valE(1,1,3,2) - D2valE(3,1,1,2)
      DvalE(2,2,3) = D2valE(1,1,3,3) - D2valE(3,1,1,3)
!
      DvalE(3,2,1) = D2valE(2,1,1,1) - D2valE(1,1,2,1)
      DvalE(3,2,2) = D2valE(2,1,1,2) - D2valE(1,1,2,2)
      DvalE(3,2,3) = D2valE(2,1,1,3) - D2valE(1,1,2,3)
!
      DvalE(1:3,2,1:3) = DvalE(1:3,2,1:3)/(-ZI*OMEGA*OMEGA_RATIO_SIGNAL*MU)

!  ...pump H-field 1st order derivatives
      DvalE(1,4,1) = D2valE(3,3,2,1) - D2valE(2,3,3,1)
      DvalE(1,4,2) = D2valE(3,3,2,2) - D2valE(2,3,3,2)
      DvalE(1,4,3) = D2valE(3,3,2,3) - D2valE(2,3,3,3)
!
      DvalE(2,4,1) = D2valE(1,3,3,1) - D2valE(3,3,1,1)
      DvalE(2,4,2) = D2valE(1,3,3,2) - D2valE(3,3,1,2)
      DvalE(2,4,3) = D2valE(1,3,3,3) - D2valE(3,3,1,3)
!
      DvalE(3,4,1) = D2valE(2,3,1,1) - D2valE(1,3,2,1)
      DvalE(3,4,2) = D2valE(2,3,1,2) - D2valE(1,3,2,2)
      DvalE(3,4,3) = D2valE(2,3,1,3) - D2valE(1,3,2,3)
!
      DvalE(1:3,4,1:3) = DvalE(1:3,4,1:3)/(-ZI*OMEGA*OMEGA_RATIO_PUMP*MU)
!
!  ...2nd order derivatives (not needed)
!
!  ...L2 components, derivatives not needed
!  ...signal EH fields
      ValQ(1:3) = ValE(1:3,1)
      ValQ(4:6) = ValE(1:3,2)
!  ...pump EH fields
      ValQ(7:9) = ValE(1:3,3)
      ValQ(10:12) = ValE(1:3,4)
!
   endif
   end select
!
end subroutine exact
