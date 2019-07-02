!------------------------------------------------------------------------------
!> Purpose : exact (manufactured, no miracle!) solution
!!
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
subroutine mfd_solutions(Xp, E,dE,d2E)
!
   use data_structure3D
   use control, only : GEOM_TOL, NEXACT
   use CommonParam
   use LaserParam
!
   implicit none
!  -----------------------------------------------------------------------------------
   real*8, dimension(3),   intent(in)  :: Xp
   VTYPE,                  intent(out) :: E   ! solution
   VTYPE,  dimension(3),   intent(out) :: dE  ! 1st derivative
   VTYPE,  dimension(3,3), intent(out) :: d2E ! 2nd derivative
!  -----------------------------------------------------------------------------------
   real*8 :: x1,x2,x3
   real*8 :: nn        ! for EXPONENTIAL solution
   real*8 :: cn,dn     ! for singular solution
   real*8 :: np_x,np_y,np_z,r0,k0,w0,phase,amplitude,om
    VTYPE  :: f_x,f_y,f_z,df_x,df_y,df_z,ddf_x,ddf_y,ddf_z
!..for bessel function modes
   real*8 :: zeta,rho,xi,theta,rcore,rcladding,alpha,alpha_scale
   real*8 :: order, bessJ, bessY, dbessJ, dbessY,bessJc, bessYc, dbessJc, dbessYc
   real*8 :: bessI, bessK, dbessI, dbessK, bessIc, bessKc, dbessIc, dbessKc
   real*8 :: d2bessJ,d2bessY,d2bessI,d2bessK
   real*8, dimension(2,2) :: hess1,hess2
   real*8 :: angular,angular_x,angular_y,angular_xx,angular_xy,angular_yy
   real*8 :: Jm,Jm_x,Jm_y,Jm_xx,Jm_xy,Jm_yy
   real*8 :: Km,Km_x,Km_y,Km_xx,Km_xy,Km_yy
!
   integer    :: icomp
   VTYPE :: c2z,uz,uz_x,uz_y,uz_z,uz_xx,uz_xy,uz_xz,uz_yy,uz_yx,uz_yz
   VTYPE :: uz_zz,uz_zy,uz_zx
   VTYPE :: pz,pz_x,pz_y,pz_z,pz_xx,pz_xy,pz_xz,pz_yy,pz_yx,pz_yz
   VTYPE :: pz_zz,pz_zy,pz_zx, zbeta,zdbeta,zd2beta
!
!  -----------------------------------------------------------------------------------
!
!..initialize variables
   E = ZERO; dE = ZERO; d2E = ZERO;
   icomp = ICOMP_EXACT
   f_x = ZERO; f_y = ZERO; f_z = ZERO
   df_x = ZERO; df_y = ZERO; df_z = ZERO
   ddf_x = ZERO; ddf_y = ZERO; ddf_z = ZERO
!
!-----------------------------------------------------------------------------------
!      D E C L A R E    S O L U T I O N    V A R I A B L E S                       |
!-----------------------------------------------------------------------------------
!
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
!--------------- 1st prob -------------------------------------------------------
!..a polynomial solution
   if (ISOL .eq. 1) then
!
      np_x=real(NPX,8); np_y=real(NPY,8); np_z=real(NPZ,8)
!
!  ...value
      f_x = x1**np_x
      f_y = x2**np_y
      f_z = x3**np_z
!
!  ...derivatives
      select case(int(np_x))
         case(0); df_x = 0.d0; ddf_x = 0.d0
         case(1); df_x = 1.d0; ddf_x = 0.d0
         case default
            df_x = np_x * x1**(np_x-1.d0)
            ddf_x = np_x * (np_x-1.d0) * x1**(np_x-2.d0)
      end select
!
      select case(int(np_y))
         case(0); df_y = 0.d0; ddf_y = 0.d0
         case(1); df_y = 1.d0; ddf_y = 0.d0
         case default
            df_y = np_y * x2**(np_y-1.d0)
            ddf_y = np_y * (np_y-1.d0) * x2**(np_y-2.d0)
      end select
!
      select case(int(np_z))
         case(0); df_z = 0.d0; ddf_z = 0.d0
         case(1); df_z = 1.d0; ddf_z = 0.d0
         case default
            df_z = np_z * x3**(np_z-1.d0)
            ddf_z = np_z * (np_z-1.d0) * x3**(np_z-2.d0)
      end select
!
!--------------- 0th prob -------------------------------------------------------
!..a polynomial solution with zero boundary values on unit cube
   elseif (ISOL .eq. 0) then
!
      np_x=real(NPX,8); np_y=real(NPY,8); np_z=real(NPZ,8)
!
!  ...value
      f_x = x1**(np_x-1.0d0) * (1.0d0-x1)
      f_y = x2**(np_y-1.0d0) * (1.0d0-x2)
      f_z = x3**(np_z-1.0d0) * (1.0d0-x3)
!
!  ...derivatives
      select case(int(np_x))
         case(0,1,2)
            write(*,*) 'ISOL=0 not yet available for np<=2. stop.'
            stop
         case default
            df_x = (np_x-1.0d0) * x1**(np_x-2.d0) * (1.0d0-x1) &
                 - x1**(np_x-1.d0)
            ddf_x = (np_x-2.0d0) * (np_x-1.0d0) * x1**(np_x-3.d0) * (1.0d0-x1) &
                  - 2.0d0 * (np_x-1.0d0) * x1**(np_x-2.d0)
      end select
!
      select case(int(np_y))
         case(0,1,2)
            write(*,*) 'ISOL=0 not yet available for np<=2. stop.'
            stop
         case default
            df_y = (np_y-1.0d0) * x2**(np_y-2.d0) * (1.0d0-x2) &
                 - x2**(np_y-1.d0)
            ddf_y = (np_y-2.0d0) * (np_y-1.0d0) * x2**(np_y-3.d0) * (1.0d0-x2) &
                  - 2.0d0 * (np_y-1.0d0) * x2**(np_y-2.d0)
      end select
!
      select case(int(np_z))
         case(0,1,2)
            write(*,*) 'ISOL=0 not yet available for np<=2. stop.'
            stop
         case default
            df_z = (np_z-1.0d0) * x3**(np_z-2.d0) * (1.0d0-x3) &
                 - x3**(np_z-1.d0)
            ddf_z = (np_z-2.0d0) * (np_z-1.0d0) * x3**(np_z-3.d0) * (1.0d0-x3) &
                  - 2.0d0 * (np_z-1.0d0) * x3**(np_z-2.d0)
      end select
!
!
!----------- 2nd prob -------------------------------------------------------
!..a smooth solution
   elseif (ISOL .eq. 2) then
!
      f_x= sin(OMEGA*Xp(1))
      f_y= sin(OMEGA*Xp(2))
      f_z= sin(OMEGA*Xp(3))
!
!     1st order derivatives
      df_x=(OMEGA)*cos(OMEGA*Xp(1))
      df_y=(OMEGA)*cos(OMEGA*Xp(2))
      df_z=(OMEGA)*cos(OMEGA*Xp(3))
!
!     2nd order derivatives
      ddf_x=-OMEGA**2*f_x
      ddf_y=-OMEGA**2*f_y
      ddf_z=-OMEGA**2*f_z
!
!----------- 3rd prob -------------------------------------------------------
!..a plane wave
   elseif (ISOL .eq. 3) then
!
      f_x= Xp(1)
      f_y= 1.d0
      f_z= 2.d0*cdexp(-ZI*OMEGA*Xp(3))
!
!     1st order derivatives
      df_x= 1.d0
      df_y= 0.d0
      df_z=(-ZI*OMEGA)*f_z
!
!     2nd order derivatives
      ddf_x=0.d0
      ddf_y=0.d0
      ddf_z=-OMEGA**2*f_z
!
!----------- 4th prob -------------------------------------------------------
!..an exponential
   elseif (ISOL .eq. 4) then
!
      nn = -OMEGA
!  ...value
      f_x=exp(-Xp(1)**2/2)
      f_y=exp(-Xp(2)**2/2)
      f_z=1.d0!exp(nn*Xp(3))
!
!  ...1st order derivatives
      df_x=-Xp(1)*f_x
      df_y=-Xp(2)*f_y
      df_z=0.d0!nn*f_z
!
!  ...2nd order derivatives
      ddf_x=-f_x-(Xp(1)*df_x)
      ddf_y=-f_y-(Xp(2)*df_y)
      ddf_z=0.d0!nn*df_z

!----------- 5th prob -------------------------------------------------------
!..fundamental TE10 mode for rectangular waveguide
   elseif (ISOL .eq. 5) then

      f_x=-ZI*(OMEGA/PI)*sin(PI*Xp(1))
      f_y= 1.d0
      f_z=cdexp(-ZI*OMEGA*Xp(3)*GAMMA)
!
!     1st order derivatives
      df_x=-ZI*(OMEGA/PI)*PI*cos(PI*Xp(1))
      df_y=0.d0
      df_z=(-ZI*OMEGA*GAMMA)*f_z
!
!     2nd order derivatives
      ddf_x=-PI**2*f_x
      ddf_y=0.d0
      ddf_z=(-ZI*OMEGA*GAMMA)*df_z

!
!
!----------- 6th prob -------------------------------------------------------
!..Gaussian pulse plane wave
   elseif (ISOL .eq. 6) then
      if(NEXACT.ne.0) then
         write(*,*) 'NEXACT must be 0 for ISOL=6'
         stop
      endif

      w0 = BEAM_WAIST
      k0 = OMEGA
      c2z = w0**2
      x1 = Xp(1)
      x2 = Xp(2)
      r0 = dsqrt((x1)**2+(x2)**2)
      uz = cdexp(-ZI*k0*Xp(3))*cdexp(-r0**2/c2z)/c2z
!
      uz_x = -2.d0*x1*uz/c2z
      uz_y = -2.d0*x2*uz/c2z
      uz_z = uz*(-ZI*k0)
!
      uz_xx = -2.d0*(uz+x1*uz_x)/c2z
      uz_xy = -2.d0*x1*uz_y/c2z
      uz_xz = -2.d0*x1*uz_z/c2z
!
      uz_yy = -2.d0*(uz+x2*uz_y)/c2z
      uz_yx = uz_xy
      uz_yz = -2.d0*x2*uz_z/c2z
!
      uz_zz = uz_z*(-ZI*k0)
      uz_zx = uz_xz
      uz_zy = uz_yz
!
      E=uz
!  ...1st order derivatives
      dE(1) = uz_x
      dE(2) = uz_y
      dE(3) = uz_z
!
!  ...2nd order derivatives
      d2E(1,1) = uz_xx
      d2E(1,2) = uz_xy
      d2E(1,3) = uz_xz
      d2E(2,1) = d2E(1,2)
      d2E(2,2) = uz_yy
      d2E(2,3) = uz_yz
      d2E(3,1) = d2E(1,3)
      d2E(3,2) = d2E(2,3)
      d2E(3,3) = uz_zz
!
!----------- 7th prob -------------------------------------------------------
!...... 3D Plane Wave + Gaussian beam
   elseif (ISOL .eq. 7) then
!
      if(NEXACT.ne.0) then
         write(*,*) 'NEXACT must be 0 for ISOL=7'
         stop
      endif
!.... beam parameters
          w0 = (0.65d0+1.619d0*(2.d0*PI*VNUM)**(-1.5d0)+2.879d0*(2.d0*PI*VNUM)**(-6.d0))*R_CORE
          !w0 = 0.1d0
          alpha_scale = 1.d0
          k0 = OMEGA
          c2z = w0**2
          c2z= c2z/alpha_scale
          if(NO_PROBLEM.eq.3) then
            if(COPUMP.eq.1) then
              ! copumped pumped, amplitude = 2.0d0
              amplitude = 2.0d0
            elseif(COPUMP.eq.0) then
              ! counter pumped, amplitude = 1.5d0
              amplitude = 1.50d0
            else
              write(*,*) 'error in ISOL = 7 mfd_solution: COPUMP = 0 or 1 only'
              stop
            endif
          elseif(NO_PROBLEM.eq.4) then
            amplitude = 2.25d0
          endif
          x1 = Xp(1)
          x2 = Xp(2)
          r0 = dsqrt((x1)**2+(x2)**2)
!
          if(NO_PROBLEM.eq.3) then
            uz = amplitude*cdexp(-ZI*k0*OMEGA_RATIO_SIGNAL*Xp(3))*cdexp(-r0**2/c2z)/c2z
          elseif(NO_PROBLEM.eq.4) then
            uz = amplitude*cdexp(-ZI*k0*OMEGA_RATIO_PUMP*Xp(3))*cdexp(-r0**2/c2z)/c2z
          endif
    !
          uz_x = -2.d0*x1*uz/c2z
          uz_y = -2.d0*x2*uz/c2z
          uz_z = uz*(-ZI*k0)
    !
          uz_xx = -2.d0*(uz+x1*uz_x)/c2z
          uz_xy = -2.d0*x1*uz_y/c2z
          uz_xz = -2.d0*x1*uz_z/c2z

          uz_yy = -2.d0*(uz+x2*uz_y)/c2z
          uz_yx = uz_xy
          uz_yz = -2.d0*x2*uz_z/c2z

          uz_zz = uz_z*(-ZI*k0)
          uz_zx = uz_xz
          uz_zy = uz_yz

!.....     If in the core, exact is  uz
!.....     else, only ZERO in the cladding
          if(dsqrt(Xp(1)**2+Xp(2)**2).le.(R_CORE)) then
            E= uz
!  .....1st order derivatives
            dE(1) = uz_x
            dE(2) = uz_y
            dE(3) = uz_z
    !
!  .....2nd order derivatives
            d2E(1,1) = uz_xx
            d2E(1,2) = uz_xy
            d2E(1,3) = uz_xz
            d2E(2,1) = d2E(1,2)
            d2E(2,2) = uz_yy
            d2E(2,3) = uz_yz
            d2E(3,1) = d2E(1,3)
            d2E(3,2) = d2E(2,3)
            d2E(3,3) = uz_zz
!
          else
            E= ZERO
   !  .....1st order derivatives
            dE(1) = ZERO
            dE(2) = ZERO
            dE(3) = ZERO
  !  !
  !  .....2nd order derivatives
            d2E(1,1) = ZERO
            d2E(1,2) = ZERO
            d2E(1,3) = ZERO
            d2E(2,1) = ZERO
            d2E(2,2) = ZERO
            d2E(2,3) = ZERO
            d2E(3,1) = ZERO
            d2E(3,2) = ZERO
            d2E(3,3) = ZERO
!!.....  end if for core/cladding check
        endif

!----------- 8th prob -------------------------------------------------------
!..Bessel Function Modes
   elseif (ISOL .eq. 8) then
        if(NEXACT.ne.0) then
          write(*,*) 'NEXACT must be 0 for ISOL=8'
          stop
        endif
          if(NO_PROBLEM.eq.3) then
            amplitude = 55.d0
            !amplitude = 1000.d0
          elseif(NO_PROBLEM.eq.4) then
            amplitude = 65.d0
          endif
          x1 = Xp(1)+GEOM_TOL
          x2 = Xp(2)+GEOM_TOL
          x3 = Xp(3)
          rho = dsqrt((x1)**2+(x2)**2)
          theta = datan(x2/x1)
          order = ORDER_BESSEL
          rcore = R_CORE
          rcladding = R_CLAD
          xi = 1.d0
          zeta = 1.d0
          zeta = (2.d0*PI*(12.d0/R_CORE))*dsqrt(REF_INDEX_CORE**2-BETA_PROPAGATION**2)
          alpha = 1.d0
          alpha = (2.d0*PI*(12.d0/R_CORE))*dsqrt(BETA_PROPAGATION**2-REF_INDEX_CLAD**2)

          hess1(1,1) = 1.d0/rho - (1.d0/rho**3)*x1**2
          hess1(1,2) = -x1*x2/rho**3
          hess1(2,1) = hess1(1,2)
          hess1(2,2) = 1.d0/rho - (1.d0/rho**3)*x2**2


          hess2(1,1) = (-order**2*x2**2*cos(order*theta)/rho**4)+ &
                       (2.d0*order*x2**3*sin(order*theta)/(x1*rho**4)) - &
                       (2.d0*order*x2*sin(order*theta)/(x1*rho**2))
          hess2(1,2) = (order**2*x2*x1*cos(order*theta)/rho**4) + &
                       (order*sin(order*theta)/rho**2) - &
                       (2.d0*order*x2**2*sin(order*theta)/(rho**4))
          hess2(2,1) = hess2(1,2)
          hess2(2,2) = (2.d0*order*x2*x1*sin(order*theta)/rho**4)- &
                       (order**2*x1**2*cos(order*theta)/rho**4)

          angular = cos(order*theta)
          angular_x = order*x2*sin(order*theta)/rho**2
          angular_y = -order*x1*sin(order*theta)/rho**2
          angular_xx = hess2(1,1)
          angular_xy = hess2(1,2)
          angular_yy = hess2(2,2)
!       Check if in core
          if(rho.le.rcore) then
            call dbessJY(rho*zeta, order, bessJ, bessY, dbessJ, dbessY)
            call dbessJY(rcore*zeta, order, bessJc, bessYc, dbessJc, dbessYc)
            call d2bessJY(rho*zeta, order, d2bessJ, d2bessY)
            Jm = bessJ
            Jm_x = dbessJ*zeta*x1/rho
            Jm_y = dbessJ*zeta*x2/rho
            Jm_xx = d2bessJ*(zeta*x1/rho)**2 + dbessJ*zeta*hess1(1,1)
            Jm_xy = d2bessJ*(zeta)**2*(x1/rho)*(x2/rho) + dbessJ*zeta*hess1(1,2)
            Jm_yy = d2bessJ*(zeta*x2/rho)**2 + dbessJ*zeta*hess1(2,2)
!
            uz = amplitude*(xi/bessJc)*(bessJ)*angular*cdexp(-ZI*BETA_PROPAGATION*x3)
            uz_x = (xi/bessJc)*(Jm_x*angular+Jm*angular_x)
            uz_y = (xi/bessJc)*(Jm_y*angular+Jm*angular_y)
            uz_xx = (xi/bessJc)*(Jm_xx*angular+2.d0*Jm_x*angular_x+Jm*angular_xx)
            uz_xy = (xi/bessJc)*(Jm_xy*angular+Jm_x*angular_y+Jm_y*angular_x+Jm*angular_xy)
            uz_yy = (xi/bessJc)*(Jm_yy*angular+2.d0*Jm_y*angular_y+Jm*angular_yy)
          else
            call dbessIK(rho*alpha, order, bessI, bessK, dbessI, dbessK)
            call dbessIK(rcore*alpha, order, bessIc, bessKc, dbessIc, dbessKc)
            call d2bessIK(rho*alpha, order, d2bessI, d2bessK)

            Km = bessK
            Km_x = dbessK*alpha*x1/rho
            Km_y = dbessK*alpha*x2/rho
            Km_xx = d2bessK*(alpha*x1/rho)**2 + dbessK*alpha*hess1(1,1)
            Km_xy = d2bessK*(alpha)**2*(x1/rho)*(x2/rho) + dbessK*alpha*hess1(1,2)
            Km_yy = d2bessK*(alpha*x2/rho)**2 + dbessK*alpha*hess1(2,2)
!
            uz = amplitude*(xi/bessKc)*(bessK)*angular*cdexp(-ZI*BETA_PROPAGATION*x3)
            uz_x = (xi/bessKc)*(Km_x*angular+Km*angular_x)
            uz_y = (xi/bessKc)*(Km_y*angular+Km*angular_y)
            uz_xx = (xi/bessKc)*(Km_xx*angular+2.d0*Km_x*angular_x+Km*angular_xx)
            uz_xy = (xi/bessKc)*(Km_xy*angular+Km_x*angular_y+Km_y*angular_x+Km*angular_xy)
            uz_yy = (xi/bessKc)*(Km_yy*angular+2.d0*Km_y*angular_y+Km*angular_yy)

          endif    !
          uz_z = uz*(-ZI*BETA_PROPAGATION)
          uz_xz = uz_x*(-ZI*BETA_PROPAGATION)
          uz_yz = uz_y*(-ZI*BETA_PROPAGATION)
          uz_zz = uz*(-BETA_PROPAGATION**2)
    !
          E=uz
    !  .....1st order derivatives
          dE(1) = uz_x
          dE(2) = uz_y
          dE(3) = uz_z
!
!  .....2nd order derivatives
          d2E(1,1) = uz_xx
          d2E(1,2) = uz_xy
          d2E(1,3) = uz_xz
          d2E(2,1) = d2E(1,2)
          d2E(2,2) = uz_yy
          d2E(2,3) = uz_yz
          d2E(3,1) = d2E(1,3)
          d2E(3,2) = d2E(2,3)
          d2E(3,3) = uz_zz

!----------- 9th prob -------------------------------------------------------
!  ...fundamental TE10 mode for rectangular waveguide with exponential
!  ... growth
   elseif (ISOL .eq. 9) then
        f_x=2.d0*ZI*sin(PI*Xp(1))
        f_y= 1.d0
        f_z=cdexp((EXP_COEFF-ZI*OMEGA)*Xp(3))
  !
  !     1st order derivatives
        df_x=PI*2.d0*ZI*cos(PI*Xp(1))
        df_y=0.d0
        df_z=(EXP_COEFF-ZI*OMEGA)*f_z
  !
  !     2nd order derivatives
        ddf_x=-PI**2*f_x
        ddf_y=0.d0
        ddf_z=(EXP_COEFF-ZI*OMEGA)*df_z

!----------- 10th prob -------------------------------------------------------
!  ...fundamental TE10 mode for rectangular waveguide with exponential
!  ... growth with complex stretching
   elseif (ISOL .eq. 10) then
        call get_Beta(Xp, zbeta,zdbeta,zd2beta)
        f_x=2.d0*ZI*sin(PI*Xp(1))
        f_y= 1.d0
        f_z=cdexp((EXP_COEFF-ZI*OMEGA)*zbeta)
  !
  !     1st order derivatives
        df_x=PI*2.d0*ZI*cos(PI*Xp(1))
        df_y=0.d0
        df_z=(EXP_COEFF-ZI*OMEGA)*zdbeta*f_z
  !
  !     2nd order derivatives
        ddf_x=-PI**2*f_x
        ddf_y=0.d0
        ddf_z=(EXP_COEFF-ZI*OMEGA)*zd2beta*f_z + &
              (EXP_COEFF-ZI*OMEGA)*zdbeta*df_z
!----------- 11th prob -------------------------------------------------------
!  ...f(x,y,z) = sin(OMEGA*x)sin(OMEGA*y)
   elseif (ISOL .eq. 11) then

      f_x= sin(OMEGA*Xp(1))
      f_y= sin(OMEGA*Xp(2))
      f_z= 1.d0
!
!     1st order derivatives
      df_x=OMEGA*cos(Xp(1))
      df_y=OMEGA*cos(Xp(2))
      df_z=0.d0
!
!     2nd order derivatives
      ddf_x=-OMEGA**2*f_x
      ddf_y=-OMEGA**2*f_y
      ddf_z=0.d0
!
!..endif ISOL
   endif
!
!
!.... set values as tensor products for all cases except
!.... Gaussian pulse and fiber beam
   select case(ISOL)
      case(0,1,2,3,4,5,9,10,11)
!     ...value
         E=f_x*f_y*f_z
!     ...1st order derivatives
         dE(1)=  df_x *   f_y *   f_z
         dE(2)=   f_x *  df_y *   f_z
         dE(3)=   f_x *   f_y *  df_z
!
!     ...2nd order derivatives
         d2E(1,1) = ddf_x *   f_y *   f_z
         d2E(1,2) =  df_x *  df_y *   f_z
         d2E(1,3) =  df_x *   f_y *  df_z
         d2E(2,1) =  d2E(1,2)
         d2E(2,2) =   f_x * ddf_y *   f_z
         d2E(2,3) =   f_x *  df_y *  df_z
         d2E(3,1) =  d2E(1,3)
         d2E(3,2) =  d2E(2,3)
         d2E(3,3) =   f_x *   f_y * ddf_z
!
      case default
         write(*,*) 'invalid ISOL param. stop.'
         stop
   end select
end subroutine mfd_solutions



!!!!!!!!!!!!!!!!!!ROUTINES FOR BESSEL FUNCTIONS FROM JAKE!!!!!!!!!!!!!!!!!
!
! File:          BesselFunctions.f95
!
! Author:        Jacob Grosek
! Start Date:    March 31, 2017
! Last Modified: March 31, 2017
! Air Force Research Laboratory at Kirtland AFB, Albuquerque, NM
! Directed Energy Directorate, Laser Division, Modeling & Simulation Program
!
! Description of the File:
!
! This file contains Fortran Bessel function algorithms, and their
! derivatives.  These Bessel functions only accept real values, and they
! only output real values as well.
!
! This code is modified from "Numerical Recipes in Fortran 77: The Art
! of Scientific Computing," Cambridge University Press (1992)
!
! File Parameters/Variables:
!
! x
!   - (input, real, [ ]) the given point at which the Bessel
!     functions will be calculated
! order
!   - (input, real, [ ]) the order of the Bessel functions
! bessJ
!   - (output, real, [ ]) the value of the Bessel J function
! bessY
!   - (output, real, [ ]) the value of the Bessel Y function
! bessJder
!   - (output, real, [ ]) the value of the derivative of the
!     Bessel J function
! bessYder
!   - (output, real, [ ]) the value of the derivative of the
!     Bessel Y function
!
! File Dependencies:
!
! This subroutine calls upon the subroutine beschbgamma defined in
! this same file.  This subroutine returns the Bessel functions
! bessJ = J_order and bessY = Y_order, and their respective
! derivatives bessJder = J_order′ and bessYder = Y_order′, for
! positive x and of an order greater or equal to 0.  The relative
! accuracy is within one or two significant digits of eps, except
! near a zero of one of the functions, where eps controls its
! absolute accuracy. The parameter "floatingpointmin" is a number
! close to the machine’s smallest floating-point number. All
! internal arithmetic is accomplished with real numbers. In oder to
! convert the entire routine to double precision, use the real(8)
! declaration and decrease eps to 10**−16.  Also convert the
! subroutine beschbgamma.
!

subroutine dbessJY(x, order, bessJ, bessY, bessJder, bessYder)

  implicit none

  integer, parameter :: maxit = 10000
  real(8) :: bessJ,bessJder, bessY, bessYder, x, order
  real(8), parameter :: eps = 1.d-16, floatingpointmin = 1.e-30, &
    pi = 3.141592653589793d0, xmin = 2.d0, &
    emc = 0.577215664901533d0
  integer :: i, signi, l, nl
  real(8) :: a, b, br, bi, c, cr, ci, d, del, del1, den, &
    di, dlr, dli, dr, e, f, fact, fact2, fact3, ff, gam, &
    gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, &
    bessJl, bessJl1, bessJmu, bessJder1, bessJderl, &
    bessJtemp, bessY1, bessYmu, bessYmup, bessYtemp, summ, &
    summ1, temp, w, x2, xi, xi2, xmu, xmu2, gampl1, negate

  if (order < 0.d0) then
    order  = dabs(order)
    negate = (-1.d0)**int(order)
  else
    negate = 1.d0
  endif ! end if statement
  if (dabs(x) <= 1.d-6) then
    if ((dabs(x) <= 1.d-17) .and. (order == 0.d0)) then
      bessJ    = 1.d0 *  negate
      bessJder = 0.d0
      bessY    = -1.d30 * negate
      bessYder = 1.d30 * negate
    elseif ((dabs(x) <= 1.d-17) .and. (order /= 0.d0)) then
      bessJ    = 0.d0
      bessJder = negate * 2.d0 * order * &
        (x / 2.d0)**(order - 1.d0) / gamma(order + 1.d0)
      bessY    = -1.d30 * negate
      bessYder = 1.d30 * negate
    elseif ((dabs(x) > 1.d-17) .and. (order == 0.d0)) then
      bessJder = negate * 2.d0 * order * &
        (x / 2.d0)**(order - 1.d0) / gamma(order + 1.d0)
      bessY    = negate * (2.d0 / pi) * (dlog(x / 2.d0) + emc)
      bessYder = negate * 2.d0 / (pi * x)
    elseif ((dabs(x) > 1.d-17) .and. (order /= 0.d0)) then
      bessJ    = negate * (x / 2.d0)**order / &
        gamma(order + 1.d0)
      bessJder = negate * 2.d0 * order * &
        (x / 2.d0)**(order - 1.d0) / gamma(order + 1.d0)
      bessY    = -negate * (gamma(order) / pi) * &
        (2.d0 / x)**order
      bessYder = negate * (order * gamma(order) / (pi * x)) * &
        (2.d0 / x)**order
    endif ! end if statement
    return
  endif ! end if statement

  ! Here "nl" is the number of downward recurrences of the J’s and
  ! upward recurrences of Y’s. "xmu" lies between −1/2 and 1/2 for
  ! x < xmin, while it is chosen so that x is greater than the
  ! turning point for x >= xmin:
  if (x < xmin) then
    nl = int(order + 0.5d0)
  else
    nl = max(0, int(order - x + 1.5d0))
  endif ! end if statement
  xmu  = order - nl
  xmu2 = xmu * xmu
  if (x /= 0.d0) then
    xi = 1.d0 / x
  else
    xi = 1.d30
  endif ! end if statement
  xi2  = 2.d0 * xi
  w    = xi2 / pi
    ! the Wronskian

  ! Here the first continued fraction is found by the modified
  ! Lentz's method.  The variable "signi" keeps track of sign
  ! changes in the denominator:
  signi = 1
  h     = order * xi
  if (h < floatingpointmin) then
    h = floatingpointmin
  endif ! end if statement
  b = xi2 * order
  d = 0.d0
  c = h
  do i = 1, maxit
    b = b + xi2
    d = b - d
    if (abs(d) < floatingpointmin) then
      d = floatingpointmin
    endif ! end if statement
    c = b - 1.d0 / c
    if (abs(c) < floatingpointmin) then
      c = floatingpointmin
    endif ! end if statement
    d   = 1.d0 / d
    del = c * d
    h   = del * h
    if (d < 0.d0) then
      signi = -signi
    endif ! end if statement
    if (dabs(del - 1.d0) < eps) then
      goto 1
    endif ! end if statement
  enddo
  print *, "x is too large in dbessJY; try asymptotic expansion"
  print *, "x = ", x
  1 continue

  ! Here J_order and J_order' are initialized for downward
  ! recurrence:
  bessJl    = signi * floatingpointmin
  bessJderl = h * bessJl
  bessJl1   = bessJl
    ! store values for future rescaling
  bessJder1 = bessJderl
  fact      = order * xi
  do l = nl, 1, -1
    bessJtemp = fact * bessJl + bessJderl
    fact      = fact - xi
    bessJderl = fact * bessJtemp - bessJl
    bessJl    = bessJtemp
  enddo
  if (bessJl == 0.d0) then
    bessJl = eps
  endif ! end if statement

  ! The subroutine already has calculated unnormalized J_order and
  ! J_order':
  f = bessJderl / bessJl
  if (x < xmin) then ! use series
    x2   = 0.5d0 * x
    pimu = pi * xmu
    if (abs(pimu) < eps) then
      fact = 1.d0
    else
      fact = pimu / dsin(pimu)
    endif ! end if statement
    d = -dlog(x2)
    e = xmu * d
    if (abs(e) < eps) then
      fact2 = 1.d0
    else
      fact2 = dsinh(e) / e
    endif ! end if statement
    call beschbgamma(xmu, gam1, gam2, gampl, gammi)
      ! Chebyshev evaluation of Gamma1 and Gamma2
    ff    = 2.d0 / pi * fact * (gam1 * dcosh(e) + gam2 * fact2 * d)
    e     = dexp(e)
    p     = e  / (gampl * pi)
    q     = 1.d0 / (e * pi * gammi)
    pimu2 = 0.5d0 * pimu
    if (abs(pimu2) < eps) then
      fact3 = 1.d0
    else
      fact3 = dsin(pimu2) / pimu2
    endif ! end if statement
    r     = pi * pimu2 * fact3 * fact3
    c     = 1.d0
    d     = -x2 * x2
    summ  = ff + r * q
    summ1 = p
    do i = 1, maxit
      ff    = (i * ff + p + q) / (i * i - xmu2)
      c     = c * d / i
      p     = p / (i - xmu)
      q     = q / (i + xmu)
      del   = c * (ff + r * q)
      summ  = summ + del
      del1  = c * p - i * del
      summ1 = summ1 + del1
      if (abs(del) < (1.d0 + abs(summ)) * eps) then
        goto 2
      endif ! end if statement
    enddo
    print *, "bessY series failed to converge"
2   continue
    bessYmu  = -summ
    bessY1   = -summ1 * xi2
    bessYmup = xmu * xi * bessYmu - bessY1
    bessJmu  = w / (bessYmup - f * bessYmu)
  else ! evaluate the second continued fraction by Lentz's method
    a    = 0.25d0 - xmu2
    p    = -0.5d0 * xi
    q    = 1.d0
    br   = 2.d0 * x
    bi   = 2.d0
    fact = a * xi / (p * p + q * q)
    cr   = br + q * fact
    ci   = bi + p * fact
    den  = br * br + bi * bi
    dr   = br / den
    di   = -bi / den
    dlr  = cr * dr - ci * di
    dli  = cr * di + ci * dr
    temp = p * dlr - q * dli
    q    = p * dli+q * dlr
    p    = temp
    do i = 2, maxit
      a  = a + 2 * (i - 1)
      bi = bi + 2.d0
      dr = a * dr + br
      di = a * di + bi
      if ((abs(dr) + abs(di)) < floatingpointmin) then
        dr = floatingpointmin
      endif ! end if statement
      fact = a / (cr * cr +  ci * ci)
      cr   = br + cr * fact
      ci   = bi - ci * fact
      if ((abs(cr) + abs(ci)) < floatingpointmin) then
        cr = floatingpointmin
      endif ! end if statement
      den  = dr * dr + di * di
      dr   = dr / den
      di   = -di / den
      dlr  = cr * dr - ci * di
      dli  = cr * di + ci * dr
      temp = p * dlr - q * dli
      q    = p * dli + q * dlr
      p    = temp
      if ((abs(dlr - 1.d0) + abs(dli)) < eps) then
        goto 3
      endif ! end if statement
    enddo
    print *, "the second continued fraction failed in dbessJY"
3   continue
    gam    = (p - f) / q
    bessJmu  = dsqrt(w / ((p - f) * gam + q))
    bessJmu  = sign(bessJmu, bessJl)
    bessYmu  = bessJmu * gam
    bessYmup = bessYmu * (p + q / gam)
    bessY1   = xmu * xi * bessYmu - bessYmup
  endif ! end if statement
  fact     = bessJmu / bessJl
  bessJ    = negate * bessJl1 * fact
    ! scale original J_order and J_order′
  bessJder = bessJder1 * fact

  ! Here is the upward recurrence of Y_order
  do i = 1, nl
    bessYtemp = (xmu + i) * xi2 * bessY1 - bessYmu
    bessYmu   = bessY1
    bessY1    = bessYtemp
  enddo
  bessY    = negate * bessYmu
  bessYder = order * xi * bessYmu - bessY1

  return
end subroutine dbessJY

! This subroutine evaluates Gamma1 and Gamma2 by a Chebyshev
! expansion for |x| <= 1/2. Also returns 1 / Gamma(1 + x) and
! 1 / Gamma(1 − x). If converting to real, set nuse1 = 7 and
! nuse2 = 8.  This subroutine calls the chebev function, which is
! contained in this file.
! Subroutine Inputs/Output:
! x     - (input, real, [ ]) the point at which the gamma functions
!         are evaluated
! gam1  - (output, real, [ ]) the function value of Gamma1
! gam2  - (output, real, [ ]) the function value of Gamma2
! gampl - (output, real, [ ]) the function value of 1 / Gamma(1 + x)
! gammi - (output, real, [ ]) the function value of 1 / Gamma(1 − x)
subroutine beschbgamma(x, gam1, gam2, gampl, gammi)

  implicit none

  integer, parameter :: nuse1 = 7, nuse2 = 8
  real(8) :: gam1, gam2, gammi, gampl, x
  real(8) :: xx, c1(7), c2(8)
  real(8), external :: chebev
  save c1, c2
  data c1 /-1.142022680371168d0, 6.5165112670737d-3, &
    3.087090173086d-4, -3.4706269649d-6, 6.9437664d-9, &
    3.67795d-11, -1.356d-13/
  data c2 /1.843740587300905d0, -7.68528408447867d-2, &
    1.2719271366546d-3, -4.9717367042d-6, -3.31261198d-8, &
    2.423096d-10, -1.702d-13, -1.49d-15/

  xx    = 8.d0 * x * x - 1.d0
    ! x is multiplied by 2 in order to change its range to
    ! -1 to 1, and then another transformation is applied in order
    ! to evaluate using the even Chebyshev series; since the
    ! function is even it would be wasteful to call the chebev
    ! function with all the odd coefficients being zero; an
    ! approximation of an even function on [-1,1] will only
    ! involve even Chebyshev polynomials.  Once x is in the range
    ! -1 to 1, then the transformation being applied is
    ! T_2n(x) = T_n(2 * x**2 - 1)
  gam1  = chebev(-1.d0, 1.d0, c1, nuse1, xx)
  gam2  = chebev(-1.d0, 1.d0, c2, nuse2, xx)
  gampl = gam2 - x * gam1
  gammi = gam2 + x * gam1

  return
end subroutine beschbgamma

! This subroutine calls upon the subroutine beschbgamma defined in
! this same file.  This subroutine returns the modified Bessel
! functions bessI = I_order and bessK = K_order, and their
! respective derivatives bessIder = I_order′ and
! bessKder = K_order′, for positive x and for an order greater or
! equal to 0.  The relative accuracy is within one or two
! significant digits of eps. The parameter "floatingpointmin" is
! a number close to the machine’s smallest floating point number.
! All internal arithmetic is accomplished with real numbers. In
! order to convert the entire routine to double precision, use
! the real(8) declaration and decrease eps to 10**−16.  Also
! convert the subroutine beschbgamma.
! Subroutine Inputs/Output:
! x        - (input, real, [ ]) the given point at which the Bessel
!         functions will be calculated
! order    - (input, real, [ ]) the order of the Bessel functions
! bessI    - (output, real, [ ]) the value of the Bessel I function
! bessK    - (output, real, [ ]) the value of the Bessel K function
! bessIder - (output, real, [ ]) the value of the derivative of the
!         Bessel I function
! bessKder - (output, real, [ ]) the value of the derivative of the
!         Bessel K function
subroutine dbessIK(x, order, bessI, bessK, bessIder, bessKder)

  implicit none

  integer, parameter :: maxit = 10000
  real(8) :: bessI, bessIder, bessK, bessKder, x, order
  real(8), parameter :: eps = 1.d-16, floatingpointmin = 1.e-30, &
    pi = 3.141592653589793d0, xmin = 2.d0
  integer :: i, l, nl
  real(8) :: a, a1, b, c, d, del, del1, delh, dels, e, f, &
    fact, fact2, ff, gam1, gam2, gammi, gampl, h, p, pimu, &
    q, q1, q2, qnew, bessIl, bessIl1, bessImu, bessIder1, &
    bessIderl, bessItemp, bessK1, bessKmu, bessKmup, &
    bessKtemp, s, summ, summ1, x2, xi, xi2, xmu, xmu2

  if ((x <= 0.d0) .or. (order < 0.d0)) then
    print *, "Bad arguments were inputed into dbessIK."
    print *, "x = ", x
    print *, "order = ", order
  endif ! end if statement

  ! Here "nl" is the number of downward recurrences of the I’s and
  ! upward recurrences of K’s. "xmu" lies between −1/2 and 1/2:
  nl   = int(order + 0.5d0)
  xmu  = order - nl
  xmu2 = xmu * xmu
  xi   = 1.d0 / x
  xi2  = 2.d0 * xi
  h    = order * xi

  ! Here the first continued fraction is found by the modified
  ! Lentz's method:
  if (h < floatingpointmin) then
    h = floatingpointmin
  endif ! end if statement
  b = xi2 * order
  d = 0.d0
  c = h
  do i = 1, maxit
    b   = b + xi2
    d   = 1.d0 / (b + d)
      ! denominators cannot be zero here, so there is no need for
      ! special precautions
    c   = b + 1.d0 / c
    del = c * d
    h   = del * h
    if (abs(del - 1.d0) < eps) then
      goto 10
    endif ! end if statement
  enddo
  print *, "x is too large in dbessIK; try an asymptotic expansion"
  print *, "x = ", x
10  continue

  ! Here I_order and I_order' are initialized for downward
  ! recurrence:
  bessIl    = floatingpointmin
  bessIderl = h * bessIl
  bessIl1   = bessIl
  bessIder1 = bessIderl
    ! store values for future rescaling
  fact      = order * xi
  do l = nl, 1, -1
    bessItemp = fact * bessIl + bessIderl
    fact      = fact - xi
    bessIderl = fact * bessItemp + bessIl
    bessIl    = bessItemp
  enddo

  ! The subroutine already has calculated unnormalized I_order and
  ! I_order':
  f = bessIderl / bessIl
  if (x <= xmin) then ! use series
    x2   = 0.5d0 * x
    pimu = pi * xmu
    if (abs(pimu) < eps) then
      fact = 1.d0
    else
      fact = pimu / dsin(pimu)
    endif ! end if statement
    d = -dlog(x2)
    e = xmu * d
    if (abs(e) < eps) then
      fact2 = 1.d0
    else
      fact2 = dsinh(e) / e
    endif ! end if statement
    call beschbgamma(xmu, gam1, gam2, gampl, gammi)
      ! Chebyshev evaluation of Gamma1 and Gamma2
    ff    = fact * (gam1 * dcosh(e) + gam2 * fact2 * d)
    summ  = ff
    e     = dexp(e)
    p     = 0.5d0 * e / gampl
    q     = 0.5d0 / (e * gammi)
    c     = 1.d0
    d     = x2 * x2
    summ1 = p
    do i = 1, maxit
      ff    = (i * ff + p + q) / (i * i - xmu2)
      c     = c * d / i
      p     = p / (i - xmu)
      q     = q / (i + xmu)
      del   = c * ff
      summ  = summ + del
      del1  = c * (p - i * ff)
      summ1 = summ1 + del1
      if (abs(del) < (abs(summ) * eps)) then
        goto 12
      endif ! end if statement
    enddo
    print *, "bessK series failed to converge"
12    continue
    bessKmu = summ
    bessK1  = summ1 * xi2
  else ! evaluate the second continued fraction by Steed's algorithm;
     ! note that there can be no zero denominators
    b    = 2.d0 * (1.d0 + x)
    d    = 1.d0 / b
    delh = d
    h    = delh
    q1   = 0.d0       ! initializing for recurrences
    q2   = 1.d0
    a1   = 0.25d0 - xmu2
    c    = a1
    q    = c
    a    = -a1
    s    = 1.d0 + q * delh
    do i = 2, maxit
      a    = a - 2 * (i - 1)
      c    = -a * c / i
      qnew = (q1 - b * q2) / a
      q1   = q2
      q2   = qnew
      q    = q + c * qnew
      b    = b + 2.d0
      d    = 1.d0 / (b + a * d)
      delh = (b * d - 1.d0) * delh
      h    = h + delh
      dels = q * delh
      s    = s + dels
      ! Here because the second continued fraction converges
      ! faster, the convergence of the sum is the only test
      ! needed:
      if (abs(dels / s) < eps) then
        goto 13
      endif ! end if statement
    enddo
    print *, "second continued fraction failed to converge in dbessIK"
13    continue
    h       = a1 * h
    bessKmu = dsqrt(pi / (2.d0 * x)) * dexp(-x) / s
      ! the exp(-x) factor has been omitted in order to rescale
      ! the modified Bessel functions by exp(x) for x >= xmin
    bessK1  = bessKmu * (xmu + x + 0.5d0 - h) * xi
  endif ! end if statement
  bessKmup = xmu * xi * bessKmu - bessK1
  bessImu  = xi  / (f * bessKmu - bessKmup)
    ! I_order is determined from the Wronskian
  bessI    = (bessImu * bessIl1) / bessIl
    ! rescaling of I_order and I_order'
  bessIder = (bessImu * bessIder1) / bessIl

  ! Here are the upward recurrences of K_order
  do i = 1, nl
    bessKtemp = (xmu + i) * xi2 * bessK1 + bessKmu
    bessKmu   = bessK1
    bessK1    = bessKtemp
  enddo
  bessK    = bessKmu
  bessKder = order * xi * bessKmu - bessK1

  return
end subroutine dbessIK

! This function evaluates the Chebyshev polynomial
! sum_(k = 1)^(m) c_k * T_(k - 1)(y) - c_1 / 2 at the point
! y = [x - (b + a) / 2] / [(b - a) / 2], where c(m) is an array of
! the Chebyshev coefficients.
! Subroutine Inputs/Output:
! a - (input, real, [ ]) lower endpoint of the interval in which
!     the Chebyshev polynomial is evaluated
! b - (input, real, [ ]) upper endpoint of the interval in which the
!     Chebyshev polynomial is evaluated
! c - (input, real vector (m), [ ]) array of the Chebyshev
!       coefficients
! m - (input, integer, [ ]) number of Chebyshev coefficients to be
!     calculated
! x - (input, real, [ ]) the point that decides where the Chebyshev
!     polynomial is evaluated
function chebev(a, b, c, m, x)

  integer, intent(in) :: m
  real(8), intent(in) :: a, b, x, c(m)
  real(8) :: chebev
  real(8) :: d, dd, sv, y, y2
  integer :: j

  if (((x - a) * (x - b)) > 0.d0) then
    print *, "x must be greater than a and b"
    print *, "x = ", x
    print *, "a = ", a, " b = ", b
  endif ! end if statement
  d  = 0.d0
  dd = 0.d0
  y  = (2.d0 * x - a - b) / (b - a)
    ! change of variable
  y2 = 2.d0 * y

  ! Here is Clenshaw's recurrence algorithm:
  do j = m, 2, -1
    sv = d
    d  = y2 * d - dd + c(j)
    dd = sv
  enddo
  chebev = y * d - dd + 0.5d0 * c(1)

  return
end function chebev

subroutine d2bessJY(x, order, d2bessJ, d2bessY)
  implicit none
  real*8, intent(in) :: x, order
  real*8, intent(out) :: d2bessJ, d2bessY
  real *8  :: bessJ_1, bessY_1, dbessJ_1, dbessY_1
  real *8  :: bessJ_2, bessY_2, dbessJ_2, dbessY_2
  if(order.lt.1.d0) then
    write(*,*) 'error from d2bessJY: order must be >=1 '
    stop
  endif
  call dbessJY(x, order-1.d0, bessJ_1, bessY_1, dbessJ_1, dbessY_1)
  call dbessJY(x, order+1.d0, bessJ_2, bessY_2, dbessJ_2, dbessY_2)
  d2bessJ = 0.5d0*(dbessJ_1-dbessJ_2)
  d2bessY = 0.5d0*(dbessY_1-dbessY_2)

end subroutine d2bessJY


subroutine d2bessIK(x, order, d2bessI, d2bessK)
  implicit none
  real*8, intent(in) :: x, order
  real*8, intent(out) :: d2bessI, d2bessK
  real *8  :: bessI_1, bessK_1, dbessI_1, dbessK_1
  real *8  :: bessI_2, bessK_2, dbessI_2, dbessK_2
  if(order.lt.1.d0) then
    write(*,*) 'error from d2bessIK: order must be >=1 '
    stop
  endif
  call dbessIK(x, order-1.d0, bessI_1, bessK_1, dbessI_1, dbessK_1)
  call dbessIK(x, order+1.d0, bessI_2, bessK_2, dbessI_2, dbessK_2)
  d2bessI =  0.5d0*(dbessI_1+dbessI_2)
  d2bessK = -0.5d0*(dbessK_1+dbessK_2)
  !write(*,*)'from d2bessIK: dbessI_1, dbessI_2 = ',dbessI_1, dbessI_2
end subroutine d2bessIK
