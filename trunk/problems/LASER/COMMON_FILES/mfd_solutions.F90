!--------------------------------------------------------------------------------
!> Purpose : exact (manufactured) solution
!
!> @param[in]  Xp  - a point in physical space
!> @param[in]  Fld - 0: pump field, 1: signal field
!> @param[out] E   - value of the solution (one electric field component)
!> @param[out] dE  - corresponding first derivatives
!> @param[out] d2E - corresponding second derivatives
!--------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine mfd_solutions(Xp,Fld, E,dE,d2E)
!
   use data_structure3D
   use control, only : GEOM_TOL, NEXACT
   use commonParam
   use laserParam
!
   implicit none
!--------------------------------------------------------------------------------
   real(8), dimension(3),  intent(in)  :: Xp
   integer,                intent(in)  :: Fld
   VTYPE,                  intent(out) :: E   ! solution
   VTYPE,  dimension(3),   intent(out) :: dE  ! 1st derivative
   VTYPE,  dimension(3,3), intent(out) :: d2E ! 2nd derivative
!--------------------------------------------------------------------------------
   integer :: modified
!
   real(8) :: x1,x2,x3
   real(8) :: nn        ! for EXPONENTIAL solution
   real(8) :: cn,dn     ! for singular solution
   real(8) :: np_x,np_y,np_z,r0,k0,w0,phase,ampl,om
   VTYPE   :: f_x,f_y,f_z,df_x,df_y,df_z,ddf_x,ddf_y,ddf_z
!..for TE modes
   real(8) :: gammaTE10
   real(8) :: gammaTE20
   real(8) :: r
!..for LP Modes
   real(8) :: gamm, beta, k, ca, cb, cc, r_x, r_y
   real(8) :: BESSEL_dJ1, BESSEL_K0, BESSEL_dK0, BESSEL_K1, BESSEL_dK1
!..for bessel function modes
   real(8) :: zeta,rho,xi,theta,rcore,rcladding,alpha,alpha_scale
   real(8) :: order, bessJ, bessY, dbessJ, dbessY,bessJc, bessYc, dbessJc, dbessYc
   real(8) :: bessI, bessK, dbessI, dbessK, bessIc, bessKc, dbessIc, dbessKc
   real(8) :: d2bessJ,d2bessY,d2bessI,d2bessK
   real(8) :: hess1(2,2),hess2(2,2)
   real(8) :: angular,angular_x,angular_y,angular_xx,angular_xy,angular_yy
   real(8) :: Jm,Jm_x,Jm_y,Jm_xx,Jm_xy,Jm_yy
   real(8) :: Km,Km_x,Km_y,Km_xx,Km_xy,Km_yy
!
   VTYPE :: c2z,uz,uz_x,uz_y,uz_z,uz_xx,uz_xy,uz_xz,uz_yy,uz_yx,uz_yz
   VTYPE :: uz_zz,uz_zy,uz_zx
   VTYPE :: pz,pz_x,pz_y,pz_z,pz_xx,pz_xy,pz_xz,pz_yy,pz_yx,pz_yz
   VTYPE :: pz_zz,pz_zy,pz_zx, zbeta,zdbeta,zd2beta
!
   VTYPE :: E01,E11,E21,E02
   VTYPE :: dE01(3),dE11(3),dE21(3),dE02(3)
!
!--------------------------------------------------------------------------------
!
!..initialize variables
   E = ZERO; dE = ZERO; d2E = ZERO;
   f_x = ZERO; f_y = ZERO; f_z = ZERO
   df_x = ZERO; df_y = ZERO; df_z = ZERO
   ddf_x = ZERO; ddf_y = ZERO; ddf_z = ZERO
!
!--------------------------------------------------------------------------------
!      D E C L A R E    S O L U T I O N    V A R I A B L E S                    |
!--------------------------------------------------------------------------------
!
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
   if (Fld.ne.0 .and. Fld.ne.1) then
      write(*,*) 'mfd_solutions: Fld_flag = ',Fld,'. stop.'
      return
   endif
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
!--------------- 2nd prob -------------------------------------------------------
!..a smooth solution
   elseif (ISOL .eq. 2) then
!
      f_x = sin(OMEGA*x1)
      f_y = sin(OMEGA*x2)
      f_z = sin(OMEGA*x3)
!
!  ...1st order derivatives
      df_x = OMEGA*cos(OMEGA*x1)
      df_y = OMEGA*cos(OMEGA*x2)
      df_z = OMEGA*cos(OMEGA*x3)
!
!  ...2nd order derivatives
      ddf_x = -OMEGA**2*f_x
      ddf_y = -OMEGA**2*f_y
      ddf_z = -OMEGA**2*f_z
!
!--------------- 3rd prob -------------------------------------------------------
!..a plane wave
   elseif (ISOL .eq. 3) then
!
      !f_x = ZI + cos(x1)
      f_x = 1.d0
      !f_y = ZI + sin(x2)
      f_y = 1.d0
      f_z = 2.d0*exp(-ZI*OMEGA*x3)
!
!  ...1st order derivatives
      !df_x = -sin(x1)
      df_x = 0.d0
      !df_y = cos(x2)
      df_y = 0.d0
      df_z = (-ZI*OMEGA)*f_z
!
!  ...2nd order derivatives
      !ddf_x = -cos(x1)
      ddf_x = 0.d0
      !ddf_y = -sin(x2)
      ddf_y = 0.d0
      ddf_z = -OMEGA**2*f_z
!
!--------------- 4th prob -------------------------------------------------------
!..an exponential
   elseif (ISOL .eq. 4) then
!
      nn = -OMEGA
!  ...value
      f_x=exp(-Xp(1)**2/2)
      f_y=exp(-Xp(2)**2/2)
      f_z=1.d0 !exp(nn*Xp(3))
!
!  ...1st order derivatives
      df_x=-Xp(1)*f_x
      df_y=-Xp(2)*f_y
      df_z=0.d0!nn*f_z
!
!  ...2nd order derivatives
      ddf_x=-f_x-(Xp(1)*df_x)
      ddf_y=-f_y-(Xp(2)*df_y)
      ddf_z=0.d0 !nn*df_z
!
!--------------- 5th prob -------------------------------------------------------
!..fundamental TE10 mode for rectangular waveguide (signal)
   elseif (ISOL .eq. 5) then
!
!     ! shifted source (s.t. fields are non-zero at z=0)
      x1 = Xp(1); x2 = Xp(2); x3 = Xp(3) - .5d0;
!
!        nn = 1: TE10 mode only
!        nn = 2: TE20 mode only
!        nn = 3: TE10 + TE20 mode
!        note OMEGA must be set correctly in set_env routine
         nn = 1
!
      if (nn .eq. 1 .or. nn .eq. 3) then
!     ...fundamental mode TE10
         gammaTE10 = sqrt(1.d0-(PI**2)/(OMEGA**2))
!
         f_x=-ZI*(OMEGA/PI)*sin(PI*x1)
         f_y=1.d0
         f_z=exp(-ZI*OMEGA*x3*gammaTE10)
!
!     ...1st order derivatives
         df_x=-ZI*OMEGA*cos(PI*x1)
         df_y=0.d0
         df_z=-(ZI*OMEGA*gammaTE10)*f_z
!
!     ...2nd order derivatives
         ddf_x=-PI**2*f_x
         ddf_y=0.d0
         ddf_z=-(ZI*OMEGA*gammaTE10)*df_z
!
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
      endif
      if (nn .eq. 2 .or. nn .eq. 3) then
!
!     ...next higher mode (TE20, same cutoff as TE01)
         gammaTE20 = sqrt(1.d0-((2.d0*PI)**2)/(OMEGA**2))
!
         f_x=-ZI*(OMEGA/(2.d0*PI))*sin(2.d0*PI*x1)
         f_y=1.d0
         f_z=exp(-ZI*OMEGA*x3*gammaTE20)
!
!     ...1st order derivatives
         df_x=-ZI*OMEGA*cos(2.d0*PI*x1)
         df_y=0.d0
         df_z=-(ZI*OMEGA*gammaTE20)*f_z
!
!     ...2nd order derivatives
         ddf_x=-(2.d0*PI)**2*f_x
         ddf_y=0.d0
         ddf_z=-(ZI*OMEGA*gammaTE20)*df_z
!
         E=E+f_x*f_y*f_z
!     ...1st order derivatives
         dE(1)=dE(1)+  df_x *   f_y *   f_z
         dE(2)=dE(2)+   f_x *  df_y *   f_z
         dE(3)=dE(3)+   f_x *   f_y *  df_z
!
!     ...2nd order derivatives
         d2E(1,1) =d2E(1,1)+ ddf_x *   f_y *   f_z
         d2E(1,2) =d2E(1,2)+  df_x *  df_y *   f_z
         d2E(1,3) =d2E(1,3)+  df_x *   f_y *  df_z
         d2E(2,1) = d2E(1,2)
         d2E(2,2) =d2E(2,2)+   f_x * ddf_y *   f_z
         d2E(2,3) =d2E(2,3)+   f_x *  df_y *  df_z
         d2E(3,1) = d2E(1,3)
         d2E(3,2) = d2E(2,3)
         d2E(3,3) =d2E(3,3)+   f_x *   f_y * ddf_z
      endif
!
!--------------- 9th prob -------------------------------------------------------
!..fundamental TE10 mode for rectangular waveguide (signal) with exponential growth
   elseif (ISOL .eq. 9) then
!  ...shifted source (s.t. fields are non-zero at z=0)
      x1 = Xp(1); x2 = Xp(2); x3 = Xp(3);! - .5d0;
!
      f_x=-ZI*(OMEGA/PI)*sin(PI*x1)
      f_y=1.d0
      f_z=exp((EXP_COEFF-ZI*OMEGA*GAMMA)*x3)
!
!  ...1st order derivatives
      df_x=-ZI*OMEGA*cos(PI*x1)
      df_y=0.d0
      df_z=(EXP_COEFF-ZI*OMEGA*GAMMA)*f_z
!
!  ...2nd order derivatives
      ddf_x=-PI**2*f_x
      ddf_y=0.d0
      ddf_z=(EXP_COEFF-ZI*OMEGA*GAMMA)*df_z
!
!--------------- 10th prob -------------------------------------------------------
!..fundamental TE10 mode for rectangular waveguide with exponential growth
!..with complex stretching
   elseif (ISOL .eq. 10) then
!
      x1 = Xp(1); x2 = Xp(2); x3 = Xp(3);
!
      call get_Beta(Xp,Fld, zbeta,zdbeta,zd2beta)
!
      f_x=-ZI*(OMEGA/PI)*sin(PI*x1)
      f_y= 1.d0
      f_z=exp((EXP_COEFF-ZI*OMEGA*GAMMA)*zbeta)
!
!  ...1st order derivatives
      df_x=-ZI*OMEGA*cos(PI*x1)
      df_y=0.d0
      df_z=(EXP_COEFF-ZI*OMEGA*GAMMA)*zdbeta*f_z
!
!  ...2nd order derivatives
      ddf_x=-PI**2*f_x
      ddf_y=0.d0
      ddf_z=(EXP_COEFF-ZI*OMEGA*GAMMA)*zd2beta*f_z + &
            (EXP_COEFF-ZI*OMEGA*GAMMA)*zdbeta*df_z
!
!--------------- 11th prob -------------------------------------------------------
!..f(x,y,z) = sin(OMEGA*x)sin(OMEGA*y)
   elseif (ISOL .eq. 11) then

      f_x= sin(OMEGA*x1)
      f_y= sin(OMEGA*x2)
      f_z= 1.d0
!
!     1st order derivatives
      df_x=OMEGA*cos(x1)
      df_y=OMEGA*cos(x2)
      df_z=0.d0
!
!     2nd order derivatives
      ddf_x=-OMEGA**2*f_x
      ddf_y=-OMEGA**2*f_y
      ddf_z=0.d0
!
!--------------- 12th prob -------------------------------------------------------
!..f(x,y,z) = cos(x)cos(y)
   elseif (ISOL .eq. 12) then
!
!  ...TM02 mode in circular waveguide, for a = sqrt(2), omega=4.5
      r = sqrt(x1*x1+x2*x2)
      if (r .lt. 1.0d-13) then
         E = 0.d0
      else
         if (ICOMP_TS .eq. 1) then
            E = ZI*BESSEL_J1(3.9d0*r)*x1/r
         elseif (ICOMP_TS .eq. 2) then
            E = ZI*BESSEL_J1(3.9d0*r)*x2/r
         endif
      endif
!     ...TODO DERIVATIVES dE(1:3)
!
!--------------- 13th prob -------------------------------------------------------
!..Fundamental mode LP01 in dielectric waveguide
   elseif (ISOL .eq. 13) then
!
!  ...LP01 in dielectric waveguide, a = sqrt(2), omega=25.7
!      k    = 37.2854d0
!      gamm =  1.16301d0
!      beta =  1.23370d0
!      ampl =  1.0d0
!
!  ...LMA fiber
!  ...LP01 (signal) in dielectric waveguide, a = 0.9*sqrt(2), omega=2*pi/0.1064=59.0525
      if (Fld .eq. 1) then
         ampl =  1.0d0
         ! Signal laser frequency (active gain)
         if (LAMBDA_SIGNAL .eq. 1064.0d-9/L_0) then
            if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.4512d0 .and. CLAD_NX.eq.1.4500d0) .or.   &
                (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.4512d0 .and. CLAD_NY.eq.1.4500d0)) then
               k    = 85.6833d0
               gamm =  1.53131d0
               beta =  3.12978d0
            else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.1520d0 .and. CLAD_NX.eq.1.1500d0) .or.  &
                     (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.1520d0 .and. CLAD_NY.eq.1.1500d0)) then
               k    = 68.0103d0
               gamm =  1.57221d0
               beta =  3.68554d0
            else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.6510d0 .and. CLAD_NX.eq.1.6500d0) .or.  &
                     (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.6510d0 .and. CLAD_NY.eq.1.6500d0)) then
               k    = 97.4838d0
               gamm =  1.52302d0
               beta =  3.03177d0
            else
               write(*,*) 'mfd_solutions: ISOL 13, unexpected case 1. stop.'
               stop
            endif
         ! Signal Raman frequency
         else if (LAMBDA_SIGNAL .eq. 1116.0d-9/L_0) then
            if (ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.4512d0 .and. CLAD_NX.eq.1.4500d0) then
               k    = 81.6899d0
               gamm =  1.51632d0
               beta =  2.95571d0
            else
               write(*,*) 'mfd_solutions: ISOL 13, unexpected case 2. stop.'
               stop
            endif
         else
            write(*,*) 'mfd_solutions: ISOL 13, unexpected case 3. stop.'
            stop
         endif
      endif
!
!  ...LP01 (pump) in dielectric waveguide, a = 0.9*sqrt(2), omega=2*pi/0.0976
      if (Fld .eq. 0) then
         ampl =  2.0d0
         ! Pump frequency (active gain)
         if (LAMBDA_PUMP .eq. 976.0d-9/L_0) then
            if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.4512d0 .and. CLAD_NX.eq.1.4500d0) .or.   &
                (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.4512d0 .and. CLAD_NY.eq.1.4500d0)) then
               k    = 93.4108d0
               gamm =  1.55709d0
               beta =  3.46466d0
            else
               write(*,*) 'mfd_solutions: ISOL 13, unexpected case 4. stop.'
               stop
            endif
         ! Pump frequency (Raman gain)
         else if (LAMBDA_PUMP .eq. 1064.0d-9/L_0) then
            if (ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.4512d0 .and. CLAD_NX.eq.1.4500d0) then
               k    = 85.6833d0
               gamm =  1.53131d0
               beta =  3.12978d0
            else
               write(*,*) 'mfd_solutions: ISOL 13, unexpected case 5. stop.'
               stop
            endif
         else
            write(*,*) 'mfd_solutions: ISOL 13, unexpected case 6. stop.'
            stop
         endif
      endif
!
      call get_LP01(Xp,ampl,k,gamm,beta, E,dE)
!
!--------------- 14th prob -------------------------------------------------------
!..LP11 mode in dielectric waveguide
   elseif (ISOL .eq. 14 .or. ISOL .eq. 140) then
!
      if (Fld .ne. 1) then
         write(*,*) 'mfd_solutions: unexpected Fld_flag. stop.'
         stop
      endif
!
!  ...LMA fiber
!  ...LP11 (signal) in dielectric waveguide, a = 0.9*sqrt(2), omega=2*pi/0.1064=59.0525
      ampl =  1.0d0
      if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.4512d0 .and. CLAD_NX.eq.1.4500d0) .or.   &
          (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.4512d0 .and. CLAD_NY.eq.1.4500d0)) then
         k    = 85.6630d0
         gamm =  2.41319d0
         beta =  2.51336d0
      else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.1520d0 .and. CLAD_NX.eq.1.1500d0) .or.  &
               (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.1520d0 .and. CLAD_NY.eq.1.1500d0)) then
         k    = 67.9830d0
         gamm = 2.48683d0
         beta = 3.14177d0
      else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.6510d0 .and. CLAD_NX.eq.1.6500d0) .or.  &
               (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.6510d0 .and. CLAD_NY.eq.1.6500d0)) then
         k    = 97.4662d0
         gamm =  2.39795d0
         beta =  2.40022d0
      else
         write(*,*) 'mfd_solutions: ISOL 14, unexpected case. stop.'
         stop
      endif
      select case(ISOL)
         case(14) ; call get_LP11a(Xp,ampl,k,gamm,beta, E,dE)
         case(140); call get_LP11b(Xp,ampl,k,gamm,beta, E,dE) ! rotated by 90 degrees
      end select
!
!--------------- 15th prob -------------------------------------------------------
!..LP21 mode in dielectric waveguide
   elseif (ISOL .eq. 15 .or. ISOL .eq. 150) then
!
      if (Fld .ne. 1) then
         write(*,*) 'mfd_solutions: unexpected Fld_flag. stop.'
         stop
      endif
!
!  ...LMA fiber
!  ...LP21 (signal) in dielectric waveguide, a = 0.9*sqrt(2), omega=2*pi/0.1064=59.0525
      ampl =  1.0d0
      if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.4512d0 .and. CLAD_NX.eq.1.4500d0) .or.   &
          (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.4512d0 .and. CLAD_NY.eq.1.4500d0)) then
         k    = 85.6380d0
         gamm =  3.17859d0
         beta =  1.42726d0
      else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.1520d0 .and. CLAD_NX.eq.1.1500d0) .or.  &
               (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.1520d0 .and. CLAD_NY.eq.1.1500d0)) then
         k    = 67.9484d0
         gamm = 3.29870d0
         beta = 2.27456d0
      else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.6510d0 .and. CLAD_NX.eq.1.6500d0) .or.  &
               (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.6510d0 .and. CLAD_NY.eq.1.6500d0)) then
         k    = 97.4447d0
         gamm =  3.15230d0
         beta =  1.25468d0
      else
         write(*,*) 'mfd_solutions: ISOL 15, unexpected case. stop.'
         stop
      endif
      select case(ISOL)
         case(15) ; call get_LP21a(Xp,ampl,k,gamm,beta, E,dE)
         case(150); call get_LP21b(Xp,ampl,k,gamm,beta, E,dE) ! rotated by 45 degrees
      end select
!
!--------------- 16th prob -------------------------------------------------------
!..LP02 mode in dielectric waveguide
   elseif (ISOL .eq. 16) then
!
      if (Fld .ne. 1) then
         write(*,*) 'mfd_solutions: unexpected Fld_flag. stop.'
         stop
      endif
!
!  ...LMA fiber
!  ...LP02 (signal) in dielectric waveguide, a = 0.9*sqrt(2), omega=2*pi/0.1064=59.0525
      ampl =  1.0d0
      if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.4512d0 .and. CLAD_NX.eq.1.4500d0) .or.   &
          (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.4512d0 .and. CLAD_NY.eq.1.4500d0)) then
         k    = 85.6322d0
         gamm =  3.33123d0
         beta =  1.02145d0
      else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.1520d0 .and. CLAD_NX.eq.1.1500d0) .or.  &
               (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.1520d0 .and. CLAD_NY.eq.1.1500d0)) then
         k    = 67.9384d0
         gamm = 3.50008d0
         beta = 1.95051d0
      else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.6510d0 .and. CLAD_NX.eq.1.6500d0) .or.  &
               (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.6510d0 .and. CLAD_NY.eq.1.6500d0)) then
         k    = 97.4401d0
         gamm =  3.28999d0
         beta =  0.82896d0
      else
         write(*,*) 'mfd_solutions: ISOL 16, unexpected case. stop.'
         stop
      endif
      call get_LP02(Xp,ampl,k,gamm,beta, E,dE)
!
!--------------- 17th prob -------------------------------------------------------
!..Mixed LP01/LP11/LP21/LP02 mode in dielectric waveguide
   elseif (ISOL .eq. 17) then
!
      if (Fld .ne. 1) then
         write(*,*) 'mfd_solutions: unexpected Fld_flag. stop.'
         stop
      endif
!
      E01 = 0.d0; dE01 = 0.d0
      E11 = 0.d0; dE11 = 0.d0
      E21 = 0.d0; dE21 = 0.d0
      E02 = 0.d0; dE02 = 0.d0
!
      if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.4512d0 .and. CLAD_NX.eq.1.4500d0) .or.   &
          (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.4512d0 .and. CLAD_NY.eq.1.4500d0)) then
!
!     ...LMA fiber, a = 0.9*sqrt(2), omega=2*pi/0.1064=59.0525
!     ...LP01 (signal)
         k    = 85.6833d0
         gamm =  1.53131d0
         beta =  3.12978d0
         ampl =  1.0d0 * 0.725d0 ! adjust so that the input is ca. 25 W
!     ...power input
         !ampl = ampl * sqrt(2.0d0) !  50 W
         !ampl = ampl * sqrt(4.0d0) ! 100 W
         !ampl = ampl * sqrt(8.0d0) ! 200 W
!     ...seed ratio
         ampl = ampl * sqrt(0.9d0) ! 90% LP01
         !ampl = ampl * sqrt(0.8d0) ! 80% LP01
         !ampl = ampl * sqrt(0.7d0) ! 70% LP01
         !ampl = ampl * sqrt(0.6d0) ! 60% LP01
         !ampl = ampl * sqrt(0.5d0) ! 50% LP01
         call get_LP01(Xp,ampl,k,gamm,beta, E01,dE01)
!
!     ...LP11 (signal)
         k    = 85.6630d0
         gamm =  2.41319d0
         beta =  2.51336d0
         ampl =  1.0d0 * 1.577d0 ! adjust so that the input is ca. 25 W
!     ...power input
         !ampl = ampl * sqrt(2.0d0) !  50 W
         !ampl = ampl * sqrt(4.0d0) ! 100 W
         !ampl = ampl * sqrt(8.0d0) ! 200 W
!     ...seed ratio
         ampl = ampl * sqrt(0.1d0) ! 10% LP11
         !ampl = ampl * sqrt(0.2d0) ! 20% LP11
         !ampl = ampl * sqrt(0.3d0) ! 30% LP11
         !ampl = ampl * sqrt(0.4d0) ! 40% LP11
         !ampl = ampl * sqrt(0.5d0) ! 50% LP11
         call get_LP11a(Xp,ampl,k,gamm,beta, E11,dE11)
         !call get_LP11b(Xp,ampl,k,gamm,beta, E11,dE11)
!
!     ...LP21 (signal)
         k    = 85.6380d0
         gamm =  3.17859d0
         beta =  1.42726d0
         ampl =  1.0d0 * 1.982d0 ! adjust so that the input is ca. 25 W
!     ...power input
         !ampl = ampl * sqrt(2.0d0) !  50 W
         !ampl = ampl * sqrt(4.0d0) ! 100 W
         !ampl = ampl * sqrt(8.0d0) ! 200 W
!     ...seed ratio
         !ampl = ampl * sqrt(0.1d0) ! 10% LP21
         !ampl = ampl * sqrt(0.2d0) ! 20% LP21
         !ampl = ampl * sqrt(0.3d0) ! 30% LP21
         !ampl = ampl * sqrt(0.4d0) ! 40% LP21
         !ampl = ampl * sqrt(0.5d0) ! 50% LP21
         !call get_LP21a(Xp,ampl,k,gamm,beta, E21,dE21)
         !call get_LP21b(Xp,ampl,k,gamm,beta, E21,dE21)
!
!     ...LP02 (signal)
         k    = 85.6322d0
         gamm =  3.33123d0
         beta =  1.02145d0
         ampl =  1.0d0 * 1.315d0    ! adjust so that the input is ca. 25 W
!     ...power input
         !ampl = ampl * sqrt(2.0d0) !  50 W
         !ampl = ampl * sqrt(4.0d0) ! 100 W
         !ampl = ampl * sqrt(8.0d0) ! 200 W
!     ...seed ratio
         !ampl = ampl * sqrt(0.1d0) ! 10% LP02
         !ampl = ampl * sqrt(0.2d0) ! 20% LP02
         !ampl = ampl * sqrt(0.3d0) ! 30% LP02
         !ampl = ampl * sqrt(0.4d0) ! 40% LP02
         !ampl = ampl * sqrt(0.5d0) ! 50% LP02
         !call get_LP02(Xp,ampl,k,gamm,beta, E02,dE02)
!
      else if ((ICOMP_EXACT.eq.1 .and. CORE_NX.eq.1.1520d0 .and. CLAD_NX.eq.1.1500d0) .or.   &
               (ICOMP_EXACT.eq.2 .and. CORE_NY.eq.1.1520d0 .and. CLAD_NY.eq.1.1500d0)) then
!
!     ...LMA fiber, a = 0.9*sqrt(2), omega=2*pi/0.1064=59.0525
!     ...LP01 (signal)
         k    = 68.0103d0
         gamm = 1.57221d0
         beta = 3.68554d0
         ampl =  1.0d0
         call get_LP01(Xp,ampl,k,gamm,beta, E01,dE01)
!
!     ...LP11 (signal)
         k    = 67.9830d0
         gamm = 2.48683d0
         beta = 3.14177d0
         ampl =  1.0d0
!        call get_LP11a(Xp,ampl,k,gamm,beta, E11,dE11)
!        call get_LP11b(Xp,ampl,k,gamm,beta, E11,dE11)
!
!     ...LP21 (signal)
         k    = 67.9484d0
         gamm = 3.29870d0
         beta = 2.27456d0
         ampl =  1.0d0
!        call get_LP21a(Xp,ampl,k,gamm,beta, E21,dE21)
!        call get_LP21a(Xp,ampl,k,gamm,beta, E21,dE21)
!
!     ...LP02 (signal)
         k    = 67.9384d0
         gamm = 3.50008d0
         beta = 1.95051d0
         ampl =  1.0d0
         call get_LP02(Xp,ampl,k,gamm,beta, E02,dE02)
!
      endif
!
      E  =  E01+ E11+ E21+ E02
      dE = dE01+dE11+dE21+dE02
!
!--------------- 18th prob -------------------------------------------------------
!..LP12 mode in dielectric waveguide
   elseif (ISOL .eq. 18) then
!
      if (Fld .ne. 1) then
         write(*,*) 'mfd_solutions: unexpected Fld_flag. stop.'
         stop
      endif
!
!  ...LMA fiber
!  ...LP12 (signal) in dielectric waveguide, a = 0.9*sqrt(2), omega=2*pi/0.1064=59.0525
      k    = 85.5944d0
      gamm =  4.19098d0
      beta =  0.25d0
      ampl =  1.0d0
      call get_LP11a(Xp,ampl,k,gamm,beta, E,dE)
      !call get_LP11b(Xp,ampl,k,gamm,beta, E,dE)
!
!
!--------------- 19th prob -------------------------------------------------------
!.."Plane wave" in core of dielectric waveguide
   elseif (ISOL .eq. 19) then
!
      r = sqrt(x1*x1+x2*x2)
      if (r .eq. 0.d0) then
         r_x = 1.d0
         r_y = 1.d0
      else
         r_x = x1/r
         r_y = x2/r
      endif
      !
      ampl = 1.d0
      ca = R_CORE/0.9d0
      !
      if (r .le. ca) then
         E = ampl*exp(-ZI*OMEGA*x3)
         dE(1) = 0.d0
         dE(2) = 0.d0
         dE(3) = -ZI*OMEGA*E
      else
         E = ampl*exp(-((r-ca)**2.d0))*exp(-ZI*OMEGA*x3)
         dE(1) = -2.d0*(r-ca)*r_x*E
         dE(2) = -2.d0*(r-ca)*r_y*E
         dE(3) = -ZI*OMEGA*E
      endif
!
!
!--------------- 20th prob -------------------------------------------------------
!..Birefringent fiber with LP01 (E_x) and LP01 (E_y)
   elseif (ISOL .eq. 20) then
!
!  ...LMA fiber
!  ...LP01 (signal) in dielectric waveguide, a = 0.9*sqrt(2), omega=2*pi/0.1064=59.0525
      select case(Fld)
!     ...signal field
         case(1)
            if (ICOMP_TS .eq. 1) then
!              with CORE_NX = 1.4512
!              with CLAD_NX = 1.4500
               ampl =  0.9d0
               if (CORE_NX .eq. 1.4512d0 .and. CLAD_NX .eq. 1.4500d0) then
                  k    = 85.6833d0
                  gamm =  1.53131d0
                  beta =  3.12978d0
               else
                  write(*,*) 'mfd_solutions. signal: unexpected CLAD_NY. stop.'
                  stop
               endif
            else if (ICOMP_TS .eq. 2) then
!              with CORE_NY = 1.6510
!              with CLAD_NY = 1.6500
               ampl =  0.3d0
               if (CORE_NY .eq. 1.4512d0 .and. CLAD_NY .eq. 1.4500d0) then
                  k    = 85.6833d0
                  gamm =  1.53131d0
                  beta =  3.12978d0
               elseif (CORE_NY .eq. 1.6510d0 .and. CLAD_NY .eq. 1.6500d0) then
                  k    = 97.4838d0
                  gamm =  1.52302d0
                  beta =  3.03177d0
               else
                  write(*,*) 'mfd_solutions. signal: unexpected CLAD_NY. stop.'
                  stop
               endif
            endif
            call get_LP01(Xp,ampl,k,gamm,beta, E,dE)
!
!     ...pump field (LP01)
         case(0)
            if (ICOMP_TS .eq. 1) then
!              with CORE_NX = 1.4512
!              with CLAD_NX = 1.4500
               ampl =  2.0d0
               k    = 93.4108d0
               gamm =  1.55709d0
               beta =  3.46466d0
            else if (ICOMP_TS .eq. 2) then
!              with CORE_NY = 1.6510
!              with CLAD_NY = 1.6500
               ampl =  0.3d0
               if (CORE_NY .eq. 1.4512d0 .and. CLAD_NY .eq. 1.4500d0) then
                  k    = 93.4108d0
                  gamm =  1.55709d0
                  beta =  3.46466d0
               elseif (CLAD_NY .eq. 1.6500d0) then
                  k    = 106.275d0
                  gamm =  1.54933d0
                  beta =  3.35860d0
               else
                  write(*,*) 'mfd_solutions. pump: unexpected CLAD_NY. stop.'
                  stop
               endif
            endif
            call get_LP01(Xp,ampl,k,gamm,beta, E,dE)
!
         case default
            write(*,*) 'mfd_solutions: Fld_flag invalid. stop.'
            stop
      end select
!
!--------------- 21st prob -------------------------------------------------------
!..Birefringent fiber with LP01 (E_x) and LP11 (E_y)
   elseif (ISOL .eq. 21) then
!
!  ...LMA fiber
!  ...a = 0.9*sqrt(2)
!  ...omega_s=2*pi/0.1064=59.0525
!  ...omega_p=2*pi/0.0976
!
      if (R_CORE .ne. 0.9d0*sqrt(2.0d0)      .or.  &
          LAMBDA_SIGNAL .ne. 1064.0d-9/L_0   .or.  &
          LAMBDA_PUMP   .ne.  976.0d-9/L_0 ) then
         write(*,*) 'mfd_solutions: unexpected parameters. stop.'
         stop
      endif
!
      select case(Fld)
!     ...signal field
         case(1)
            ! compute signal power oscillation (1% up/down)
            ! ..if 'modified' activated
            cc = 1.0d0
            modified = 0
            if (modified .eq. 1) then
               ca = sqrt(0.99d0)
               cb = sqrt(1.01d0)
               if (TIMESTEP .le. 100) then
                  cc = 1.0d0
               else
                  ! oscillate at f = 1/(4*DELTA_T)
                  ! if DELTA_T=0.1ms, then f = 2500Hz
                  if (MOD(TIMESTEP,4) .eq. 1) then
                     cc = cb
                  else if (MOD(TIMESTEP,4) .eq. 3) then
                     cc = ca
                  else
                     cc = 1.0d0
                  endif
               endif
            endif
            ! ampl 0.9/0.6 -> ~90% power in LP01, ~10% power in LP11
            if (ICOMP_TS .eq. 1) then
               if (CORE_NX .eq. 1.4512d0 .and. CLAD_NX .eq. 1.4500d0) then
                  ampl =  0.9d0 * cc ! oscillating signal power
                  k    = 85.6833d0
                  gamm =  1.53131d0
                  beta =  3.12978d0
               else
                  write(*,*) 'mfd_solutions: unexpected CORE_NX,CLAD_NX. stop.'
                  stop
               endif
               call get_LP01(Xp,ampl,k,gamm,beta, E,dE)
            else if (ICOMP_TS .eq. 2) then
               if (CORE_NY .eq. 1.6510d0 .and. CLAD_NY .eq. 1.6500d0) then
                  ampl =  0.6d0 * cc ! oscillating signal power
                  k    = 97.4662d0
                  gamm =  2.39795d0
                  beta =  2.40022d0
               else
                  write(*,*) 'mfd_solutions: unexpected CORE_NY,CLAD_NY. stop.'
                  stop
               endif
               call get_LP11a(Xp,ampl,k,gamm,beta, E,dE)
               !call get_LP11b(Xp,ampl,k,gamm,beta, E,dE)
            endif
!
!     ...pump field LP01 (E_x)
         case(0)
            ! compute increasing pump power
            ! ..if 'modified' activated
            cc = 1.0d0
            modified = 1
            if (modified .eq. 1) then
               ca = 1.0d0/sqrt(2.0d0)
               cb = sqrt(2.0d0)
               if (TIMESTEP .le. 100) then
                  cc = ca
               else if (TIMESTEP .le. 150) then
                  cc = ca + (cb-ca)*(TIMESTEP-100.d0)/50.d0
               else
                  cc = cb
               endif
            endif
!
            ! ampl 2.0/3.0 -> ~66% power in LP01, ~33% in LP11
            if (ICOMP_TS .eq. 1) then
               if (CORE_NX .eq. 1.4512d0 .and. CLAD_NX .eq. 1.4500d0) then
                  ampl =  2.5d0 * cc ! increasing pump power
                  k    = 93.4108d0
                  gamm =  1.55709d0
                  beta =  3.46466d0
               else
                  write(*,*) 'mfd_solutions: unexpected CORE_NX,CLAD_NX. stop.'
                  stop
               endif
               call get_LP01(Xp,ampl,k,gamm,beta, E,dE)
            else if (ICOMP_TS .eq. 2) then
               if (CORE_NY .eq. 1.6510d0 .and. CLAD_NY .eq. 1.6500d0) then
                  ampl =   3.0d0 * cc ! increasing pump power
                  k    = 106.258d0
                  gamm =   2.44593d0
                  beta =   2.77454d0
               else
                  write(*,*) 'mfd_solutions. unexpected CLAD_NY. stop.'
                  stop
               endif
               call get_LP11a(Xp,ampl,k,gamm,beta, E,dE)
               !call get_LP11b(Xp,ampl,k,gamm,beta, E,dE)
            endif
         case default
            write(*,*) 'mfd_solutions: Fld_flag invalid. stop.'
            stop
      end select
!..endif ISOL
   endif
!
!
!---------------------------------------------------------------------------------
!
!.... set values as tensor products for all cases except
!.... Gaussian pulse and fiber beam
   select case(ISOL)
      case(0,1,2,3,4,9,10,11)
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
      case(5,8,12,13,14,140,15,150,16,17,18,19,20,21)
!     ...already computed fields
!
      case default
         write(*,*) 'mfd_solutions: invalid ISOL param. stop.'
         stop
   end select
end subroutine mfd_solutions

function BESSEL_dJ1(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   fval = BESSEL_J0(x) - BESSEL_J1(x)/x
end function

function BESSEL_J2(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   fval = BESSEL_JN(2, x)
end function

function BESSEL_dJ2(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessJY(x,2.d0, a,b,fval,c)
end function

function BESSEL_K0(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,0.d0, a,fval,b,c)
end function

function BESSEL_dK0(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,0.d0, a,b,c,fval)
end function

function BESSEL_K1(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,1.d0, a,fval,b,c)
end function

function BESSEL_dK1(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,1.d0, a,b,c,fval)
end function

function BESSEL_K2(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,2.d0, a,fval,b,c)
end function

function BESSEL_dK2(x) result(fval)
   real(8), intent(in) :: x
   real(8) :: fval
   real(8) :: a,b,c
   call dbessIK(x,2.d0, a,b,c,fval)
end function
!
!------------------------------------------------------
! subroutine get_LP01
!------------------------------------------------------
subroutine get_LP01(Xp,ampl,k,gamm,beta, E,dE)
!
   use commonParam
   use laserParam
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: ampl, k, gamm, beta
   VTYPE  , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, r_x, r_y, ca, cb
   real(8) :: BESSEL_K0, BESSEL_K1
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
!..radial coordinate
   r = sqrt(x1*x1+x2*x2)
!
   if (r .eq. 0.d0) then
      r_x = 1.d0
      r_y = 1.d0
   else
      r_x = x1/r
      r_y = x2/r
   endif
!
   if (r .le. R_CORE) then
      ca = ampl/BESSEL_J0(gamm*R_CORE)
      E = ca*BESSEL_J0(gamm*r)
      cb = -ca*gamm*BESSEL_J1(gamm*r)
   else
      ca = ampl/BESSEL_K0(beta*R_CORE)
      E = ca*BESSEL_K0(beta*r)
      cb = -ca*beta*BESSEL_K1(beta*r)
   endif
   dE(1) = cb*r_x*exp(-ZI*k*x3)
   dE(2) = cb*r_y*exp(-ZI*k*x3)
   E = E*exp(-ZI*k*x3)
   dE(3) = -ZI*k*E
!
end subroutine get_LP01
!
!------------------------------------------------------
! subroutine get_LP11a
!------------------------------------------------------
subroutine get_LP11a(Xp,ampl,k,gamm,beta, E,dE)
!
   use commonParam
   use laserParam
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: ampl, k, gamm, beta
   VTYPE  , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, ca, cb, cc
   real(8) :: BESSEL_dJ1, BESSEL_K1, BESSEL_dK1
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
!..shift source away from zero
   if (abs(x1) .lt. GEOM_TOL) then
      x1 = x1+GEOM_TOL
   endif
   if (abs(x2) .lt. GEOM_TOL) then
      x2 = x2+GEOM_TOL
   endif
   r = sqrt(x1*x1+x2*x2)
!
   if (r .le. R_CORE) then
      ca = ampl/BESSEL_J1(gamm*R_CORE)
      E = ca*(x1/r)*BESSEL_J1(gamm*r)
      cb = ca*(((x1/r)**(2.d0))*gamm*BESSEL_dJ1(gamm*r)+ &
               ((x2/r)**(2.d0))*BESSEL_J1(gamm*r)/r)
      cc = ca*(x2/r)*(x1/r)*(gamm*BESSEL_dJ1(gamm*r)-BESSEL_J1(gamm*r)/r)
   else
      ca = ampl/BESSEL_K1(beta*R_CORE)
      E = ca*(x1/r)*BESSEL_K1(beta*r)
      cb = ca*(((x1/r)**(2.d0))*beta*BESSEL_dK1(beta*r)+ &
               ((x2/r)**(2.d0))*BESSEL_K1(beta*r)/r)
      cc = ca*(x2/r)*(x1/r)*(beta*BESSEL_dK1(beta*r)-BESSEL_K1(beta*r)/r)
   endif
!
   dE(1) = cb*exp(-ZI*k*x3)
   dE(2) = cc*exp(-ZI*k*x3)
   E = E*exp(-ZI*k*x3)
   dE(3) = -ZI*k*E
!
   if (R_CLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP11a
!
!------------------------------------------------------
! subroutine get_LP11b (rotated LP11 mode)
!------------------------------------------------------
subroutine get_LP11b(Xp,ampl,k,gamm,beta, E,dE)
!
   use commonParam
   use laserParam
   use control, only : GEOM_TOL
!
   implicit none
!
   real*8, intent(in)  :: Xp(3)
   real*8, intent(in)  :: ampl, k, gamm, beta
   VTYPE , intent(out) :: E, dE(3)
!
   real*8 :: x1, x2, x3, r, ca, cb, cc
   real*8 :: BESSEL_dJ1, BESSEL_K1, BESSEL_dK1
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
!..shift source away from zero
   if (abs(x1) .lt. GEOM_TOL) then
      x1 = x1+GEOM_TOL
   endif
   if (abs(x2) .lt. GEOM_TOL) then
      x2 = x2+GEOM_TOL
   endif
   r = sqrt(x1*x1+x2*x2)
!
   if (r .le. R_CORE) then
      ca = ampl/BESSEL_J1(gamm*R_CORE)
      E = ca*(x2/r)*BESSEL_J1(gamm*r)
      cb = ca*(x2/r)*(x1/r)*(gamm*BESSEL_dJ1(gamm*r)-BESSEL_J1(gamm*r)/r)
      cc = ca*(((x2/r)**(2.d0))*gamm*BESSEL_dJ1(gamm*r) + &
               ((x1/r)**(2.d0))*BESSEL_J1(gamm*r)/r)
   else
      ca = ampl/BESSEL_K1(beta*R_CORE)
      E = ca*(x2/r)*BESSEL_K1(beta*r)
      cb = ca*(x2/r)*(x1/r)*(beta*BESSEL_dK1(beta*r)-BESSEL_K1(beta*r)/r)
      cc = ca*(((x2/r)**(2.d0))*beta*BESSEL_dK1(beta*r) + &
               ((x1/r)**(2.d0))*BESSEL_K1(beta*r)/r)
   endif
!
   dE(1) = cb*exp(-ZI*k*x3)
   dE(2) = cc*exp(-ZI*k*x3)
   E = E*exp(-ZI*k*x3)
   dE(3) = -ZI*k*E
!
   if (R_CLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP11b
!
!------------------------------------------------------
! subroutine get_LP02
!------------------------------------------------------
subroutine get_LP02(Xp,ampl,k,gamm,beta, E,dE)
!
   use commonParam
   use laserParam
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: ampl, k, gamm, beta
   VTYPE  , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, r_x, r_y, ca, cb
   real(8) :: BESSEL_K0, BESSEL_K1
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
!..radial coordinate
   r = sqrt(x1*x1+x2*x2)
!
   if (r .eq. 0.d0) then
      r_x = 1.d0
      r_y = 1.d0
   else
      r_x = x1/r
      r_y = x2/r
   endif
!
   if (r .le. R_CORE) then
      ca = ampl/BESSEL_J0(gamm*R_CORE)
      E = ca*BESSEL_J0(gamm*r)
      cb = -ca*gamm*BESSEL_J1(gamm*r)
   else
      ca = ampl/BESSEL_K0(beta*R_CORE)
      E = ca*BESSEL_K0(beta*r)
      cb = -ca*beta*BESSEL_K1(beta*r)
   endif
   dE(1) = cb*r_x*exp(-ZI*k*x3)
   dE(2) = cb*r_y*exp(-ZI*k*x3)
   E = E*exp(-ZI*k*x3)
   dE(3) = -ZI*k*E
!
   if (R_CLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP02
!
!------------------------------------------------------
! subroutine get_LP21a
!------------------------------------------------------
subroutine get_LP21a(Xp,ampl,k,gamm,beta, E,dE)
!
   use commonParam
   use laserParam
   use control, only : GEOM_TOL
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   real(8), intent(in)  :: ampl, k, gamm, beta
   VTYPE  , intent(out) :: E, dE(3)
!
   real(8) :: x1, x2, x3, r, r_x, r_y, ca, cb, cc
   real(8) :: BESSEL_J2, BESSEL_dJ2, BESSEL_K2, BESSEL_dK2
!
   real(8) :: cos_t,cos_2t
   real(8) :: sin_t,sin_2t
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
!..shift source away from zero
   if (abs(x1) .lt. GEOM_TOL) then
      x1 = x1+GEOM_TOL
   endif
   if (abs(x2) .lt. GEOM_TOL) then
      x2 = x2+GEOM_TOL
   endif
   r = sqrt(x1*x1+x2*x2)
!
   cos_t  = x1/r
   sin_t  = x2/r
   cos_2t = cos_t**(2.d0) - sin_t**(2.d0)
   sin_2t = 2 * sin_t * cos_t
!
   if (r .le. R_CORE) then
      ca = ampl/BESSEL_J2(gamm*R_CORE)
      E = ca*cos_2t*BESSEL_J2(gamm*r)
      cb = ca*(cos_t*cos_2t*gamm*BESSEL_dJ2(gamm*r) + &
               2.d0*sin_t*sin_2t*BESSEL_J2(gamm*r)/r)
      cc = ca*(sin_t*cos_2t*gamm*BESSEL_dJ2(gamm*r) - &
               2.d0*sin_2t*cos_t*BESSEL_J2(gamm*r)/r)
   else
      ca = ampl/BESSEL_K2(beta*R_CORE)
      E = ca*cos_2t*BESSEL_K2(beta*r)
      cb = ca*(cos_t*cos_2t*beta*BESSEL_dK2(beta*r) + &
               2.d0*sin_t*sin_2t*BESSEL_K2(beta*r)/r)
      cc = ca*(sin_t*cos_2t*beta*BESSEL_dK2(beta*r) - &
      2.d0*sin_2t*cos_t*BESSEL_K2(beta*r)/r)
   endif
!
   dE(1) = cb*exp(-ZI*k*x3)
   dE(2) = cc*exp(-ZI*k*x3)
   E = E*exp(-ZI*k*x3)
   dE(3) = -ZI*k*E
!
   if (R_CLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP21a

!------------------------------------------------------
! subroutine get_LP21b (rotated LP21 mode)
!------------------------------------------------------
subroutine get_LP21b(Xp,ampl,k,gamm,beta, E,dE)
!
   use commonParam
   use laserParam
   use control, only : GEOM_TOL
!
   implicit none
!
   real*8, intent(in)  :: Xp(3)
   real*8, intent(in)  :: ampl, k, gamm, beta
   VTYPE , intent(out) :: E, dE(3)
!
   real*8 :: x1, x2, x3, r, r_x, r_y, ca, cb, cc
   real*8 :: BESSEL_J2, BESSEL_dJ2, BESSEL_K2, BESSEL_dK2
!
   real*8 :: cos_t,cos_2t
   real*8 :: sin_t,sin_2t
!
!------------------------------------------------------
!
!..Cartesian coordinates
   x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)
!
!..shift source away from zero
   if (abs(x1) .lt. GEOM_TOL) then
      x1 = x1+GEOM_TOL
   endif
   if (abs(x2) .lt. GEOM_TOL) then
      x2 = x2+GEOM_TOL
   endif
   r = sqrt(x1*x1+x2*x2)
!
   cos_t  = x1/r
   sin_t  = x2/r
   cos_2t = cos_t**(2.d0) - sin_t**(2.d0)
   sin_2t = 2 * sin_t * cos_t
!
   if (r .le. R_CORE) then
      ca = ampl/BESSEL_J2(gamm*R_CORE)
      E = ca*sin_2t*BESSEL_J2(gamm*r)
      cb = ca*(cos_t*sin_2t*gamm*BESSEL_dJ2(gamm*r) - &
               2.d0*sin_t*cos_2t*BESSEL_J2(gamm*r)/r)
      cc = ca*(sin_t*sin_2t*gamm*BESSEL_dJ2(gamm*r) + &
               2.d0*cos_2t*cos_t*BESSEL_J2(gamm*r)/r)
   else
      ca = ampl/BESSEL_K2(beta*R_CORE)
      E = ca*sin_2t*BESSEL_K2(beta*r)
      cb = ca*(cos_t*cos_2t*beta*BESSEL_dK2(beta*r) - &
               2.d0*sin_t*cos_2t*BESSEL_K2(beta*r)/r)
      cc = ca*(sin_t*sin_2t*beta*BESSEL_dK2(beta*r) + &
               2.d0*cos_2t*cos_t*BESSEL_K2(beta*r)/r)
   endif
!
   dE(1) = cb*exp(-ZI*k*x3)
   dE(2) = cc*exp(-ZI*k*x3)
   E = E*exp(-ZI*k*x3)
   dE(3) = -ZI*k*E
!
   if (R_CLAD - r < GEOM_TOL) then
      E = 0.d0
      dE(1:3) = 0.d0
   endif
!
end subroutine get_LP21b
