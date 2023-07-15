!-----------------------------------------------------------------------------------
!> @brief      Compute relevant quantities for solutions in elastiticy problems
!!
!> @param[in]  X        - a point in physical space
!> @param[out] U        - value of the solution displacement
!> @param[out] GradU    - corresponding first derivatives
!> @param[out] Epsilon  - value of the solution strain
!> @param[out] Sigma    - value of the solution stress
!> @param[out] DivSigma - divergence of the solution stress
!!
!> @date       July 2023
!-----------------------------------------------------------------------------------
   subroutine elast_solution(X, U,GradU,Epsilon,Sigma,DivSigma)
!
      use data_structure3D
      use common_prob_data, only : IEXACT_PROB,PI,NP1,NP2,NP3
      use isotropic_elast_material
!      use transverse_isotropic_elast_material
!
      implicit none
!
      real(8), intent(in)  :: X(3)
      real(8), intent(out) :: U(3)
      real(8), intent(out) :: GradU(3,3)
      real(8), intent(out) :: Epsilon(3,3)
      real(8), intent(out) :: Sigma(3,3)
      real(8), intent(out) :: DivSigma(3)
!
!  ...local variables
      character(len=6) :: prob
      integer :: itest, icomp
      integer :: i, j, k, l, m
      real(8) :: x1, x2, x3
      real(8) :: tol, diffmax, dmax
      real(8) :: TmpTensor(3,3)
      real(8) :: Grad2U(3,3,3), GradSigma(3,3,3)
      real(8) :: C(3,3,3,3), A(3,3,3,3)
!
      integer :: iprint = 0
!
!------------------------------------------------------------------------------
!
      itest  = 1
      tol = 1.0d-14

!  ...get stiffness and compliance tensors
      call getC(X, C)
      call getA(X, A)

!  ...initialize variables
      U = 0.d0; GradU = 0.d0; Epsilon = 0.d0; Sigma = 0.d0; DivSigma = 0.d0
      Grad2U = 0.d0; GradSigma = 0.d0
!
!-----------------------------------------------------------------------------------
!      D E C L A R E    S O L U T I O N    V A R I A B L E S                       |
!-----------------------------------------------------------------------------------
!
!  NOTES:
!
!     - The first 4 solutions are calculated from explicit formulas for the
!       displacement variable.
!
!     - The first singular solution (case 5) has all values computed independently.
!       This is the plane strain extended L-shaped domain problem.
!
!     - The next 4 solutions are calculated from explicit formulas for the stress.
!       The displacement variable is unknown in this case so only the residual would
!       be able to be used to determine convergence rates. However these solutions can
!       be used to ensure that the stress boundary conditions are being imposed
!       properly
!
!     - The final solution (case 10) is a singular solution with unknown displacement.
!       This is the plane stress extended L-shaped domain problem. Again, only the
!       residual can be used to accurately calculate convergence rates.
!
!-----------------------------------------------------------------------------------
!
      x1 = X(1); x2 = X(2); x3 = X(3)
!
      select case(IEXACT_PROB)
!  ...polynomial solution (0/4)
      case(0)
!     ...displacement
         U(1) = x2
         U(2) = x3
         U(3) = x1
!
!     ...1st derivatives
         GradU(1,2) = 1.d0
         GradU(2,3) = 1.d0
         GradU(3,1) = 1.d0
!
!  ...polynomial solution (2/4)
      case(1)
!     ...displacement
         U(1:3) = x1*x2*x3
!
!     ...1st derivatives
         GradU(1:3,1) = x2*x3
         GradU(1:3,2) = x1*x3
         GradU(1:3,3) = x1*x2
!
!     ...2nd derivatives
         Grad2U(1:3,1,2) = x3; Grad2U(1:3,1,3) = x2; Grad2U(1:3,2,3) = x1
         Grad2U(1:3,2,1) = Grad2U(1:3,1,2)
         Grad2U(1:3,3,1) = Grad2U(1:3,1,3)
         Grad2U(1:3,3,2) = Grad2U(1:3,2,3)
!
!  ...polynomial solution (2/4)
      case(2)
!     ...displacement
         U(1) = x1**NP1
         U(2) = x2**NP2
         U(3) = x3**NP3
!
!     ...1st derivatives
         GradU(1,1) = NP1*x1**(NP1-1)
         GradU(2,2) = NP2*x2**(NP2-1)
         GradU(3,3) = NP3*x3**(NP3-1)
!
!     ...2nd derivatives
         Grad2U(1,1,1) = NP1*(NP1-1)*x1**(NP1-2)
         Grad2U(2,2,2) = NP2*(NP2-1)*x2**(NP2-2)
         Grad2U(3,3,3) = NP3*(NP3-1)*x3**(NP3-2)
!
!  ...exponential solution (1/2)
      case(3)
!     ...displacement
         U(1) = 1.d0
         U(2) = 1.d0
         U(3) = dexp(x3)
!
!     ...1st derivatives
         GradU(3,3) = dexp(x3)
!
!     ...2nd derivatives
         Grad2U(3,3,3) = dexp(x3)
!
!  ...sinusoidal solution (1/2)
      case(4)
!     ...displacement
         U(1:3) = dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
!
!     ...1st derivatives
         GradU(1:3,1) = PI*dcos(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
         GradU(1:3,2) = PI*dsin(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
         GradU(1:3,3) = PI*dsin(PI*x1)*dsin(PI*x2)*dcos(PI*x3)
!
!     ...2nd derivatives
         Grad2U(1:3,1,1) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
         Grad2U(1:3,1,2) = PI**2*dcos(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
         Grad2U(1:3,1,3) = PI**2*dcos(PI*x1)*dsin(PI*x2)*dcos(PI*x3)
         Grad2U(1:3,2,1) = Grad2U(1:3,1,2)
         Grad2U(1:3,2,2) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
         Grad2U(1:3,2,3) = PI**2*dsin(PI*x1)*dcos(PI*x2)*dcos(PI*x3)
         Grad2U(1:3,3,1) = Grad2U(1:3,1,3)
         Grad2U(1:3,3,2) = Grad2U(1:3,2,3)
         Grad2U(1:3,3,3) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
!
!  ...singular solution (1/2)
!     NOTE:: use NEXACT=2 and L-shape geometry
      case(5)
         prob = 'strain'
         call singular_solution(X(1:2),prob, U,GradU,Epsilon,Sigma)
!
!  ...polynomial solution (3/4)
      case(6)
!     ...stress
         Sigma(1:3,1:3) = x1*x2*x3
!
!     ...1st derivatives
         GradSigma(1:3,1:3,1) = x2*x3
         GradSigma(1:3,1:3,2) = x1*x3
         GradSigma(1:3,1:3,3) = x1*x2
!
!  ...polynomial solution (4/4)
      case(7)
!     ...stress
         Sigma(1,1) = x1**NP1
         Sigma(2,2) = x2**NP2
         Sigma(3,3) = x3**NP3
!
!     ...1st derivatives
         GradSigma(1,1,1) = NP1*x1**(NP1-1)
         GradSigma(2,2,2) = NP2*x2**(NP2-1)
         GradSigma(3,3,3) = NP3*x3**(NP3-1)
!
!  ...exponential solution (2/2)
      case(8)
!     ...stress
         Sigma(1:3,1:3) = 1.d0
         Sigma(3,3) = dexp(x3)
!
!     ...1st derivatives
         GradSigma(3,3,3) = dexp(x3)
!
!  ...sinusoidal solution (2/2)
      case(9)
!     ...stress
         Sigma(1:3,1:3) = dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
!
!     ...1st order derivatives
         GradSigma(1:3,1:3,1) = PI*dcos(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
         GradSigma(1:3,1:3,2) = PI*dsin(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
         GradSigma(1:3,1:3,3) = PI*dsin(PI*x1)*dsin(PI*x2)*dcos(PI*x3)
!
!  ...singular solution (2/2)
!     NOTE: use NEXACT=2 and L-shape geometry
      case(10)
         prob = 'stress'
         call singular_solution(X(1:2),prob, U,GradU,Epsilon,Sigma)
      end select
!
!-----------------------------------------------------------------------------------
!      C A L C U L A T E    D E R I V E D    V A R I A B L E S                     |
!-----------------------------------------------------------------------------------
!
      select case(IEXACT_PROB)
!
!  ...Here the strain and stress have to be calculated
      case(0,1,2,3,4)
         do icomp = 1,3
            do k = 1,3
!           ...calculate strain
               Epsilon(icomp,k) = 0.5d0*(GradU(icomp,k)+GradU(k,icomp))
               do l = 1,3
!              ...calculate Cauchy stress
                  Sigma(icomp,1:3) = Sigma(icomp,1:3) + C(icomp,1:3,k,l)*GradU(k,l)
                  do m = 1,3
                     GradSigma(icomp,1:3,m) = GradSigma(icomp,1:3,m) + C(icomp,1:3,k,l)*Grad2U(k,l,m)
                  enddo
               enddo
            enddo
!        ...calculate divergence of the stress
            do m = 1,3
               DivSigma(icomp) = DivSigma(icomp) + GradSigma(icomp,m,m)
            enddo
         enddo
!
!  ...Here the strain and body forces have to be calculated
      case(6,7,8,9)
         do icomp = 1,3
            do k = 1,3
!           ...calculate divergence of the stress
               DivSigma(icomp) = DivSigma(icomp) + GradSigma(icomp,k,k)
               do l = 1,3
!              ...calculate strain
                  Epsilon(icomp,1:3) = Epsilon(icomp,1:3) + A(icomp,1:3,k,l)*Sigma(k,l)
               enddo
            enddo
         enddo
      end select
!
!-----------------------------------------------------------------------------------
!     T E S T S    A N D    P R I N T    S T A T E M E N T S                       |
!-----------------------------------------------------------------------------------
!
#if DEBUG_MODE
   if (itest.eq.1) then
!
!   SYMMETRY
!
      diffmax = 0.d0; dmax = 0.d0
      do i=1,3
         do j=i,3
            diffmax = max(diffmax,abs(Epsilon(i,j)-Epsilon(j,i)))
            dmax = max(dmax,abs(Epsilon(i,j)))
         enddo
      enddo
      if (diffmax/dmax.gt.tol) then
         write(*,7000) diffmax, dmax
 7000    format('elast_solution: diffmax, dmax FOR Epsilon = ',2e12.5)
         call pause
      endif
!
      diffmax = 0.d0; dmax = 0.d0
      do i=1,3
         do j=i,3
            diffmax = max(diffmax,abs(Sigma(i,j)-Sigma(j,i)))
            dmax = max(dmax,abs(Sigma(i,j)))
         enddo
      enddo
      if (diffmax/dmax.gt.tol) then
         write(*,7001) diffmax, dmax
 7001    format('elast_solution: diffmax, dmax FOR Sigma = ',2e12.5)
         call pause
      endif
!
!   C:grad(u) = sigma
!
      if (IEXACT_PROB.le.5) then
         tmpTensor = 0.d0
         do icomp=1,3
            do k=1,3
               do l=1,3
                  tmpTensor(icomp,1:3) = tmpTensor(icomp,1:3) + C(icomp,1:3,k,l)*GradU(k,l)
               enddo
            enddo
         enddo
         diffmax = 0.d0; dmax = 0.d0
         do i=1,3
            do j=i,3
               diffmax = max(diffmax,abs(tmpTensor(i,j)-Sigma(i,j)))
               dmax = max(dmax,abs(Sigma(i,j)))
            enddo
         enddo
         if (diffmax/dmax.gt.tol) then
            write(*,7002) diffmax, dmax
   7002     format('elast_solution: diffmax, dmax FOR C:grad(u) - Sigma = ',2e12.5)
            call pause
         endif
      endif
!
!   C:epsilon = sigma
!
      tmpTensor = 0.d0
      do icomp=1,3
         do k=1,3
            do l=1,3
               tmpTensor(icomp,1:3) = tmpTensor(icomp,1:3) + C(icomp,1:3,k,l)*Epsilon(k,l)
            enddo
         enddo
      enddo
      diffmax = 0.d0; dmax = 0.d0
      do i=1,3
         do j=i,3
            diffmax = max(diffmax,abs(tmpTensor(i,j)-Sigma(i,j)))
            dmax = max(dmax,abs(Sigma(i,j)))
         enddo
      enddo
      if (diffmax/dmax.gt.tol) then
         write(*,7003) diffmax, dmax
 7003    format('elast_solution: diffmax, dmax FOR C:epsilon - Sigma = ',2e12.5)
         call pause
      endif
!
!   A:sigma = epsilon
!
      tmpTensor = 0.d0
      do icomp=1,3
         do k=1,3
            do l=1,3
               tmpTensor(icomp,1:3) = tmpTensor(icomp,1:3) + A(icomp,1:3,k,l)*Sigma(k,l)
            enddo
         enddo
      enddo
!
      diffmax = 0.d0; dmax = 0.d0
      do i=1,3
         do j=i,3
            diffmax = max(diffmax,abs(tmpTensor(i,j)-Epsilon(i,j)))
            dmax = max(dmax,abs(Sigma(i,j)))
         enddo
      enddo
!
      if (diffmax/dmax.gt.tol) then
         write(*,7004) diffmax, dmax
 7004    format('elast_solution: diffmax, dmax FOR A:Sigma - Epsilon = ',2e12.5)
         call pause
      endif
   endif
#endif
!
   end subroutine elast_solution





!-----------------------------------------------------------------------------------
!> @brief      Evaluate a manufactured solution for the L-shape linear elasticity
!              problem
!!
!> @param[in]  X       - a point in physical space
!> @param[out] U       - value of the solution displacement
!> @param[out] GradU   - corresponding first derivatives
!> @param[out] Epsilon - value of the solution strain
!> @param[out] Sigma   - value of the solution stress
!!
!> @date       July 2023
!-----------------------------------------------------------------------------------
   subroutine singular_solution(Xtmp,prob, U,GradU,Epsilon,Sigma)
!
      use common_prob_data, only : PI, EPS
      use isotropic_elast_material, only : LAMBDA,MU
!
      implicit none
!
      character(len=6), intent(in)  :: prob
      real(8),          intent(in)  :: Xtmp(2)
      real(8),          intent(out) :: U(3)
      real(8),          intent(out) :: GradU(3,3)
      real(8),          intent(out) :: Epsilon(3,3)
      real(8),          intent(out) :: Sigma(3,3)
!
  
      real(8) :: x, y, r, drdx, drdy, poisson, sig
      real(8) :: theta, dthetadx,dthetady
      real(8) :: a, b, alpha, ak,
      real(8) :: F, Fd, Fdd, G, Gd, Gdd
      real(8) :: ur, utheta, urdr, urdtheta, uthetadr, uthetadtheta
      real(8) :: trr, tthetatheta, trtheta
      real(8) :: tmp1, tmp2, v11tmp, v12tmp, v22tmp
      real(8) :: u11tmp, u12tmp
      real(8) :: u21tmp, u22tmp
!
      real*8, parameter :: inf = 1.0d32
!
      integer :: iprint = 0
!
!-----------------------------------------------------------------------------------
!
!  ...initialize variables (these stay equal to 0 if r < EPS)
      U(:) = 0.d0
      GradU(:,:) = 0.d0
      Epsilon(:,:) = 0.d0
      Sigma(:,:) = 0.d0
!
! NOTE: this is always a solution of the homogeneous equation so DivSigma remains
!       equal to 0
!
!-----------------------------------------------------------------------------------
!      D E C L A R E    C O O R D I N A T E S    A N D    P A R A M E T E R S      |
!-----------------------------------------------------------------------------------
!
!  ...rotate coordinates 45 degrees counter-clockwise
      x = (Xtmp(1)+Xtmp(2))/dsqrt(2.d0)
      y = (-Xtmp(1)+Xtmp(2))/dsqrt(2.d0)
!
      r = dsqrt(x**2+y**2)
!
      if (r.lt.EPS) return
!
      theta = datan2(y,x)
!
      drdx = x/r ; dthetadx =-y/r**2
      drdy = y/r ; dthetady = x/r**2
!
      poisson = LAMBDA/(2.d0*(LAMBDA+MU))
      select case(prob)
      case('stress')
         sig = poisson/(1.d0+poisson)
         a = 0.6037781005214915d0
      case('strain')
         sig = poisson
         a = 0.5945633983288462d0
      case default
         write(*,*) 'singular_solution : NOT A VALID SOLUTION TYPE'
         stop 1
      end select

      alpha = 3.d0*PI/4.d0
      ak = (4.d0*(1-sig)-(a+1.d0))*dsin((a-1.d0)*alpha) &
            /((a+1.d0)*dsin((a+1.d0)*alpha))
!
!-----------------------------------------------------------------------------------
!      A S S E M B L E    S O L U T I O N                                          |
!-----------------------------------------------------------------------------------
!
      F   =  ak*dsin((a+1.d0)*theta) + dsin((a-1.d0)*theta)
      Fd  =  (a+1.d0)*ak*dcos((a+1.d0)*theta) &
           + (a-1.d0)*dcos((a-1.d0)*theta)
      Fdd = -(a+1.d0)**2*ak*dsin((a+1.d0)*theta) &
           - (a-1.d0)**2*dsin((a-1.d0)*theta)
      G   = -4.d0/(a-1.d0)*dcos((a-1.d0)*theta)
      Gd  =  4.d0*dsin((a-1.d0)*theta)
      Gdd =  4.d0*(a-1.d0)*dcos((a-1.d0)*theta)
!
      b      = 1.d0/(2.d0*MU)
      ur     = b*r**a*(-(a+1.d0)*F+(1-sig)*Gd)
      utheta = b*r**a*(-Fd+(1.d0-sig)*(a-1.d0)*G)
!
      urdr         = a*b*r**(a-1.d0)*(-(a+1.d0)*F+(1-sig)*Gd)
      urdtheta     = b*r**a*(-(a+1.d0)*Fd+(1-sig)*Gdd)
      uthetadr     = a*b*r**(a-1.d0)*(-Fd+(1.d0-sig)*(a-1.d0)*G)
      uthetadtheta = b*r**a*(-Fdd+(1.d0-sig)*(a-1.d0)*Gd)
!
      if (iprint.eq.1) then
         write(*,7001) ur,utheta,theta
7001     format('singular_solution : ur,utheta,theta = ',3e12.5)
      endif
!
!  ...trr,tthetatheat,trtheta are the cylindrical stress components
      trr         = r**(a-1.d0)*(Fdd+(a+1.d0)*F)
      tthetatheta = a*(a+1.d0)*r**(a-1.d0)*F
      trtheta     =-a*r**(a-1.d0)*Fd
!
!  ...displacements in the rotated system of coordinates
      tmp1 = ur*dcos(theta) - utheta*dsin(theta)
      tmp2 = ur*dsin(theta) + utheta*dcos(theta)
!
!  ...transform displacements to the rotated system of coordinates
      U(1) = (tmp1-tmp2)/dsqrt(2.d0)
      U(2) = (tmp1+tmp2)/dsqrt(2.d0)
!
!  ...and their derivatives in the rotated system of coordinates
      u11tmp = (urdr*dcos(theta)-uthetadr*dsin(theta))*drdx  &
             + (urdtheta*dcos(theta)-ur*dsin(theta)-uthetadtheta*dsin(theta)  &
             - utheta*dcos(theta))*dthetadx
      u12tmp = (urdr*dcos(theta)-uthetadr*dsin(theta))*drdy  &
             + (urdtheta*dcos(theta)-ur*dsin(theta)-uthetadtheta*dsin(theta)  &
             - utheta*dcos(theta))*dthetady
      u21tmp = (urdr*dsin(theta)+uthetadr*dcos(theta))*drdx  &
             + (urdtheta*dsin(theta)+ur*dcos(theta)+uthetadtheta*dcos(theta)  &
             - utheta*dsin(theta))*dthetadx
      u22tmp = (urdr*dsin(theta)+uthetadr*dcos(theta))*drdy  &
             + (urdtheta*dsin(theta)+ur*dcos(theta)+uthetadtheta*dcos(theta)  &
             - utheta*dsin(theta))*dthetady
!
!  ...transform displacement derivatives to the original system
      GradU(1,1) = 0.5d0*(u11tmp-u21tmp-u12tmp+u22tmp)
      GradU(1,2) = 0.5d0*(u11tmp-u21tmp+u12tmp-u22tmp)
      GradU(2,1) = 0.5d0*(u11tmp+u21tmp-u12tmp-u22tmp)
      GradU(2,2) = 0.5d0*(u11tmp+u21tmp+u12tmp+u22tmp)
!
      Epsilon(1,1) = GradU(1,1)
      Epsilon(2,2) = GradU(2,2)
      Epsilon(1,2) = 0.5*(GradU(1,2)+GradU(2,1)) ! may be improved
      Epsilon(2,1) = Epsilon(1,2)
!
!  ...transform stresses from the cylindrical to the rotated Cartesian
!     system of coordinates
      v11tmp = trr*(dcos(theta)**2)  &
               - 2.d0*dsin(theta)*dcos(theta)*trtheta  &
               + tthetatheta*(dsin(theta)**2)
      v22tmp = trr*(dsin(theta)**2)  &
               + 2.d0*dsin(theta)*dcos(theta)*trtheta  &
               + tthetatheta*(dcos(theta)**2)
      v12tmp = trr*dsin(theta)*dcos(theta)  &
               + (dcos(theta)**2-dsin(theta)**2)*trtheta  &
               - tthetatheta*dsin(theta)*dcos(theta)
!
!  ...transform stresses to the original system of coordinates
      Sigma(1,1) = 0.5d0*(v11tmp-2*v12tmp+v22tmp)
      Sigma(1,2) = 0.5d0*(v11tmp-v22tmp)
      Sigma(2,2) = 0.5d0*(v11tmp+2*v12tmp+v22tmp)
      Sigma(2,1) = Sigma(1,2)
!
      select case(prob)
      case('stress')
         U(3) = inf ! need formula for this component
         GradU(3,3) =-LAMBDA/(2.0d0*MU*(2.0d0*MU+3.0d0*LAMBDA))  &
                     *(Sigma(1,1)+Sigma(2,2))
         Epsilon(3,3) = GradU(3,3)
      case('strain')
         Sigma(3,3) = LAMBDA*(Epsilon(1,1)+Epsilon(2,2))
      case default
         write(*,*) 'singular_solution : NOT A VALID SOLUTION TYPE'
         stop 1
      end select
!
   end subroutine singular_solution
