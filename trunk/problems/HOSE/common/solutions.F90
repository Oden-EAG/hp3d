!-----------------------------------------------------------------------------------
!> Purpose : compute all relevant quantities of the exact solutions for the
!            linear elasticty problems
!!
!! @param[in]  X        - a point in physical space
!! @param[out] U        - value of the solution displacement
!! @param[out] GradU    - corresponding first derivatives
!! @param[out] Epsilon  - value of the solution strain
!! @param[out] Sigma    - value of the solution stress
!! @param[out] DivSigma - divergence of the solution stress
!-----------------------------------------------------------------------------------
!
subroutine elast_solution(X, U,GradU,Epsilon,Sigma,DivSigma)
  use data_structure3D
  use common_prob_data, only : IEXACT_PROB,PI,NP1,NP2,NP3
  use sheathed_isotropic_materials
  implicit none
!-----------------------------------------------------------------------------------
  real*8, dimension(3),   intent(in)  :: X
  real*8, dimension(3),   intent(out) :: U
  real*8, dimension(3,3), intent(out) :: GradU
  real*8, dimension(3,3), intent(out) :: Epsilon
  real*8, dimension(3,3), intent(out) :: Sigma
  real*8, dimension(3),   intent(out) :: DivSigma
!-----------------------------------------------------------------------------------
! local variables
  integer :: iprint,itest,icomp,i,j,k,l,m
  integer :: ndom
  real*8  :: x1,x2,x3,tol,diffmax,dmax,u_0
  real*8  :: r
  real*8, dimension(3,3)     :: TmpTensor
  real*8, dimension(3,3,3)   :: Grad2U,GradSigma
  real*8, dimension(3,3,3,3) :: C,A
! !------------------------------------------------------------------------------
! ! workspace for singular solution
!   real*8 :: SolV1(2),Soldiv1,SolV2(2),Soldiv2,SolQ1,SolQ2,SolQ3, &
!             Solu11,Solu12,Solu21,Solu22
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
  iprint = 1
  itest  = 0
  tol = 1.0d-14
!
! get stiffness and compliance tensors
  ndom=1
  r = dsqrt(X(2)**2+X(3)**2)
  if (r.lt.R_middle) ndom=2
  if (NU_rubber.eq.0.5d0) then
    !  Stiffness tensor is not defined as is.
    iprint = 0
    C=0.d0
  else
    call getC(X,ndom, C)
  endif    
  call getA(X,ndom, A)

! initialize variables
  U = 0.d0; GradU = 0.d0; Epsilon = 0.d0; Sigma = 0.d0; DivSigma = 0.d0
  Grad2U = 0.d0; GradSigma = 0.d0
!
!-----------------------------------------------------------------------------------
!      D E C L A R E    S O L U T I O N    V A R I A B L E S                       |
!-----------------------------------------------------------------------------------
!
!  NOTES::
!
!    The first 4 solutions are calculated from explicit formulas for the
!  displacement variable.
!    The first singular solution (case 5) has all values computed independently.
!  This is the plane strain extended L-shaped domain problem.
!    The next 4 solutions are calculated from explicit formulas for the stress.
!  The displacement variable is unknown in this case so only the residual would
!  be able to be used to determine convergence rates. However these solutions can
!  be used to ensure that the stress boundary conditions are being imposed
!  properly
!    The final solution (case 10) is a singular solution with unknown displacement.
!  This is the plane stress extended L-shaped domain problem. Again, only the
!  residual can be used to accurately calculate convergence rates.
!
!-----------------------------------------------------------------------------------
!
  x1 = X(1); x2 = X(2); x3 = X(3)
  u_0 = 1.0d-2
!
  select case(IEXACT_PROB)
! POLYNOMIAL SOLUTION (0/4)
   case(0)
!   displacement
    U(1) = u_0*x2
    U(2) = u_0*x3
    U(3) = u_0*x1
!
!   1st order derivatives
    GradU(1,2) = u_0
    GradU(2,3) = u_0
    GradU(3,1) = u_0
!
! POLYNOMIAL SOLUTION (1/4)
  case(1)
!   displacement
    U(1) = u_0*x2*x3
    U(2) = u_0*x1*x3
    U(3) = u_0*x1*x2
!
!   1st order derivatives
    GradU(1,1) = 0.d0;   GradU(1,2) = u_0*x3; GradU(1,3) = u_0*x2
    GradU(2,1) = u_0*x3; GradU(2,2) = 0.d0;   GradU(2,3) = u_0*x1
    GradU(3,1) = u_0*x2; GradU(3,2) = u_0*x1; GradU(3,3) = 0.d0
!
!   2nd order derivatives
    Grad2U(1,2,3) = u_0; Grad2U(1,3,2) = u_0;
    Grad2U(2,1,3) = u_0; Grad2U(2,3,1) = u_0;
    Grad2U(3,1,2) = u_0; Grad2U(3,2,1) = u_0;
!
! POLYNOMIAL SOLUTION (2/4)
  case(2)
!   displacement
    U(1) = x1**NP1
    U(2) = x2**NP2
    U(3) = x3**NP3
!
!   1st order derivatives
    GradU(1,1) = NP1*x1**(NP1-1)
    GradU(2,2) = NP2*x2**(NP2-1)
    GradU(3,3) = NP3*x3**(NP3-1)
!
!   2nd order derivatives
    Grad2U(1,1,1) = NP1*(NP1-1)*x1**(NP1-2)
    Grad2U(2,2,2) = NP2*(NP2-1)*x2**(NP2-2)
    Grad2U(3,3,3) = NP3*(NP3-1)*x3**(NP3-2)
!
! EXPONENTIAL SOLUTION (1/2)
  case(3)
!   displacement
    U(1) = 1.d0
    U(2) = 1.d0
    U(3) = dexp(x2)
!
!   1st order derivatives
    GradU(3,2) = dexp(x2)
!
!   2nd order derivatives
    Grad2U(3,2,2) = dexp(x2)
!
! SINUSOIDAL SOLUTION (1/2)
  case(4)
!   displacement
    U(1:3) = dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
!
!   1st order derivatives
    GradU(1:3,1) = PI*dcos(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
    GradU(1:3,2) = PI*dsin(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
    GradU(1:3,3) = PI*dsin(PI*x1)*dsin(PI*x2)*dcos(PI*x3)
!
!   2nd order derivatives
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
! Exact Solution for uniform pressure
  case(5)
    call exact_solution(X,ndom, U,GradU,Epsilon,Sigma)
!
! POLYNOMIAL SOLUTION (3/4)
  case(6)
!   stress
    Sigma(1:3,1:3) = x1*x2*x3
!
!   1st order derivatives
    GradSigma(1:3,1:3,1) = x2*x3
    GradSigma(1:3,1:3,2) = x1*x3
    GradSigma(1:3,1:3,3) = x1*x2
!
! POLYNOMIAL SOLUTION (4/4)
  case(7)
!   stress
    Sigma(1,1) = x1**NP1
    Sigma(2,2) = x2**NP2
    Sigma(3,3) = x3**NP3
!
!   1st order derivatives
    GradSigma(1,1,1) = NP1*x1**(NP1-1)
    GradSigma(2,2,2) = NP2*x2**(NP2-1)
    GradSigma(3,3,3) = NP3*x3**(NP3-1)
!
! EXPONENTIAL SOLUTION (2/2)
  case(8)
!   stress
    Sigma(1:3,1:3) = 1.d0
    Sigma(3,3) = dexp(x3)
!
!   1st order derivatives
    GradSigma(3,3,3) = dexp(x3)
!
! SINUSOIDAL SOLUTION (2/2)
  case(9)
!   stress
    Sigma(1:3,1:3) = dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
!
!   1st order derivatives
    GradSigma(1:3,1:3,1) = PI*dcos(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
    GradSigma(1:3,1:3,2) = PI*dsin(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
    GradSigma(1:3,1:3,3) = PI*dsin(PI*x1)*dsin(PI*x2)*dcos(PI*x3)
  end select
!
!-----------------------------------------------------------------------------------
!      C A L C U L A T E    D E R I V E D    V A R I A B L E S                     |
!-----------------------------------------------------------------------------------
!
  select case(IEXACT_PROB)
!
! Here the strain and stress have to be calculated
  case(0,1,2,3,4)
    do icomp=1,3
      do k=1,3
        !  calculate strain
        Epsilon(icomp,k) = 0.5d0*(GradU(icomp,k)+GradU(k,icomp))
        do l=1,3
          !  calculate Cauchy stress
          Sigma(icomp,1:3) = Sigma(icomp,1:3) + C(icomp,1:3,k,l)*GradU(k,l)
          do m=1,3
            GradSigma(icomp,1:3,m) = GradSigma(icomp,1:3,m) + C(icomp,1:3,k,l)*Grad2U(k,l,m)
          enddo
        enddo
      enddo
      !  calculate divergence of the stress
      do m=1,3
        DivSigma(icomp) = DivSigma(icomp) + GradSigma(icomp,m,m)
      enddo
    enddo
!
! Here the strain and body forces have to be calculated
  case(6,7,8,9)
    do icomp=1,3
      do k=1,3
        !  calculate divergence of the stress
        DivSigma(icomp) = DivSigma(icomp) + GradSigma(icomp,k,k)
        do l=1,3
          !  calculate strain
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
 7000 format('elast_solution: diffmax, dmax FOR Epsilon = ',2e12.5)
      call pause
    endif

    diffmax = 0.d0; dmax = 0.d0
    do i=1,3
      do j=i,3
        diffmax = max(diffmax,abs(Sigma(i,j)-Sigma(j,i)))
        dmax = max(dmax,abs(Sigma(i,j)))
      enddo
    enddo
    if (diffmax/dmax.gt.tol) then
      write(*,7001) diffmax, dmax
 7001 format('elast_solution: diffmax, dmax FOR Sigma = ',2e12.5)
      call pause
    endif
!
!   C:grad(u) = sigma
!
    if (IEXACT_PROB.le.5) then
      tmpTensor = 0.d0
      do icomp=1,3
        do k=1,3; do l=1,3
          tmpTensor(icomp,1:3) = tmpTensor(icomp,1:3) + C(icomp,1:3,k,l)*GradU(k,l)
        enddo; enddo
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
   7002 format('elast_solution: diffmax, dmax FOR C:grad(u) - Sigma = ',2e12.5)
        call pause
      endif
    endif
!
!   C:epsilon = sigma
!
    tmpTensor = 0.d0
    do icomp=1,3
      do k=1,3; do l=1,3
        tmpTensor(icomp,1:3) = tmpTensor(icomp,1:3) + C(icomp,1:3,k,l)*Epsilon(k,l)
      enddo; enddo
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
 7003 format('elast_solution: diffmax, dmax FOR C:epsilon - Sigma = ',2e12.5)
      call pause
    endif
!
!   A:sigma = epsilon
!
    tmpTensor = 0.d0
    do icomp=1,3
      do k=1,3; do l=1,3
        tmpTensor(icomp,1:3) = tmpTensor(icomp,1:3) + A(icomp,1:3,k,l)*Sigma(k,l)
      enddo; enddo
    enddo

    diffmax = 0.d0; dmax = 0.d0
    do i=1,3
      do j=i,3
        diffmax = max(diffmax,abs(tmpTensor(i,j)-Epsilon(i,j)))
        dmax = max(dmax,abs(Sigma(i,j)))
      enddo
    enddo
    if (diffmax/dmax.gt.tol) then
      write(*,7004) diffmax, dmax
 7004 format('elast_solution: diffmax, dmax FOR A:Sigma - Epsilon = ',2e12.5)
      call pause
    endif
  endif

end subroutine elast_solution

!-----------------------------------------------------------------------------------
!> Purpose : evaluate an exact solution for the sheathed materials problem
!!
!! @param[in]  X       - a point in physical space
!! @param[out] U       - value of the solution displacement
!! @param[out] GradU   - corresponding first derivatives
!! @param[out] Epsilon - value of the solution strain
!! @param[out] Sigma   - value of the solution stress
!
!-----------------------------------------------------------------------------------
!
  subroutine exact_solution(X,Ndom, U,GradU,Epsilon,Sigma)
!
  use common_prob_data, only : EPS
  use sheathed_isotropic_materials
  implicit none
!-----------------------------------------------------------------------------------
  real*8, dimension(3),   intent(in)  :: X
  integer,                intent(in)  :: Ndom
  real*8, dimension(3),   intent(out) :: U
  real*8, dimension(3,3), intent(out) :: GradU
  real*8, dimension(3,3), intent(out) :: Epsilon
  real*8, dimension(3,3), intent(out) :: Sigma
!-----------------------------------------------------------------------------------
  integer :: iprint,itest,i,j,k,l,icomp
  real*8 :: y,z,r,c1,c2,p,denom,theta,drdy,drdz,dthetady,dthetadz,   &
            ur,urdr,trr,tthetatheta,diffmax,dmax,tol
  real*8, dimension(3,3,3,3) :: A
  real*8, dimension(3,3)     :: tmpTensor
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
  iprint = 0
  itest  = 0
  tol = 1.d-14

! initialize variables (these stay equal to 0 if r < EPS)
  U = 0.d0; GradU = 0.d0; Epsilon = 0.d0; Sigma = 0.d0

! NOTE: this is always a solution of the homogeneous equation so DivSigma remains
!       equal to 0
!
!-----------------------------------------------------------------------------------
!      D E C L A R E    C O O R D I N A T E S    A N D    P A R A M E T E R S      |
!-----------------------------------------------------------------------------------
!
  y = X(2)
  z = X(3)
  r = dsqrt(y**2+z**2)
!
  if (r.lt.EPS) return
!
  theta = datan2(z,y)
!
  drdy = y/r ; dthetady =-z/r**2
  drdz = z/r ; dthetadz = y/r**2
!
!-----------------------------------------------------------------------------------
!      A S S E M B L E    S O L U T I O N                                          |
!-----------------------------------------------------------------------------------
!
!
  denom = R_inside**2*( R_outside**2*(MU_rubber-MU_steel)*(LAMBDA_steel+MU_steel)  &
                      + R_middle**2*MU_steel*(LAMBDA_steel+MU_steel+MU_rubber) )  &
        - R_middle**2*MU_rubber*( R_middle**2*MU_steel + R_outside**2*(LAMBDA_steel+MU_steel) )
  if (r.le.R_middle) then
    c1 = R_inside**2*R_middle**2  &
       *( P_inner*( R_middle**2*MU_steel + R_outside**2*(LAMBDA_steel+MU_steel) )  &
        - P_outer*R_outside**2*(LAMBDA_steel+2.d0*MU_steel) )
    p  = P_inner*R_inside**2*( R_outside**2*(MU_rubber-MU_steel)*(LAMBDA_steel+MU_steel)  &
                             + R_middle**2*MU_steel*(LAMBDA_steel+MU_rubber+MU_steel ) )  &
       - P_outer*R_middle**2*R_outside**2*MU_rubber*(LAMBDA_steel+2.d0*MU_steel)
    c1 = 0.5d0*c1/denom
    p  = p/denom

    ur          = c1/r
    urdr        =-c1/r**2
    trr         =-2.d0*MU_rubber*c1/r**2+p
    tthetatheta = 2.d0*MU_rubber*c1/r**2+p
    Sigma(1,1)  = p
  else
    c1 = R_middle**2*R_outside**2  &
       *( P_inner*R_inside**2*(LAMBDA_steel+MU_steel)  &
        - P_outer*( R_inside**2*(LAMBDA_steel+MU_steel+MU_rubber) - R_middle**2*MU_rubber ) )
    c2 = P_outer*R_outside**2*( R_inside**2*(MU_rubber-MU_steel) - R_middle**2*MU_rubber )  &
       + P_inner*R_inside**2*R_middle**2*MU_steel
    c1 = 0.5d0*c1/denom
    c2 = 0.5d0*c2/denom

    ur          = c1/r + c2*r
    urdr        = c2-c1/r**2
    trr         = 2.d0*MU_steel*(c2-c1/r**2) + 2.d0*LAMBDA_steel*c2
    tthetatheta = 2.d0*MU_steel*(c2+c1/r**2) + 2.d0*LAMBDA_steel*c2
    Sigma(1,1)  = 2.d0*LAMBDA_steel*c2
  endif

!
! displacements in Cartesian coordinates
  U(2) = ur*dcos(theta)
  U(3) = ur*dsin(theta)
!
! derivatives in Cartesian coordinates
  GradU(2,2) = urdr*dcos(theta)*drdy - ur*dsin(theta)*dthetady
  GradU(2,3) = urdr*dcos(theta)*drdz - ur*dsin(theta)*dthetadz
  GradU(3,2) = urdr*dsin(theta)*drdy + ur*dcos(theta)*dthetady
  GradU(3,3) = urdr*dsin(theta)*drdz + ur*dcos(theta)*dthetadz
!
! strain in Cartesian coordinates
  Epsilon(2,2) = GradU(2,2)
  Epsilon(3,3) = GradU(3,3)
  Epsilon(2,3) = 0.5d0*(GradU(2,3)+GradU(3,2)) ! may be improved
  Epsilon(3,2) = Epsilon(2,3)
!
! transform stresses from cylindrical to Cartesian coordinates
  Sigma(2,2) = trr*dcos(theta)**2 + tthetatheta*dsin(theta)**2
  Sigma(3,3) = trr*dsin(theta)**2 + tthetatheta*dcos(theta)**2
  Sigma(2,3) = (trr-tthetatheta)*dsin(theta)*dcos(theta)
  Sigma(3,2) = Sigma(2,3)

  !
  !   A:sigma = epsilon
  !
  if (itest.eq.1) then
    call getA(X,Ndom, A)

    tmpTensor = 0.d0
    do icomp=1,3
      do k=1,3; do l=1,3
        tmpTensor(icomp,1:3) = tmpTensor(icomp,1:3) + A(icomp,1:3,k,l)*Sigma(k,l)
      enddo; enddo
    enddo

    diffmax = 0.d0; dmax = 0.d0
    do i=1,3
      do j=i,3
        diffmax = max(diffmax,abs(tmpTensor(i,j)-Epsilon(i,j)))
        dmax = max(dmax,abs(Sigma(i,j)))
      enddo
    enddo
    if (diffmax/dmax.gt.tol) then
      write(*,7004) diffmax, dmax
 7004 format('exact_solution: diffmax, dmax FOR A:Sigma - Epsilon = ',2e12.5)
      call pause
    endif

    write(*,*) '+/- Sigma_n(2:3) = ', Sigma(2,2)*dcos(theta)+Sigma(2,3)*dsin(theta), &
                                      Sigma(2,3)*dcos(theta)+Sigma(3,3)*dsin(theta)
    write(*,*) 'y,z,r = ', y,z,r
    write(*,*) 'Sigma(2,2) = ', Sigma(2,2)
    write(*,*) 'Sigma(3,3) = ', Sigma(3,3)
    write(*,*) 'Sigma(2,3) = ', Sigma(2,3)
    write(*,*) 'tmpTensor = ', tmpTensor
    write(*,*) 'Epsilon = ', Epsilon
    call pause
  endif

!
  end subroutine exact_solution