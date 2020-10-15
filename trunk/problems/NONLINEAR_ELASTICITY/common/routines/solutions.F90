!-----------------------------------------------------------------------------------
!> Purpose : compute all relevant quantities of the exact solutions for the
!            linear elasticty problems
!
!! @param[in]  Mdle     - element middle node number
!! @param[in]  Factor   - load factor to multiply prescribed displacements
!! @param[in]  X        - a point in physical space
!! @param[out] U        - value of the solution displacement
!! @param[out] GradU    - corresponding first derivatives
!! @param[out] PK1      - value of the first Piola-Kirchhoff stress
!! @param[out] DivPK1   - divergence of the first Piola-Kirchhoff stress
!-----------------------------------------------------------------------------------
!
subroutine hyperelast_solution(Mdle, Factor, X, U,GradU,PK1,DivPK1)
  use data_structure3D
  use common_prob_data, only : IEXACT_PROB,PI,NP1,NP2,NP3
  use hyperelasticity

  implicit none
!-----------------------------------------------------------------------------------
  integer,                intent(in)  :: Mdle
  real*8,                 intent(in)  :: Factor
  real*8, dimension(3),   intent(in)  :: X
  real*8, dimension(3),   intent(out) :: U
  real*8, dimension(3,3), intent(out) :: GradU
  real*8, dimension(3,3), intent(out) :: PK1
  real*8, dimension(3),   intent(out) :: DivPK1
!-----------------------------------------------------------------------------------
! local variables
  character*6 :: prob
  integer :: iprint,itest,icomp,i,j,k,l,m,     imat,np
  real*8  :: x1,x2,x3,tol,diffmax,dmax ,                               &
             ! pres , e , nu , grad_pres , pres_lambda ,                 &
             pres_tilde , rad, the, phi, mu_m, mu_s, e_s, e_m, nu_s,   &
             nu_m , lambda_s, alpha , W , p(MAX_NR_P)
  real*8, dimension(3)       :: sph
  real*8, dimension(3,3)     :: TmpTensor,    transf, Ftensor, Epsilon(3,3)
  real*8, dimension(3,3,3)   :: Grad2U,GradPK1
  real*8, dimension(3,3,3,3) :: C,A
!------------------------------------------------------------------------------
! workspace for singular solution
  real*8 :: SolV1(2),Soldiv1,SolV2(2),Soldiv2,SolQ1,SolQ2,SolQ3, &
            Solu11,Solu12,Solu21,Solu22
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
  iprint = 0
  itest  = 0
  tol = 1.0d-10

! get stiffness and compliance tensors
  ! call getC(X, C)
  ! call getA(X, A)

! initialize variables
  U = 0.d0; GradU = 0.d0; PK1 = 0.d0; DivPK1 = 0.d0
  Grad2U = 0.d0; GradPK1 = 0.d0
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
!
  select case(IEXACT_PROB)

!   case(999)
!     ! s: sphere
!     nu_s = PN(1)
!     e_s  = EE(1)
!     mu_s = 0.5d0*e_s/(1.d0+nu_s)
!     lambda_s = (e_s*nu_s)/(1.d0-nu_s-2.d0*nu_s**2)
!     ! m: matrix
!     nu_m = PN(2)
!     e_m  = EE(2)
!     mu_m = e_m/3.d0
!     ! hydrostatic pressure
!     pres = 10.d0
!     ! constant for manufactured solution
!     alpha = pres * (1.d0-2.d0*nu_s)/(32.d0*e_s+(1.d0-2.d0*nu_s)*128*e_m/3.d0)

!     pres_tilde = pres - 128.d0 * mu_m * alpha

!     call coord_cart2sphere(X,sph)
!     rad= sph(1) ! radius
!     the= sph(2) ! inclination
!     phi= sph(3) ! azimuthal
    
!     transf(1:3,1) = (/ sin(the)*cos(phi), sin(the)*sin(phi),  cos(the) /)
!     transf(1:3,2) = (/ cos(the)*cos(phi), cos(the)*sin(phi), -sin(the) /)
!     transf(1:3,3) = (/         -sin(phi),          cos(phi),      0.d0 /)

!     ! select case(imat)
!     ! case (1) ! sphere
!     if (rad.le.0.25d0) then
!       u(1) = - pres_tilde * rad *nu_s / ((1+nu_s)*lambda_s)
!       u    = matmul(transf,u)
!       Epsilon(1,1) = - pres_tilde *nu_s / ((1+nu_s)*lambda_s)
!       Epsilon(2,2) = Epsilon(1,1)
!       Epsilon(3,3) = Epsilon(1,1)         
!       PK1   = -pres_tilde * del
!       ! Now transform to Cartesian
!       Epsilon = matmul(transf,matmul(Epsilon,transpose(transf)))
!       GradU   = Epsilon   
!       PK1   = matmul(transf,matmul(PK1,transpose(transf)))
!     else
!     ! case (2) ! matrix (incompressible)
!       ! if (rad.lt.0.2d0) then
!       !   write(*,*) 'INCORRECT MATERIAL'
!       ! endif
!       ! if (imat.ne.2) write(*,*) 'INCORRECT MATERIAL. imat,rad =',imat,x
!       u(1) = -alpha/(2.d0*rad**2)
!       u    = matmul(transf,u)
!       Epsilon(1,1) = alpha/(rad**3)
!       Epsilon(2,2) = -alpha/(2.d0*rad**3)
!       Epsilon(3,3) = Epsilon(2,2)
!       PK1   = - pres * del + 2 * mu_m * Epsilon
!       ! Now transform to Cartesian
!       Epsilon = matmul(transf,matmul(Epsilon,transpose(transf)))
!       GradU   = Epsilon      
!       PK1   = matmul(transf,matmul(PK1,transpose(transf)))
!     ! end select
!     endif


!   case(123)

!     nu = PN(iMat)
!     e  = EE(iMat)
!     mu     = 0.5d0*e/(1.d0+nu)

!     pres = 10.d0
!     grad_pres = 0.d0

!     if (abs(nu-0.5d0).lt.tol) then
!       pres_lambda = 0.d0
!     else
!       lambda = (e*nu)/(1.d0-nu-2.d0*nu**2)    
!       pres_lambda = pres / lambda
!     endif

!     U(1) = 3.d0*x2
!     U(2) = 3.d0*x3
!     U(3) = 3.d0*x1

!     U(:) = U(:) - pres_lambda * X(:) / 3.d0

! !   1st order derivatives
!     GradU(1,2) = 3.d0
!     GradU(2,3) = 3.d0
!     GradU(3,1) = 3.d0

!     GradU = GradU - pres_lambda * del / 3.d0


!     Epsilon(1,2)= 1.5d0
!     Epsilon(1,3)= 1.5d0
!     Epsilon(2,1)= 1.5d0  
!     Epsilon(2,3)= 1.5d0
!     Epsilon(3,1)= 1.5d0
!     Epsilon(3,2)= 1.5d0

!     Epsilon = Epsilon - pres_lambda * del / 3.d0

!     PK1 = - pres * del + 2.d0 * mu * Epsilon

!     DivPK1 = - grad_pres !+ 2.d0 * mu * Grad2U

! constant displacement
    case(99)
    U(1) = 0.01d0
    U(2) = 0.02d0
    U(3) = 0.03d0
! POLYNOMIAL SOLUTION (0/4)
   case(0)
!   displacement
    U(1) = x2
    U(2) = x3
    U(3) = x1
!
!   1st order derivatives
    GradU(1,2) = 1.d0
    GradU(2,3) = 1.d0
    GradU(3,1) = 1.d0
!
! POLYNOMIAL SOLUTION (1/4)
  case(1)
!   displacement
    U(1:3) = x1*x2*x3
!
!   1st order derivatives
    GradU(1:3,1) = x2*x3
    GradU(1:3,2) = x1*x3
    GradU(1:3,3) = x1*x2
!
!   2nd order derivatives
    Grad2U(1:3,1,2) = x3; Grad2U(1:3,1,3) = x2; Grad2U(1:3,2,3) = x1
    Grad2U(1:3,2,1) = Grad2U(1:3,1,2)
    Grad2U(1:3,3,1) = Grad2U(1:3,1,3)
    Grad2U(1:3,3,2) = Grad2U(1:3,2,3)
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
! model linear solutions
! equibiaxial strain directions 1 2
  case(11)
!   displacement
    U(1) = x1
    U(2) = x2
!
!   1st order derivatives
    GradU(1,1) = 1.d0
! 
! simple shear 12
  case(12)
!   displacement
    U(1) = x2
!
!   1st order derivatives
    GradU(1,2) = 1.d0
! 
! uniaxial strain direction 3
  case(33)
!   displacement
    U(3) = x3
!
!   1st order derivatives
    GradU(3,3) = 1.d0
! 
! compression
  case(321)
!   displacement
    U(1) = -0.1d0*x1
    U(2) = -0.1d0*x2
    U(3) = -0.1d0*x3
!
!   1st order derivatives
    GradU(1,1) = -0.1d0
    GradU(2,2) = -0.1d0
    GradU(3,3) = -0.1d0
! 
! EXPONENTIAL SOLUTION (1/2)
  case(3)
!   displacement
    U(1) = 0.d0
    U(2) = 0.d0
    U(3) = dexp(x3)-1.d0
!
!   1st order derivatives
    GradU(3,3) = dexp(x3)
!
!   2nd order derivatives
    Grad2U(3,3,3) = dexp(x3)
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
! SINGULAR SOLUTION (1/2)
!  NOTE:: use NEXACT=2 and L-shape geometry
  case(5)
!   PLANE STRAIN PROBLEM
    prob = 'strain'

    if (MATERIALS(Imat)%CONSTIT.ne.LINEAR) then
      write(*,*) 'hyperelast_solution: singular solution supports linear materials only!'
      stop
    endif

    np = MATERIALS(Imat)%NR_PARAMS
    p = 0.d0
    p(1:np) = MATERIALS(Imat)%PARAMS

    call singular_solution(X(1:2),p(1),p(2),prob, U,GradU,Epsilon,PK1)

  end select
!
!-----------------------------------------------------------------------------------
!      C A L C U L A T E    D E R I V E D    V A R I A B L E S                     |
!-----------------------------------------------------------------------------------
!
  select case(IEXACT_PROB)
!
! Here the strain and stress have to be calculated

  case(0,1,2,3,4,30,11,12,33,321)


    U = Factor * U

    GradU = Factor * GradU

    Grad2U= Factor * Grad2U

    Ftensor = DEL + GradU
    
    call find_material(Mdle,imat)

    call eval_strain_energy_W_F(imat,X,Ftensor,W,PK1,C)
    do icomp=1,3
      do k=1,3
        !  calculate strain
        Epsilon(icomp,k) = 0.5d0*(GradU(icomp,k)+GradU(k,icomp))
        do l=1,3
          !  calculate Cauchy stress
          ! PK1(icomp,1:3) = PK1(icomp,1:3) + C(icomp,1:3,k,l)*GradU(k,l)
          do m=1,3
            GradPK1(icomp,1:3,m) = GradPK1(icomp,1:3,m) + C(icomp,1:3,k,l)*Grad2U(k,l,m)
          enddo
        enddo
      enddo
      !  calculate divergence of the stress
      do m=1,3
        DivPK1(icomp) = DivPK1(icomp) + GradPK1(icomp,m,m)
      enddo
    enddo
!
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
        diffmax = max(diffmax,abs(PK1(i,j)-PK1(j,i)))
        dmax = max(dmax,abs(PK1(i,j)))
      enddo
    enddo
    if (diffmax/dmax.gt.tol) then
      write(*,7001) diffmax, dmax
 7001 format('elast_solution: diffmax, dmax FOR PK1 = ',2e12.5)
      call pause
    endif
!
!   C:grad(u) = PK1
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
          diffmax = max(diffmax,abs(tmpTensor(i,j)-PK1(i,j)))
          dmax = max(dmax,abs(PK1(i,j)))
        enddo
      enddo
      if (diffmax/dmax.gt.tol) then
        write(*,7002) diffmax, dmax
   7002 format('elast_solution: diffmax, dmax FOR C:grad(u) - PK1 = ',2e12.5)
        call pause
      endif
    endif
!
!   C:epsilon = PK1
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
        diffmax = max(diffmax,abs(tmpTensor(i,j)-PK1(i,j)))
        dmax = max(dmax,abs(PK1(i,j)))
      enddo
    enddo
    if (diffmax/dmax.gt.tol) then
      write(*,7003) diffmax, dmax
 7003 format('elast_solution: diffmax, dmax FOR C:epsilon - PK1 = ',2e12.5)
      call pause
    endif
!
!   A:PK1 = epsilon
!
    tmpTensor = 0.d0
    do icomp=1,3
      do k=1,3; do l=1,3
        tmpTensor(icomp,1:3) = tmpTensor(icomp,1:3) + A(icomp,1:3,k,l)*PK1(k,l)
      enddo; enddo
    enddo

    diffmax = 0.d0; dmax = 0.d0
    do i=1,3
      do j=i,3
        diffmax = max(diffmax,abs(tmpTensor(i,j)-Epsilon(i,j)))
        dmax = max(dmax,abs(PK1(i,j)))
      enddo
    enddo
    if (diffmax/dmax.gt.tol) then
      write(*,7004) diffmax, dmax
 7004 format('elast_solution: diffmax, dmax FOR A:PK1 - Epsilon = ',2e12.5)
      call pause
    endif
  endif

 !  call singular_sol_df(x1,x2,                       &
 !                      solV1,soldiv1,solV2,soldiv2, &
 !                      solQ1,solQ2,solQ3,           &
 !                      solu11,solu12,solu21,solu22)

 !  write(*,*) 'solV1 = ', solV1
 !  write(*,*) 'solV2 = ', solV2
 !  write(*,*) 'solQ1 = ', solQ1
 !  write(*,*) 'solQ2 = ', solQ2
 !  write(*,*) 'solQ3 = ', solQ3
 !  write(*,*) 'solu11 = ', solu11
 !  write(*,*) 'solu12 = ', solu12
 !  write(*,*) 'solu21 = ', solu21
 !  write(*,*) 'solu22 = ', solu22

 !    write(*,*) ''
 !    write(*,*) '             PK1           '
 !    write(*,1000) PK1
 ! 1000 format(3(/,3e14.7))

 !    write(*,*) ''
 !    write(*,*) '             grad(u)          '
 !    write(*,1000) GradU

 !    write(*,*) GradU(1,2)-solu12

end subroutine hyperelast_solution

!-----------------------------------------------------------------------------------
!> Purpose : evaluate a manufactured solution for the L-shape linear elasticity
!            problem
!!
!! @param[in]  X       - a point in physical space
!! @param[out] U       - value of the solution displacement
!! @param[out] GradU   - corresponding first derivatives
!! @param[out] Epsilon - value of the solution strain
!! @param[out] PK1   - value of the solution stress
!
!-----------------------------------------------------------------------------------
!
  subroutine singular_solution(Xtmp,Lambda,Mu,prob, U,GradU,Epsilon,PK1)
!
  use common_prob_data, only : PI, EPS
  ! use isotropic_elast_material, only : LAMBDA,MU
  implicit none
!-----------------------------------------------------------------------------------
  real*8, dimension(2),   intent(in)  :: Xtmp
  real*8,                 intent(in)  :: Lambda,Mu
  character*6,            intent(in)  :: prob
  real*8, dimension(3),   intent(out) :: U
  real*8, dimension(3,3), intent(out) :: GradU
  real*8, dimension(3,3), intent(out) :: Epsilon
  real*8, dimension(3,3), intent(out) :: PK1
!-----------------------------------------------------------------------------------
  integer :: iprint
  real*8 :: x,y,r,poisson,sig,theta,drdx,drdy,dthetadx,dthetady,   &
            a,alpha,ak,F,Fd,Fdd,G,Gd,Gdd,b,ur,utheta,urdr,urdtheta,  &
            uthetadr,uthetadtheta,trr,tthetatheta,trtheta,tmp1,tmp2, &
            u11tmp,u12tmp,u21tmp,u22tmp,v11tmp,v12tmp,v22tmp
  real*8, parameter :: inf = 1.0d32
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
  iprint = 0

! initialize variables (these stay equal to 0 if r < EPS)
  U = 0.d0; GradU = 0.d0; Epsilon = 0.d0; PK1 = 0.d0

! NOTE: this is always a solution of the homogeneous equation so DivPK1 remains
!       equal to 0
!
!-----------------------------------------------------------------------------------
!      D E C L A R E    C O O R D I N A T E S    A N D    P A R A M E T E R S      |
!-----------------------------------------------------------------------------------
!
! rotate coordinates 45 degrees counter-clockwise
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
  F   = ak*dsin((a+1.d0)*theta) + dsin((a-1.d0)*theta)
  Fd  = (a+1.d0)*ak*dcos((a+1.d0)*theta) &
      + (a-1.d0)*dcos((a-1.d0)*theta)
  Fdd =-(a+1.d0)**2*ak*dsin((a+1.d0)*theta) &
      - (a-1.d0)**2*dsin((a-1.d0)*theta)
  G   =-4.d0/(a-1.d0)*dcos((a-1.d0)*theta)
  Gd  = 4.d0*dsin((a-1.d0)*theta)
  Gdd = 4.d0*(a-1.d0)*dcos((a-1.d0)*theta)
!
  b      = 1.d0/(2.d0*MU)
  ur     = b*r**a*(-(a+1.d0)*F+(1-sig)*Gd)
  utheta = b*r**a*(-Fd+(1.d0-sig)*(a-1.d0)*G)
  urdr         = a*b*r**(a-1.d0)*(-(a+1.d0)*F+(1-sig)*Gd)
  urdtheta     = b*r**a*(-(a+1.d0)*Fd+(1-sig)*Gdd)
  uthetadr     = a*b*r**(a-1.d0)*(-Fd+(1.d0-sig)*(a-1.d0)*G)
  uthetadtheta = b*r**a*(-Fdd+(1.d0-sig)*(a-1.d0)*Gd)
  if (iprint.eq.1) then
    write(*,7001) ur,utheta,theta
7001 format('singular_solution : ur,utheta,theta = ',3e12.5)
  endif
!
! trr,tthetatheat,trtheta are the cylindrical stress components
  trr         = r**(a-1.d0)*(Fdd+(a+1.d0)*F)
  tthetatheta = a*(a+1.d0)*r**(a-1.d0)*F
  trtheta     =-a*r**(a-1.d0)*Fd
!
! displacements in the rotated system of coordinates
  tmp1 = ur*dcos(theta) - utheta*dsin(theta)
  tmp2 = ur*dsin(theta) + utheta*dcos(theta)
!
! transform displacements to the rotated system of coordinates
  U(1) = (tmp1-tmp2)/dsqrt(2.d0)
  U(2) = (tmp1+tmp2)/dsqrt(2.d0)
! !
! ! infinitesimal rotation component
!   tmp1  = (urdr*dcos(theta)-uthetadr*dsin(theta))*drdy &
!         + (urdtheta*dcos(theta)-uthetadtheta*dsin(theta))*dthetady &
!         + (-ur*dsin(theta)-utheta*dcos(theta))*dthetady
!   tmp2  = (urdr*dsin(theta)+uthetadr*dcos(theta))*drdx &
!         + (urdtheta*dsin(theta)+uthetadtheta*dcos(theta))*dthetadx &
!         + (ur*dcos(theta)-utheta*dsin(theta))*dthetadx
!   SolQ3 = (tmp1-tmp2)/2.d0
!
! and their derivatives in the rotated system of coordinates
  u11tmp = (urdr*dcos(theta)-uthetadr*dsin(theta))*drdx  &
         + (urdtheta*dcos(theta)-ur*dsin(theta)-uthetadtheta*dsin(theta)  &
          -utheta*dcos(theta))*dthetadx
  u12tmp = (urdr*dcos(theta)-uthetadr*dsin(theta))*drdy  &
         + (urdtheta*dcos(theta)-ur*dsin(theta)-uthetadtheta*dsin(theta)  &
          -utheta*dcos(theta))*dthetady
  u21tmp = (urdr*dsin(theta)+uthetadr*dcos(theta))*drdx  &
         + (urdtheta*dsin(theta)+ur*dcos(theta)+uthetadtheta*dcos(theta)  &
          -utheta*dsin(theta))*dthetadx
  u22tmp = (urdr*dsin(theta)+uthetadr*dcos(theta))*drdy  &
         + (urdtheta*dsin(theta)+ur*dcos(theta)+uthetadtheta*dcos(theta)  &
          -utheta*dsin(theta))*dthetady
!
! transform displacement derivatives to the original system
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
! transform stresses from the cylindrical to the rotated Cartesian
! system of coordinates
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
! transform stresses to the original system of coordinates
  PK1(1,1) = 0.5d0*(v11tmp-2*v12tmp+v22tmp)
  PK1(1,2) = 0.5d0*(v11tmp-v22tmp)
  PK1(2,2) = 0.5d0*(v11tmp+2*v12tmp+v22tmp)
  PK1(2,1) = PK1(1,2)
!
  select case(prob)
  case('stress')
    U(3) = inf ! need formula for this component
    GradU(3,3) =-LAMBDA/(2.0d0*MU*(2.0d0*MU+3.0d0*LAMBDA))  &
                *(PK1(1,1)+PK1(2,2))
    Epsilon(3,3) = GradU(3,3)
  case('strain')
    PK1(3,3) = LAMBDA*(Epsilon(1,1)+Epsilon(2,2))
  case default
    write(*,*) 'singular_solution : NOT A VALID SOLUTION TYPE'
    stop 1
  end select
!
  end subroutine singular_solution