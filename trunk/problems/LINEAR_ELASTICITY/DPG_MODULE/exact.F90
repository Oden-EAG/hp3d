!------------------------------------------------------------------------------
!> Purpose : exact (manufactured, no miracle!) solution
!!
!! @param[in]  X       - a point in physical space
!! @param[in]  Icase   - node case (specifies what variables are supported)
!! @param[out] zvalH   - value of the H1 solution
!! @param[out] ZdvalH  - corresponding first derivatives
!! @param[out] Zd2valH - corresponding second derivatives
!! @param[out] ZvalE   - value of the H(curl) solution
!! @param[out] ZdvalE  - corresponding first derivatives
!! @param[out] Zd2valE - corresponding second derivatives
!! @param[out] ZvalV   - value of the H(div) solution
!! @param[out] ZdvalV  - corresponding first derivatives
!! @param[out] Zd2valV - corresponding second derivatives
!! @param[out] ZvalQ   - value of the H(div) solution
!! @param[out] ZdvalQ  - corresponding first derivatives
!! @param[out] Zd2valQ - corresponding second derivatives
!------------------------------------------------------------------------------
!
subroutine exact(X,Icase, ValH,DvalH,D2valH, &
                          ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, &
                          ValQ,DvalQ,D2valQ)
  use data_structure3D
  use common_prob_data, only : IEXACT_PROB,PI,NP1,NP2,NP3
  use isotropic_elast_material
  implicit none
!------------------------------------------------------------------------------
  real*8,dimension(3),             intent(in)  :: X
  integer,                         intent(in)  :: Icase
  real*8,dimension(  MAXEQNH    ), intent(out) ::   ValH
  real*8,dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
  real*8,dimension(  MAXEQNH,3,3), intent(out) :: D2valH
  real*8,dimension(3,MAXEQNE    ), intent(out) ::   ValE
  real*8,dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
  real*8,dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
  real*8,dimension(3,MAXEQNV    ), intent(out) ::   ValV
  real*8,dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
  real*8,dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
  real*8,dimension(  MAXEQNQ    ), intent(out) ::   ValQ
  real*8,dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
  real*8,dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!------------------------------------------------------------------------------
! workspace for singular solution
  real*8 :: SolV1(2),Soldiv1,SolV2(2),Soldiv2,SolQ1,SolQ2,SolQ3, &
            Solu11,Solu12,Solu21,Solu22
!------------------------------------------------------------------------------
  integer :: isol,icomp,k,l,m
  real*8  :: x1,x2,x3
  real*8, dimension(3,3,3,3) :: C
!------------------------------------------------------------------------------
! initialize exact solution
  ValH = 0.d0 ; DvalH = 0.d0 ; D2valH = 0.d0
  ValE = 0.d0 ; DvalE = 0.d0 ; D2valE = 0.d0
  ValV = 0.d0 ; DvalV = 0.d0 ; D2valV = 0.d0
  ValQ = 0.d0 ; DvalQ = 0.d0 ; D2valQ = 0.d0
!------------------------------------------------------------------------------
!
  x1 = X(1); x2 = X(2); x3 = X(3)
!
  select case(IEXACT_PROB)
! POLYNOMIAL SOLUTION (1/2)
  case(1)
!   value
    ValH(1:3) = x1*x2*x3
!
!   1st order derivatives
    DvalH(1:3,1) = x2*x3
    DvalH(1:3,2) = x1*x3
    DvalH(1:3,3) = x1*x2
!
!   2nd order derivatives
    D2valH(1:3,1,2) = x3; D2valH(1:3,1,3) = x2; D2valH(1:3,2,3) = x1
    D2valH(1:3,2,1) = D2valH(1:3,1,2)
    D2valH(1:3,3,1) = D2valH(1:3,1,3)
    D2valH(1:3,3,2) = D2valH(1:3,2,3)
!
! POLYNOMIAL SOLUTION (2/2)
  case(2)
!   value
    ValH(1) = x1**NP1
    ValH(2) = x2**NP2
    ValH(3) = x3**NP3
!
!   1st order derivatives
    DvalH(1,1) = NP1*x1**(NP1-1)
    DvalH(2,2) = NP2*x2**(NP2-1)
    DvalH(3,3) = NP3*x3**(NP3-1)
!
!   2nd order derivatives
    D2valH(1,1,1) = NP1*(NP1-1)*x1**(NP1-2)
    D2valH(2,2,2) = NP2*(NP2-1)*x2**(NP2-2)
    D2valH(3,3,3) = NP3*(NP3-1)*x3**(NP3-2)
!
! EXPONENTIAL SOLUTION (1/1)
  case(3)
!   value
    ValH(1) = 1.d0
    ValH(2) = 1.d0
    ValH(3) = dexp(x3)
!
!   1st order derivatives
    DvalH(3,3) = dexp(x3)
!
!   2nd order derivatives
    D2valH(3,3,3) = dexp(x3)
!
! SINUSOIDAL SOLUTION (1/1)
  case(4)
!
!   value
    ValH(1:3) = dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
!
!   1st order derivatives
    DvalH(1:3,1) = PI*dcos(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
    DvalH(1:3,2) = PI*dsin(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
    DvalH(1:3,3) = PI*dsin(PI*x1)*dsin(PI*x2)*dcos(PI*x3)
!
!   2nd order derivatives
    D2valH(1:3,1,1) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
    D2valH(1:3,1,2) = PI**2*dcos(PI*x1)*dcos(PI*x2)*dsin(PI*x3)
    D2valH(1:3,1,3) = PI**2*dcos(PI*x1)*dsin(PI*x2)*dcos(PI*x3)
    D2valH(1:3,2,1) = D2valH(1:3,1,2)
    D2valH(1:3,2,2) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
    D2valH(1:3,2,3) = PI**2*dsin(PI*x1)*dcos(PI*x2)*dcos(PI*x3)
    D2valH(1:3,3,1) = D2valH(1:3,1,3)
    D2valH(1:3,3,2) = D2valH(1:3,2,3)
    D2valH(1:3,3,3) =-PI**2*dsin(PI*x1)*dsin(PI*x2)*dsin(PI*x3)
!
! SINGULAR SOLUTION (USE NEXACT = 2 AND L_SHAPE GEOMETRY)
  case(5)
    call singular_sol_2D(x1,x2,                       &
                         solV1,soldiv1,solV2,soldiv2, &
                         solQ1,solQ2,solQ3,           &
                         solu11,solu12,solu21,solu22)
!   value
    ValH(1) = solQ1
    ValH(2) = solQ2
!
!   1st order derivatives
    DvalH(1,1) = solu11
    DvalH(1,2) = solu12
    DvalH(2,1) = solu21
    DvalH(2,2) = solu22
!
!   2nd order derivatives are not necessary (homogeneous solution)
  end select

  call getC(X, C)

  do icomp=1,3
    do k=1,3; do l=1,3
      ValV(1:3,icomp)      = C(icomp,1:3,k,l)*DvalH(k,l)
      do m=1,3
        DvalV(1:3,icomp,m) = C(icomp,1:3,k,l)*D2valH(k,l,m)
      enddo
    enddo; enddo
  enddo


end subroutine exact


!-----------------------------------------------------------------------
!
!   routine name       - singular_sol_2D
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 15
!   purpose            - routine evaluates a manufactured (exact)
!                        solution to the L-shape linear elasticity
!                        problem
!   arguments :
!     in:
!              Xtmp,Ytmp- coordinates of a point in 2D
!     out:
!              solV1    = (\sigma_xx, \sigma_xy)
!              Soldiv1  = \sigma_xx,x + \sigma_xy,y
!              SolV2    = (\sigma_yx, \sigma_yy)
!              Soldiv2  = \sigma_yx,x + \sigma_yy,y
!              SolQ1    = displacement u
!              SolQ2    = displacement v
!              SolQ3    = infinitezimal rotation
!                         \omega_xy = 1/2(u_y-v_x)
!              Solu11   = u,x
!              Solu12   = u,y
!              Solu21   = v,x
!              Solu22   = v,y
!
!-----------------------------------------------------------------------
!
  subroutine singular_sol_2D(Xtmp,Ytmp,                   &
                             SolV1,Soldiv1,SolV2,Soldiv2, &
                             SolQ1,SolQ2,SolQ3,           &
                             Solu11,Solu12,Solu21,Solu22)
!
  use common_prob_data, only : PI, EPS
  use isotropic_elast_material, only : LAMBDA,MU
  implicit none
!------------------------------------------------------------------------------
  real*8,               intent(in)  :: Xtmp
  real*8,               intent(in)  :: Ytmp
  real*8, dimension(2), intent(out) :: SolV1
  real*8,               intent(out) :: Soldiv1
  real*8, dimension(2), intent(out) :: SolV2
  real*8,               intent(out) :: Soldiv2
  real*8,               intent(out) :: SolQ1
  real*8,               intent(out) :: SolQ2
  real*8,               intent(out) :: SolQ3
  real*8,               intent(out) :: Solu11
  real*8,               intent(out) :: Solu12
  real*8,               intent(out) :: Solu21
  real*8,               intent(out) :: Solu22
!------------------------------------------------------------------------------
  integer :: iprint
  real*8 :: x,y,r,poisson,sigma,theta,drdx,drdy,dthetadx,dthetady,  &
            a,alpha,ak,F,Fd,Fdd,G,Gd,Gdd,b,ur,utheta,urdr,urdtheta, &
            uthetadr,uthetadtheta,trr,tthetatheta,trtheta,tmp1,tmp2
!------------------------------------------------------------------------------
  iprint = 0
!
! note that this is a solution of the homogeneous equation
  Soldiv1 = 0.d0; Soldiv2 = 0.d0
!
! rotate coordinates 45 degrees counter-clockwise
  x = (Xtmp+Ytmp)/dsqrt(2.d0)
  y = (-Xtmp+Ytmp)/dsqrt(2.d0)
!
  r = dsqrt(x**2+y**2)
  if (r.lt.EPS) return
!
  if (iprint.eq.1) then
    write(*,*) 'MU,LAMBDA = ',MU,LAMBDA
  endif
  poisson = LAMBDA/(2.d0*(LAMBDA+MU))
  sigma   = poisson/(1.d0+poisson)
  if (y.eq.0.d0) then
    if (x.lt.0.d0) then
      theta = PI
    elseif (x.gt.0.d0) then
      theta = 0.d0
    endif
  elseif (y.gt.0.d0) then
    theta = dacos(x/r)
  else
    theta =-dacos(x/r)
  endif
!
  drdx = x/r; dthetadx =-y/r**2
  drdy = y/r; dthetady = x/r**2
!
  ! a = 0.60404d0
  a = 0.6037781005214915d0
!
  alpha = 3.d0*PI/4.d0
  ak = (4.d0*(1-sigma)-(a+1.d0))*dsin((a-1.d0)*alpha) &
      /((a+1.d0)*dsin((a+1.d0)*alpha))
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
  ur     = b*r**a*(-(a+1.d0)*F+(1-sigma)*Gd)
  utheta = b*r**a*(-Fd+(1.d0-sigma)*(a-1.d0)*G)
  urdr         = a* b*r**(a-1.d0)*(-(a+1.d0)*F+(1-sigma)*Gd)
  urdtheta     = b*r**a*(-(a+1.d0)*Fd+(1-sigma)*Gdd)
  uthetadr     = a*b*r**(a-1.d0)*(-Fd+(1.d0-sigma)*(a-1.d0)*G)
  uthetadtheta = b*r**a*(-Fdd+(1.d0-sigma)*(a-1.d0)*Gd)
  if (iprint.eq.1) then
    write(*,7001) ur,utheta,theta
7001 format('dsingular: ur,utheta,theta = ',3e12.5)
  endif
!
! trr,tthetatheat,trtheta are the cylindrical stress components
  trr         = r**(a-1.d0)*(Fdd+(a+1.d0)*F)
  tthetatheta = a*(a+1.d0)*r**(a-1.d0)*F
  trtheta     =-a*r**(a-1.d0)*Fd
!
! displacements in the rotated system of coordinates
  SolQ1 = ur*dcos(theta) - utheta*dsin(theta)
  SolQ2 = ur*dsin(theta) + utheta*dcos(theta)
!
! transform displacements to the rotated system of coordinates
  SolQ1 = (SolQ1-SolQ2)/dsqrt(2.d0)
  SolQ2 = (SolQ1+SolQ2)/dsqrt(2.d0)
!
! and their derivatives in the rotated system of coordinates
  Solu11 = (urdr*dcos(theta)-uthetadr*dsin(theta))*drdx &
         + (urdtheta*dcos(theta)-ur*dsin(theta)-uthetadtheta*dsin(theta) &
          -utheta*dcos(theta))*dthetadx
  Solu12 = (urdr*dcos(theta)-uthetadr*dsin(theta))*drdy &
         + (urdtheta*dcos(theta)-ur*dsin(theta)-uthetadtheta*dsin(theta) &
          -utheta*dcos(theta))*dthetady
  Solu21 = (urdr*dsin(theta)+uthetadr*dcos(theta))*drdx &
         + (urdtheta*dsin(theta)+ur*dcos(theta)+uthetadtheta*dcos(theta) &
          -utheta*dsin(theta))*dthetadx
  Solu22 = (urdr*dsin(theta)+uthetadr*dcos(theta))*drdy &
         + (urdtheta*dsin(theta)+ur*dcos(theta)+uthetadtheta*dcos(theta) &
          -utheta*dsin(theta))*dthetady
!
! infinitesimal rotation component
  tmp1  = (urdr*dcos(theta)-uthetadr*dsin(theta))*drdy &
        + (urdtheta*dcos(theta)-uthetadtheta*dsin(theta))*dthetady &
        + (-ur*dsin(theta)-utheta*dcos(theta))*dthetady
  tmp2  = (urdr*dsin(theta)+uthetadr*dcos(theta))*drdx &
        + (urdtheta*dsin(theta)+uthetadtheta*dcos(theta))*dthetadx &
        + (ur*dcos(theta)-utheta*dsin(theta))*dthetadx
  SolQ3 = (tmp1-tmp2)/2.d0
!
! transform displacement derivatives to the original system
  Solu11 = 0.5d0*(Solu11-Solu21-Solu12+Solu22)
  Solu12 = 0.5d0*(Solu11-Solu21+Solu12-Solu22)
  Solu21 = 0.5d0*(Solu11+Solu21-Solu12-Solu22)
  Solu22 = 0.5d0*(Solu11+Solu21+Solu12+Solu22)
!
! transform stresses from the cylindrical to the rotated Cartesian
! system of coordinates
  SolV1(1) = trr*(dcos(theta)**2) &
           - 2.d0*dsin(theta)*dcos(theta)*trtheta &
           + tthetatheta*(dsin(theta)**2)
  SolV2(2) = trr*(dsin(theta)**2) &
           + 2.d0*dsin(theta)*dcos(theta)*trtheta &
           + tthetatheta*(dcos(theta)**2)
  SolV1(2) = trr*dsin(theta)*dcos(theta) &
           + (dcos(theta)**2-dsin(theta)**2)*trtheta &
           - tthetatheta*dsin(theta)*dcos(theta)
  SolV2(1) = SolV1(2)
!
! transform stresses to the original system of coordinates
  SolV1(1) = 0.5d0*(SolV1(1)-SolV2(1)-SolV1(2)+SolV2(2))
  SolV1(2) = 0.5d0*(SolV1(1)-SolV2(1)+SolV1(2)-SolV2(2))
  SolV2(1) = 0.5d0*(SolV1(1)+SolV2(1)-SolV1(2)-SolV2(2))
  SolV2(2) = 0.5d0*(SolV1(1)+SolV2(1)+SolV1(2)+SolV2(2))
!
!
  end subroutine singular_sol_2D