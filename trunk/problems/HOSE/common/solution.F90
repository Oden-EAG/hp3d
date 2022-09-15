!-----------------------------------------------------------------------------------
!> Purpose : compute all relevant quantities of the exact solution for the
!            linear elasticty problem
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
  integer :: ndom
  real*8  :: r
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
!
! determine subdomain
  ndom=1
  r = dsqrt(X(2)**2+X(3)**2)
  if (r.lt.R_middle) ndom=2

! initialize variables
  U = 0.d0; GradU = 0.d0; Epsilon = 0.d0; Sigma = 0.d0; DivSigma = 0.d0

  call exact_solution(X,ndom, U,GradU,Epsilon,Sigma)

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