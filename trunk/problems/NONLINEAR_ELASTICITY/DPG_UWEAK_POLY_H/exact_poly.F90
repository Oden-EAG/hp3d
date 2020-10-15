!------------------------------------------------------------------------------
!> Purpose : exact (manufactured, no miracle!) solution
!!
!! @param[in]  X       - a point in physical space
!! @param[in]  Fn      - a vector (usually a unit normal for a face)
!! @param[in]  Idom    - Domain index
!! @param[out] zvalH   - value of the H1 solution
!! @param[out] ZdvalH  - corresponding first derivatives
!! @param[out] Zd2valH - corresponding second derivatives
!! @param[out] ZvalE   - value of the H(curl) solution
!! @param[out] ZdvalE  - corresponding first derivatives
!! @param[out] Zd2valE - corresponding second derivatives
!! @param[out] ZvalV   - value of the H(div) solution
!! @param[out] ZdvalV  - corresponding first derivatives
!! @param[out] Zd2valV - corresponding second derivatives
!! @param[out] ZvalU   - trace of H1 solutions
!! @param[out] ZvalF   - normal trace of H(div) solutions
!! @param[out] ZvalQ   - value of the H(div) solution
!! @param[out] ZdvalQ  - corresponding first derivatives
!! @param[out] Zd2valQ - corresponding second derivatives
!------------------------------------------------------------------------------
!
subroutine exact(X,Fn,Idom, ValH,DvalH,D2valH, &
                          ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, &
                          ValU,ValF,ValQ,DvalQ,D2valQ)
  use data_structure3D_poly
  implicit none
!------------------------------------------------------------------------------
  real*8, dimension(3),             intent(in)  :: X, Fn
  integer,                          intent(in)  :: Idom
  real*8, dimension(  MAXEQNH    ), intent(out) ::   ValH
  real*8, dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
  real*8, dimension(  MAXEQNH,3,3), intent(out) :: D2valH
  real*8, dimension(3,MAXEQNE    ), intent(out) ::   ValE
  real*8, dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
  real*8, dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
  real*8, dimension(3,MAXEQNV    ), intent(out) ::   ValV
  real*8, dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
  real*8, dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
  real*8, dimension(  MAXEQNU    ), intent(out) ::   ValU
  real*8, dimension(  MAXEQNF    ), intent(out) ::   ValF
  real*8, dimension(  MAXEQNQ    ), intent(out) ::   ValQ
  real*8, dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
  real*8, dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!------------------------------------------------------------------------------
  real*8, dimension(3)   :: u
  real*8, dimension(3,3) :: gradu
  real*8, dimension(3,3) :: epsilon
  real*8, dimension(3,3) :: sigma
  real*8, dimension(3)   :: divsigma
  integer                :: imat
!------------------------------------------------------------------------------
! initialize exact solution
  ValH = 0.d0 ; DvalH = 0.d0 ; D2valH = 0.d0
  ValE = 0.d0 ; DvalE = 0.d0 ; D2valE = 0.d0
  ValV = 0.d0 ; DvalV = 0.d0 ; D2valV = 0.d0 ; ValF = 0.d0
  ValQ = 0.d0 ; DvalQ = 0.d0 ; D2valQ = 0.d0 ; ValU = 0.d0
!------------------------------------------------------------------------------

  ! call find_material(X,Fn,imat)
  select case(idom)
  case(1)
    imat = 1
  case default
    imat = 2
  end select


  call elast_solution(imat , X, u,gradu,epsilon,sigma,divsigma)

  ValQ(1:3)  = u
!          ( sigma1   sigma4  sigma5 )
!  sigma = ( sigma4   sigma2  sigma6 )
!          ( sigma5   sigma6  sigma3 )
  ValQ(4)  = sigma(1,1)
  ValQ(5)  = sigma(2,2)
  ValQ(6)  = sigma(3,3)
  ValQ(7)  = sigma(1,2)
  ValQ(8)  = sigma(1,3)
  ValQ(9)  = sigma(2,3)
!          (   0     omega3  -omega2 )    1
!  omega = (-omega3    0      omega1 ) = --- ( grad(u)-grad(u)^T)
!          ( omega2  -omega1    0    )    2
  ValQ(10)  =   gradu(2,3)-epsilon(2,3)
  ValQ(11)  =   gradu(3,1)-epsilon(3,1)
  ValQ(12)  =   gradu(1,2)-epsilon(1,2)

  ValU(1:3) = u
  
  ValF(1  ) = sigma(1,1)*Fn(1)+sigma(1,2)*Fn(2)+sigma(1,3)*Fn(3)
  ValF(2  ) = sigma(2,1)*Fn(1)+sigma(2,2)*Fn(2)+sigma(2,3)*Fn(3)
  ValF(3  ) = sigma(3,1)*Fn(1)+sigma(3,2)*Fn(2)+sigma(3,3)*Fn(3)

end subroutine exact