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
  use isotropic_elast_material_compo, only: find_material
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
  real*8, dimension(3)   :: u
  real*8, dimension(3,3) :: gradu
  real*8, dimension(3,3) :: epsilon
  real*8, dimension(3,3) :: sigma
  real*8, dimension(3)   :: divsigma
  integer :: imat
!------------------------------------------------------------------------------
! initialize exact solution
  ValH = 0.d0 ; DvalH = 0.d0 ; D2valH = 0.d0
  ValE = 0.d0 ; DvalE = 0.d0 ; D2valE = 0.d0
  ValV = 0.d0 ; DvalV = 0.d0 ; D2valV = 0.d0
  ValQ = 0.d0 ; DvalQ = 0.d0 ; D2valQ = 0.d0
!------------------------------------------------------------------------------
 call find_material(x,(/0.d0,0.d0,0.d0/),imat)
 call elast_solution(X,imat, u,gradu,epsilon,sigma,divsigma)
 !
 ValH = u
 DvalH = gradu
 ValV  = sigma  !  stress is symmetric so this works
 DvalV(1,1,1) = divsigma(1)
 DvalV(2,2,2) = divsigma(2)
 DvalV(3,3,3) = divsigma(3)
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
ValQ(10)  = gradu(2,3)-epsilon(2,3)
ValQ(11)  = -(gradu(1,3)-epsilon(1,3))
ValQ(12)  = gradu(1,2)-epsilon(1,2)
!
end subroutine exact