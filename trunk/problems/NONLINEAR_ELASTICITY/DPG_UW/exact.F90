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
subroutine exact(X,Mdle,Icase,ValH,DvalH,D2valH, &
                              ValE,DvalE,D2valE, &
                              ValV,DvalV,D2valV, &
                              ValQ,DvalQ,D2valQ)
  use data_structure3D  
  use hyperelasticity, only: LOAD_FACTOR,DEL
  implicit none
!------------------------------------------------------------------------------
  real*8,dimension(3),             intent(in)  :: X
  integer,                         intent(in)  :: Mdle
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
  real*8, dimension(3,3) :: stress
  real*8, dimension(3)   :: divstress
!------------------------------------------------------------------------------
! initialize exact solution
  ValH = 0.d0 ; DvalH = 0.d0 ; D2valH = 0.d0
  ValE = 0.d0 ; DvalE = 0.d0 ; D2valE = 0.d0
  ValV = 0.d0 ; DvalV = 0.d0 ; D2valV = 0.d0
  ValQ = 0.d0 ; DvalQ = 0.d0 ; D2valQ = 0.d0
!------------------------------------------------------------------------------

 call hyperelast_solution(Mdle,LOAD_FACTOR,X, u,gradu,stress,divstress)

 ValH(4:6)  = u(1:3)
 DvalH(4:6,:) = gradu(1:3,:)

! for the correct storage of stress as Hdiv rows, we need the transpose
 ValV(:,4:6)  = transpose(stress(:,:))

 ValQ(NRQVAR+1 :NRQVAR+3 ) = u(1:3)
! reshape stress tensor into a vector of lenght 9
 ValQ(NRQVAR+4:NRQVAR+12) = reshape(stress,(/9/))
 ! save deformation gradient in gradu; then reshape tensor into a vector of lenght 9
 ! gradu = DEL + gradu
 ValQ(NRQVAR+13:NRQVAR+21) = reshape(gradu ,(/9/)) 

end subroutine exact