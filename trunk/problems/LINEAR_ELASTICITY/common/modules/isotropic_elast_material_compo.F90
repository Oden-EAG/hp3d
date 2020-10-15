!--------------------------------------------------------------------------------
! MODULE FOR COMPUTATION OF ELASTICITY AND COMPLIANCE TENSORS (C and S=C^{-1})
! FOR A COMPOSITE OF TWO MATERIALS
!
! Calculation of C and S(here labeled A) is made w.r.t Poisson coefficcient \nu 
! and modulus of elasticity E, with the intention of covering incompressibility.
!
!--------------------------------------------------------------------------------
!
module isotropic_elast_material_compo

  real*8, parameter :: RHO   = 0.d0
  real*8, parameter :: OMEG  = 0.d0

  real*8, parameter :: PN(2) = (/ 0.195996045477014d0 , 0.3d0 /) !POISSON COEFF FOR 2 MATERIALS
  real*8, parameter :: EE(2) = (/ 294.2150271873455d0 , 1.0d0 /) !ELASTICITY MODULI FOR 2 MATERIALS

!  Kronecker's delta
  real*8, parameter, dimension(3,3) :: del = &
        reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0 /),(/3,3/))

contains

  subroutine find_material(X,Fn,IMat)
    real*8, intent(in ):: X(3),Fn(3)
    integer,intent(out):: IMat
! To do...
    Imat = 1
  end subroutine

  !------------------------------------------------------------------------------
  !> Purpose : return elasticity tensor C for linear elasticity
  !! @param[in]  iMat - Index of material (1 or 2)
  !! @param[in]  X    - physical coordinates
  !! @param[out] C    - elasticity tensor
  !! 
  !------------------------------------------------------------------------------
  !
  subroutine getC(iMat, X, C)
    implicit none
  !------------------------------------------------------------------------------
    integer,                   intent(in)  :: iMat
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: C
  !------------------------------------------------------------------------------
    real*8, dimension(3,3,3,3)             :: symm
    real*8                                 :: nu,e,lambda,mu
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  ! Set up the constants
    nu = PN(iMat)
    e  = EE(iMat)
    lambda = (e*nu)/(1.d0-nu-2.d0*nu**2)
    mu     = 0.5d0*e/(1.d0+nu)
    call getSymm(symm)
  !  define the elastictiy tensor
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
      C(i,j,k,l) = lambda * del(i,j)*del(k,l) + 2.d0 * mu * symm(i,j,k,l)
    enddo; enddo; enddo; enddo
  !
  end subroutine getC

  !------------------------------------------------------------------------------
  !> Purpose : return compliance tensor A for linear elasticity
  !! @param[in]  iMat - Index of material (1 or 2)
  !! @param[in]  X    - physical coordinates
  !! @param[out] A    - compliance tensor
  !! 
  !------------------------------------------------------------------------------
  !
  subroutine getA(iMat, X, A)
    implicit none
  !------------------------------------------------------------------------------
    integer,                   intent(in)  :: iMat
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: A
  !------------------------------------------------------------------------------
    real*8, dimension(3,3,3,3)             :: symm
    real*8                                 :: nu,e,const1,const2
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  ! Set up the constants
    nu = PN(iMat)
    e  = EE(iMat)
    const1 = -nu/e
    const2 = (1.d0+nu)/e
    call getSymm(symm)
  !  define the elastictiy tensor
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
      A(i,j,k,l) = const1 * del(i,j)*del(k,l) + const2 * symm(i,j,k,l)
    enddo; enddo; enddo; enddo
  !
  end subroutine getA

  !------------------------------------------------------------------------------
  !> Purpose : return C:C tensor for linear elasticity
  !! @param[in]  iMat - Index of material (1 or 2)
  !! @param[in]  X    - physical coordinates
  !! @param[out] CC   - tensor
  !! 
  !------------------------------------------------------------------------------
  !
  subroutine getCC(iMat, X, CC)
    implicit none
  !------------------------------------------------------------------------------
    integer,                   intent(in)  :: iMat
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: CC
  !------------------------------------------------------------------------------
    real*8, dimension(3,3,3,3)             :: symm
    real*8                                 :: nu,e,lambda,mu,const1,const2
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  ! Set up the constants
    nu = PN(iMat)
    e  = EE(iMat)
    lambda = e * nu / ( 1.d0 - nu - 2.d0*nu**2 )
    mu     = 0.5d0 * e / ( 1.d0 + nu )
    const1 = lambda * (3.d0*lambda + 4.d0*mu)
    const2 = 4.d0 * mu**2
    call getSymm(symm)
  !  define the elastictiy tensor
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
      CC(i,j,k,l) = const1 * del(i,j)*del(k,l) + const2 * symm(i,j,k,l)
    enddo; enddo; enddo; enddo
  !
  end subroutine getCC

  !------------------------------------------------------------------------------
  !> Purpose : return A:A for linear elasticity
  !! @param[in]  iMat - Index of material (1 or 2)
  !! @param[in]  X    - physical coordinates
  !! @param[out] AA   - tensor
  !------------------------------------------------------------------------------
  !
  subroutine getAA(iMat,X, AA)
    implicit none
  !------------------------------------------------------------------------------
    integer,                   intent(in)  :: iMat
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: AA
  !------------------------------------------------------------------------------
    real*8, dimension(3,3,3,3)             :: symm
    real*8                                 :: nu,e,const1,const2
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  ! Set up the constants
    nu = PN(iMat)
    e  = EE(iMat)
    const1 = nu*(nu-2.d0)/(e**2)
    const2 = ((1.d0+nu)/e)**2
    call getSymm(symm)
  !  define the tensor
  do i=1,3; do j=1,3; do k=1,3; do l=1,3
    AA(i,j,k,l) = const1 * del(i,j)*del(k,l) + const2 * symm(i,j,k,l)
  enddo; enddo; enddo; enddo
  !
  end subroutine getAA

  !------------------------------------------------------------------------------
  !> Purpose : return symmetric product tensor
  !!
  !! @param[out] Symm - tensor
  !------------------------------------------------------------------------------
  !
  subroutine getSymm(Symm)
    implicit none
  !------------------------------------------------------------------------------
    real*8,dimension(3,3,3,3), intent(out) :: Symm
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  !  define the tensor
  do i=1,3; do j=1,3; do k=1,3; do l=1,3
    Symm(i,j,k,l) = 0.5d0*(del(i,k)*del(j,l) + del(i,l)*del(j,k))
  enddo; enddo; enddo; enddo
  !
  end subroutine getSymm

end module isotropic_elast_material_compo
