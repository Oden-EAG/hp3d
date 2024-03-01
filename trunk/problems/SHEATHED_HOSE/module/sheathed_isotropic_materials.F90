!--------------------------------------------------------------------------------
!> @brief define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module sheathed_isotropic_materials
!
  implicit none
!
!  G E O M E T R Y
  real(8), parameter :: R_inside  = 0.5d0
  real(8), parameter :: R_middle  = 0.99d0
  real(8), parameter :: R_outside = 1.0d0
  real(8), parameter :: X_1       = 0.0d0
  real(8), parameter :: X_2       = 1.0d0
!
!  M A T E R I A L   D A T A
  real(8), parameter :: E_steel   = 2.0d4
  real(8), parameter :: E_rubber  = 1.0d0
  real(8), parameter :: NU_steel  = 0.285d0
  real(8), parameter :: NU_rubber = 0.5d0
!
  real(8), parameter :: LAMBDA_steel  = E_steel*NU_steel/((1+NU_steel)*(1-2*NU_steel))
  real(8), parameter :: MU_steel      = 0.5d0*E_steel/(1+NU_steel)
  real(8), parameter :: MU_rubber     = 0.5d0*E_rubber/(1+NU_rubber)
!
  real(8), parameter :: RHO   = 0.d0
  real(8), parameter :: OMEG  = 0.d0
!
!  Pressure
  real(8), parameter :: P_inner =-1.d-1
  real(8), parameter :: P_outer = 0.d0

!  Kronecker's delta
  real(8), parameter, dimension(3,3) :: del =  &
      reshape((/ 1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0 /),(/3,3/))

contains

  !------------------------------------------------------------------------------
  !> @brief return stiffness tensor for linear elasticity
  !!
  !> @param[in]  X    - physical coordinates
  !> @param[in]  Dom  - Domain number
  !> @param[out] C    - stiffness tensor
  !------------------------------------------------------------------------------
  !
  subroutine getC(X,Dom, C)
    implicit none
  !------------------------------------------------------------------------------
    real(8),dimension(3),       intent(in)  :: X
    integer,                   intent(in)  :: Dom
    real(8),dimension(3,3,3,3), intent(out) :: C
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
    real(8)  :: mu,lambda
  !------------------------------------------------------------------------------
  !
  if (Dom.eq.1) then
    !  STEEL
    mu     = MU_steel
    lambda = LAMBDA_steel
  elseif (Dom.eq.2) then
    !  RUBBER
    mu     = MU_rubber
    if (NU_rubber.eq.0.5d0) then
      write(*,*) 'Error: cannot define C in this way for an incompressible material'
      stop 1
    endif
  endif

  !  define the stiffness tensor
  do i=1,3; do j=1,3; do k=1,3; do l=1,3
    C(i,j,k,l) = mu*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
               + lambda*del(i,j)*del(k,l)
  enddo; enddo; enddo; enddo
  !
  end subroutine getC

  !------------------------------------------------------------------------------
  !> @brief return compliance tensor for linear elasticity
  !!
  !> @param[in]  X    - physical coordinates
  !> @param[in]  Dom  - Domain number
  !> @param[out] A    - compliance tensor
  !------------------------------------------------------------------------------
  !
  subroutine getA(X,Dom, A)
    implicit none
  !------------------------------------------------------------------------------
    real(8),dimension(3),       intent(in)  :: X
    integer,                   intent(in)  :: Dom
    real(8),dimension(3,3,3,3), intent(out) :: A
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
    real(8)  :: nu,e
  !------------------------------------------------------------------------------
  !
  if (Dom.eq.1) then
    !  STEEL
    nu = NU_steel
    e  = E_steel
  elseif (Dom.eq.2) then
    !  RUBBER
    nu = NU_rubber
    e  = E_rubber
  else
    write(*,*) 'getA: invalid Dom = ', Dom
    stop
  endif

  !  define the compliance tensor
  do i=1,3; do j=1,3; do k=1,3; do l=1,3
    A(i,j,k,l) = 0.5d0*(1+nu)/e*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
               - nu/e*del(i,j)*del(k,l)
  enddo; enddo; enddo; enddo
  !
  end subroutine getA

  !------------------------------------------------------------------------------
  !> @brief return A:A for linear elasticity
  !!
  !> @param[in]  X    - physical coordinates
  !> @param[out] AA   - tensor
  !------------------------------------------------------------------------------
  !
  subroutine getAA(X, AA)
    implicit none
  !------------------------------------------------------------------------------
    real(8),dimension(3),       intent(in)  :: X
    real(8),dimension(3,3,3,3), intent(out) :: AA
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
    real(8) :: mag
  !------------------------------------------------------------------------------
  !
  mag = X(2)**2+X(3)**2
  mag = dsqrt(mag)
  !  define the compliance tensor
  if (mag.gt.R_middle) then
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
      AA(i,j,k,l) = 0.5d0*((1+NU_steel)/E_steel)**2*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
                  + NU_steel*(NU_steel-2)/E_steel**2*del(i,j)*del(k,l)
    enddo; enddo; enddo; enddo
  elseif (mag.le.R_middle) then
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
      AA(i,j,k,l) = 0.5d0*((1+NU_rubber)/E_rubber)**2*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
                  + NU_rubber*(NU_rubber-2)/E_rubber**2*del(i,j)*del(k,l)
    enddo; enddo; enddo; enddo
  endif
  !
  end subroutine getAA

  !------------------------------------------------------------------------------
  !> @brief return symmetric product tensor
  !!
  !> @param[out] Symm - tensor
  !------------------------------------------------------------------------------
  !
  subroutine getSymm(Symm)
    implicit none
  !------------------------------------------------------------------------------
    real(8),dimension(3,3,3,3), intent(out) :: Symm
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

  !------------------------------------------------------------------------------
  !> @brief return skew-symmetric product tensor
  !!
  !> @param[out] Skew - tensor
  !------------------------------------------------------------------------------
  !
  subroutine getSkew(Skew)
    implicit none
  !------------------------------------------------------------------------------
    real(8),dimension(3,3,3,3), intent(out) :: Skew
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  !  define the tensor
  do i=1,3; do j=1,3; do k=1,3; do l=1,3
    Skew(i,j,k,l) = 0.5d0*(del(i,k)*del(j,l) - del(i,l)*del(j,k))
  enddo; enddo; enddo; enddo
  !
  end subroutine getSkew

  !------------------------------------------------------------------------------
  !> @brief return skew-symmetric product tensor
  !!
  !> @param[in]  X    - physical coordinates
  !> @param[out] Skew - tensor
  !------------------------------------------------------------------------------
  !
  subroutine getWeightedSkew(X, Skew)
    implicit none
  !------------------------------------------------------------------------------
    real(8),dimension(3),       intent(in)  :: X
    real(8),dimension(3,3,3,3), intent(out) :: Skew
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
    real(8) :: mag
  !------------------------------------------------------------------------------
  !
  mag = X(2)**2+X(3)**2
  mag = dsqrt(mag)
  !  define the tensor
  if (mag.gt.R_middle) then
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
      Skew(i,j,k,l) = 0.5d0/E_steel**2*(del(i,k)*del(j,l) - del(i,l)*del(j,k))
    enddo; enddo; enddo; enddo
  elseif (mag.le.R_middle) then
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
      Skew(i,j,k,l) = 0.5d0/E_rubber**2*(del(i,k)*del(j,l) - del(i,l)*del(j,k))
    enddo; enddo; enddo; enddo
  endif
  !
  end subroutine getWeightedSkew

  !------------------------------------------------------------------------------
  !> @brief return Young's modulus
  !!
  !> @param[in]  X    - physical coordinates
  !> @param[out] E    - Young's modulus
  !------------------------------------------------------------------------------
  !
  subroutine getYoungsModulus(X, E)
    implicit none
  !------------------------------------------------------------------------------
    real(8),dimension(3), intent(in)  :: X
    real(8),              intent(out) :: E
  !------------------------------------------------------------------------------
    real(8) :: mag
  !------------------------------------------------------------------------------
  !
  mag = X(2)**2+X(3)**2
  mag = dsqrt(mag)

  if (mag.gt.R_middle) then
    E = E_steel
  elseif (mag.le.R_middle) then
    E = E_rubber
  endif
  !
  end subroutine getYoungsModulus




end module sheathed_isotropic_materials
