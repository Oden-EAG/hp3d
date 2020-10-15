!--------------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module isotropic_elast_material
  !
  !  M A T E R I A L   D A T A
  !  Lame parameters for the singular solution
  real*8, parameter :: LAMBDA = 123.d0
  real*8, parameter :: MU     = 79.3d0
  ! real*8, parameter :: LAMBDA = 66.16380367543577d0
  ! real*8, parameter :: MU     = 99.3d0
  !  Custom Lame parameters
  ! real*8, parameter :: LAMBDA = 1.d0
  ! real*8, parameter :: MU     = 1.d0
!
  real*8, parameter :: RHO   = 0.d0
  real*8, parameter :: OMEG  = 0.d0

  ! real*8, parameter :: NU(2) = (/ 0.3d0 , 0.3d0 /)
  ! real*8, parameter :: EE(2) = (/ 1.0d0 , 1.0d0 /)

!  Kronecker's delta
   real*8, parameter, dimension(3,3) :: del = &
        reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0 /),(/3,3/))

contains

  !------------------------------------------------------------------------------
  !> Purpose : return stiffness tensor for linear elasticity
  !!
  !! @param[in]  X    - physical coordinates
  !! @param[out] C    - stiffness tensor
  !------------------------------------------------------------------------------
  !
  subroutine getC(X, C)
    implicit none
  !------------------------------------------------------------------------------
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: C
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  !  define the stiffness tensor
     do i=1,3; do j=1,3; do k=1,3; do l=1,3
        C(i,j,k,l) = MU*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
                   + LAMBDA*del(i,j)*del(k,l)
     enddo; enddo; enddo; enddo
  !
  end subroutine getC

  !------------------------------------------------------------------------------
  !> Purpose : return compliance tensor for linear elasticity
  !!
  !! @param[in]  X    - physical coordinates
  !! @param[out] A    - compliance tensor
  !------------------------------------------------------------------------------
  !
  subroutine getA(X, A)
    implicit none
  !------------------------------------------------------------------------------
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: A
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  !  define the elastictiy tensor
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
      A(i,j,k,l) = 1.d0/(2.d0*MU)*( 0.5d0*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
                 - LAMBDA/(2.d0*MU+3.d0*LAMBDA)*del(i,j)*del(k,l) )
    enddo; enddo; enddo; enddo
  !
  end subroutine getA

  !------------------------------------------------------------------------------
  !> Purpose : return C:C for linear elasticity
  !!
  !! @param[in]  X    - physical coordinates
  !! @param[out] CC   - tensor
  !------------------------------------------------------------------------------
  !
  subroutine getCC(X, CC)
    implicit none
  !------------------------------------------------------------------------------
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: CC
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  !  define the tensor
  do i=1,3; do j=1,3; do k=1,3; do l=1,3
    CC(i,j,k,l) = 2.0d0*MU**2*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
                + LAMBDA*(4.0d0*MU+3.0d0*LAMBDA)*del(i,j)*del(k,l)
  enddo; enddo; enddo; enddo
  !
  end subroutine getCC

  !------------------------------------------------------------------------------
  !> Purpose : return A:A for linear elasticity
  !!
  !! @param[in]  X    - physical coordinates
  !! @param[out] AA   - tensor
  !------------------------------------------------------------------------------
  !
  subroutine getAA(X, AA)
    implicit none
  !------------------------------------------------------------------------------
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: AA
  !------------------------------------------------------------------------------
    integer :: i,j,k,l
  !------------------------------------------------------------------------------
  !
  !  define the tensor
  do i=1,3; do j=1,3; do k=1,3; do l=1,3
    AA(i,j,k,l) = 1.d0/(4.d0*MU**2)*( 0.5d0*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
                - LAMBDA*(4.d0*LAMBDA+3.d0*MU)/(2.d0*MU+3.d0*LAMBDA)**2*del(i,j)*del(k,l) )
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

end module isotropic_elast_material
