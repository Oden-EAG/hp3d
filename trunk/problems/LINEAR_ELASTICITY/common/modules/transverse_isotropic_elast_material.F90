!--------------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module transverse_isotropic_elast_material
  !
  !  M A T E R I A L   D A T A
  real*8, parameter :: E1   = 1.d0
  ! real*8, parameter :: E3   = 1.d1
  real*8, parameter :: E3   = 1.d2
  ! real*8, parameter :: E3   = 1.d3
  ! real*8, parameter :: E3   = 1.d4
  real*8, parameter :: NU12 = 0.4d0
  real*8, parameter :: NU31 = 0.3d0
  real*8, parameter :: G31  = 0.4d0
!
  real*8, parameter :: RHO  = 0.d0
  real*8, parameter :: OMEG = 0.d0
!
  real*8, parameter :: E1byE3 = E1/E3
  real*8, parameter :: const1 = 1.d0 - NU12 - 2.d0*NU31**2*E1byE3
  real*8, parameter :: const2 = E1/(1.d0 + NU12)
  real*8, parameter :: C11  = (1.d0 - NU31**2*E1byE3)*const2/const1 
  real*8, parameter :: C12  = (NU12 + NU31**2*E1byE3)*const2/const1
  real*8, parameter :: C13  = NU31*E1/const1
  real*8, parameter :: C33  = (1.d0 - NU12)*E3/const1
  real*8, parameter :: C44  = G31
  !  (N11-N12)/2
  real*8, parameter :: C66  =  const2/2.d0
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
  C=0.d0
  !
  C(1,1,1,1) = C11
  C(2,2,2,2) = C11
  C(3,3,3,3) = C33
  C(1,1,2,2) = C12; C(2,2,1,1) = C12
  C(1,1,3,3) = C13; C(3,3,1,1) = C13
  C(2,2,3,3) = C13; C(3,3,2,2) = C13
  !
  C(2,3,2,3) = C44; C(2,3,3,2) = C44; C(3,2,2,3) = C44; C(3,2,3,2) = C44
  C(3,1,3,1) = C44; C(3,1,1,3) = C44; C(1,3,3,1) = C44; C(1,3,1,3) = C44
  C(1,2,1,2) = C66; C(1,2,2,1) = C66; C(2,1,1,2) = C66; C(2,1,2,1) = C66

  ! write(*,*) 'lambda = ', E1*NU12/((1.d0+NU12)*(1.d0-2.d0*NU12))
  ! write(*,*) 'mu = ', E1/(1.d0+NU12)/2.d0

  ! do i=1,3; do j=1,3; do k=1,3; do l=1,3
  !   write(*,*) ' C(',i,',',j,',',k,',',l,') = ', C(i,j,k,l)
  ! enddo; enddo; enddo; enddo

  ! call pause

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
    ! real*8,dimension(3,3,3,3) :: C
    ! real*8,dimension(3,3,3,3) :: tmp
  !------------------------------------------------------------------------------
    integer :: i,j,k,l !,m,n
    real*8 :: const
  !------------------------------------------------------------------------------
  !
  A=0.d0
  !
  A(1,1,1,1) = C11*C33 - C13**2
  A(2,2,2,2) = C11*C33 - C13**2
  A(3,3,3,3) = C11**2 - C12**2
  A(1,1,2,2) = C13**2 - C12*C33; A(2,2,1,1) = A(1,1,2,2)
  A(1,1,3,3) = (C12-C11)*C13;    A(3,3,1,1) = A(1,1,3,3)
  A(2,2,3,3) = (C12-C11)*C13;    A(3,3,2,2) = A(2,2,3,3)

  const = (C11-C12)*((C11+C12)*C33 - 2.d0*C13**2)
  A = A/const
  !
  A(2,3,2,3) = 0.25d0/C44; A(2,3,3,2) = 0.25d0/C44; A(3,2,2,3) = 0.25d0/C44; A(3,2,3,2) = 0.25d0/C44
  A(3,1,3,1) = 0.25d0/C44; A(3,1,1,3) = 0.25d0/C44; A(1,3,3,1) = 0.25d0/C44; A(1,3,1,3) = 0.25d0/C44
  A(1,2,1,2) = 0.25d0/C66; A(1,2,2,1) = 0.25d0/C66; A(2,1,1,2) = 0.25d0/C66; A(2,1,2,1) = 0.25d0/C66

  ! call getC(X, C) 

  ! tmp = 0.d0

  ! do i=1,3; do j=1,3; do k=1,3; do l=1,3
  !   do m=1,3; do n=1,3
  !     tmp(i,j,k,l) = tmp(i,j,k,l) + A(i,j,m,n)*C(m,n,k,l)
  !   enddo; enddo
  !   write(*,*) ' tmp(',i,',',j,',',k,',',l,') = ', tmp(i,j,k,l)
  !   ! write(*,*) ' C(',i,',',j,',',k,',',l,') = ', C(i,j,k,l)
  !   ! write(*,*) ' A(',i,',',j,',',k,',',l,') = ', A(i,j,k,l)
  ! enddo; enddo; enddo; enddo

  ! do i=1,3; do j=1,3; do k=1,3; do l=1,3
  !   write(*,*) ' A(',i,',',j,',',k,',',l,') = ', A(i,j,k,l)
  ! enddo; enddo; enddo; enddo

  ! call pause



  !
  end subroutine getA

  ! !------------------------------------------------------------------------------
  ! !> Purpose : return C:C for linear elasticity
  ! !!
  ! !! @param[in]  X    - physical coordinates
  ! !! @param[out] CC   - tensor
  ! !------------------------------------------------------------------------------
  ! !
  ! subroutine getCC(X, CC)
  !   implicit none
  ! !------------------------------------------------------------------------------
  !   real*8,dimension(3),       intent(in)  :: X
  !   real*8,dimension(3,3,3,3), intent(out) :: CC
  ! !------------------------------------------------------------------------------
  !   integer :: i,j,k,l
  ! !------------------------------------------------------------------------------
  ! !
  ! !  define the tensor
  ! do i=1,3; do j=1,3; do k=1,3; do l=1,3
  !   CC(i,j,k,l) = 2.0d0*MU**2*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
  !               + LAMBDA*(4.0d0*MU+3.0d0*LAMBDA)*del(i,j)*del(k,l)
  ! enddo; enddo; enddo; enddo
  ! !
  ! end subroutine getCC

  ! !------------------------------------------------------------------------------
  ! !> Purpose : return A:A for linear elasticity
  ! !!
  ! !! @param[in]  X    - physical coordinates
  ! !! @param[out] AA   - tensor
  ! !------------------------------------------------------------------------------
  ! !
  ! subroutine getAA(X, AA)
  !   implicit none
  ! !------------------------------------------------------------------------------
  !   real*8,dimension(3),       intent(in)  :: X
  !   real*8,dimension(3,3,3,3), intent(out) :: AA
  ! !------------------------------------------------------------------------------
  !   integer :: i,j,k,l
  ! !------------------------------------------------------------------------------
  ! !
  ! !  define the tensor
  ! do i=1,3; do j=1,3; do k=1,3; do l=1,3
  !   AA(i,j,k,l) = 1.d0/(4.d0*MU**2)*( 0.5d0*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
  !               - LAMBDA*(4.d0*LAMBDA+3.d0*MU)/(2.d0*MU+3.d0*LAMBDA)**2*del(i,j)*del(k,l) )
  ! enddo; enddo; enddo; enddo
  ! !
  ! end subroutine getAA

  ! !------------------------------------------------------------------------------
  ! !> Purpose : return symmetric product tensor
  ! !!
  ! !! @param[out] Symm - tensor
  ! !------------------------------------------------------------------------------
  ! !
  ! subroutine getSymm(Symm)
  !   implicit none
  ! !------------------------------------------------------------------------------
  !   real*8,dimension(3,3,3,3), intent(out) :: Symm
  ! !------------------------------------------------------------------------------
  !   integer :: i,j,k,l
  ! !------------------------------------------------------------------------------
  ! !
  ! !  define the tensor
  ! do i=1,3; do j=1,3; do k=1,3; do l=1,3
  !   Symm(i,j,k,l) = 0.5d0*(del(i,k)*del(j,l) + del(i,l)*del(j,k))
  ! enddo; enddo; enddo; enddo
  ! !
  ! end subroutine getSymm

end module transverse_isotropic_elast_material
