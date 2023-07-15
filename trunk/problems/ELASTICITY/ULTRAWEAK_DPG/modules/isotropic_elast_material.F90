!--------------------------------------------------------------------------------
!> @brief      Define material tensors for isotropic elastic material
!!
!> @date       July 2023
!--------------------------------------------------------------------------------
module isotropic_elast_material
!
!..Lame parameters
   real(8), parameter :: LAMBDA = 1.d0
   real(8), parameter :: MU     = 1.d0
!..for singular solution:
   !real(8), parameter :: LAMBDA = 123.d0
   !real(8), parameter :: MU     = 79.3d0
   !real(8), parameter :: LAMBDA = 66.16380367543577d0
   !real(8), parameter :: MU     = 99.3d0
!
   real(8), parameter :: RHO   = 0.d0
   real(8), parameter :: OMEG  = 0.d0
!
!..Kronecker's delta
   real(8), parameter, dimension(3,3) :: del = &
        reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0 /),(/3,3/))


contains



!------------------------------------------------------------------------------
!> @brief      Compute stiffness tensor for linear elasticity
!!
!> @param[in]  X    - physical coordinates
!> @param[out] C    - stiffness tensor
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine getC(X, C)
!
      implicit none
!
    real*8,dimension(3),       intent(in)  :: X
    real*8,dimension(3,3,3,3), intent(out) :: C
!
    integer :: i, j, k, l
!
!------------------------------------------------------------------------------
!
!  ...define the stiffness tensor
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  C(i,j,k,l) = MU*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
                            + LAMBDA*del(i,j)*del(k,l)
               enddo
            enddo
         enddo
      enddo
!
   end subroutine getC




!------------------------------------------------------------------------------
!> @brief      Compute compliance tensor for linear elasticity
!!
!! @param[in]  X    - physical coordinates
!! @param[out] A    - compliance tensor
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine getA(X, A)
!
      implicit none
!
    real(8), intent(in)  :: X(3)
    real(8), intent(out) :: A(3,3,3,3)
!
    integer :: i, j, k, l
!
!------------------------------------------------------------------------------
!
!  ...define the elastictiy tensor
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  A(i,j,k,l) = 1.d0/(2.d0*MU)                                    &
                             * ( 0.5d0*(del(i,k)*del(j,l) + del(i,l)*del(j,k))   &
                               - LAMBDA/(2.d0*MU + 3.d0*LAMBDA)*del(i,j)*del(k,l) )
               enddo
            enddo
         enddo
      enddo
!
   end subroutine getA




!------------------------------------------------------------------------------
!> @brief      Return C:C for linear elasticity
!!
!! @param[in]  X    - physical coordinates
!! @param[out] CC   - tensor
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine getCC(X, CC)
!
      implicit none
!
      real(8), intent(in)  :: X(3)
      real(8), intent(out) :: CC(3,3,3,3)
!
      integer :: i, j, k, l
!
!------------------------------------------------------------------------------
!
!  ...define the tensor
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  CC(i,j,k,l) = 2.d0*MU**2                                       &
                              * (del(i,k)*del(j,l) + del(i,l)*del(j,k))          &
                              + LAMBDA*(4.0d0*MU + 3.0d0*LAMBDA)*del(i,j)*del(k,l)
               enddo
            enddo
         enddo
      enddo
!
   end subroutine getCC





!------------------------------------------------------------------------------
!> @brief      Return A:A for linear elasticity
!!
!! @param[in]  X    - physical coordinates
!! @param[out] AA   - tensor
!!
!> @date       July 2023
!------------------------------------------------------------------------------
!
   subroutine getAA(X, AA)
!
      implicit none
!
      real(8), intent(in)  :: X(3)
      real(8), intent(out) :: AA(3,3,3,3)
!
      integer :: i, j, k, l
!
!------------------------------------------------------------------------------
!
!  ...define the tensor
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  AA(i,j,k,l) = 1.d0/(4.d0*MU**2)                                &
                              * ( 0.5d0*(del(i,k)*del(j,l) + del(i,l)*del(j,k))  &
                                - LAMBDA*(4.d0*LAMBDA+3.d0*MU)                   &
                                  / (2.d0*MU + 3.d0*LAMBDA)**2*del(i,j)*del(k,l) )
               enddo
            enddo
         enddo
      enddo
!
   end subroutine getAA





!------------------------------------------------------------------------------
!> @brief      Return symmetric product tensor
!!
!> @param[out] Symm - tensor
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine getSymm(Symm)
!
      implicit none
!
      real(8), intent(out) :: Symm(3,3,3,3)
!
      integer :: i, j, k, l
!
!------------------------------------------------------------------------------
!
!  ...define the tensor
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Symm(i,j,k,l) = 0.5d0*(del(i,k)*del(j,l) + del(i,l)*del(j,k))
               enddo
            enddo
         enddo
      enddo
!
   end subroutine getSymm

end module isotropic_elast_material
