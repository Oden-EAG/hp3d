!--------------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module primal_module

  use parametersDPG
  use parameters
!
!  ...MATRICES
!     stiffnes matrices (and copies) for the enriched test space
!       real*8, dimension(3*MAXbrickHH,3*MAXbrickH) :: EnrField,EnrFieldc
! !$OMP THREADPRIVATE (EnrField, EnrFieldc)
!       real*8, dimension(3*MAXbrickHH,3*MAXbrickV) :: EnrTrace,EnrTracec
! !$OMP THREADPRIVATE (EnrTrace, EnrTracec)
! !     Gram matrix for the local Riesz matrix in LAPACK format
!       real*8, dimension(3*MAXbrickHH*(3*MAXbrickHH+1)/2) :: Gram
! !$OMP THREADPRIVATE (Gram)

contains
  !------------------------------------------------------------------------------
  !> Purpose : return operator defining the semi inner product
  !!
  !! @param[out] T - tensor defining the semi inner product used in the Gram
  !!                 matrix construction
  !------------------------------------------------------------------------------
  !
  subroutine getSemiIPTensor(X, T)
    use isotropic_elast_material
    use common_prob_data
    use parameters, only : ZERO
    implicit none
  !------------------------------------------------------------------------------
    real*8, dimension(3),       intent(in)  :: X
    real*8, dimension(3,3,3,3), intent(out) :: T
  !------------------------------------------------------------------------------
    real*8, dimension(3,3,3,3) :: C
    integer :: i,j,k,l,m,n
  !
  !  Kronecker's delta
    real*8, parameter, dimension(3,3) :: delta = &
          reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0 /),(/3,3/))
  !------------------------------------------------------------------------------
    T = ZERO

    select case(TEST_NORM)
    case(ADJOINT_GRAPH)
      !   Adjoint-Graph
      call getC(X, C)
      do i=1,3; do j=1,3; do k=1,3; do l=1,3
        do m=1,3; do n=1,3
          T(i,j,k,l) = T(i,j,k,l)  &
                     + C(m,n,i,j)*C(m,n,k,l)
        enddo; enddo
      enddo; enddo; enddo; enddo
    case(MATHEMATICIANS)
      ! H1
      do i=1,3; do j=1,3; do k=1,3; do l=1,3
        T(i,j,k,l) = delta(i,k)*delta(j,l)
     enddo; enddo; enddo; enddo
    case(H1SYMM)
      ! H1_symm
     do i=1,3; do j=1,3; do k=1,3; do l=1,3
        T(i,j,k,l) = 0.5d0*(delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k))
     enddo; enddo; enddo; enddo
    end select

  end subroutine getSemiIPTensor

end module primal_module
