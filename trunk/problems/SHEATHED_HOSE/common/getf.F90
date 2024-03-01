!------------------------------------------------------------------------------
!> @brief source term
!!
!> @param[in]  Mdle  - element (middle node) number
!> @param[in]  X     - physical coordinates
!> @param[out] Fval  - rhs
!------------------------------------------------------------------------------
!
subroutine getf(Mdle,X, Fval)
  use control   , only: NEXACT
  use assembly  , only: NR_RHS
  implicit none
!------------------------------------------------------------------------------
  integer,                       intent(in)  :: Mdle
  real(8), dimension(3),          intent(in)  :: X
  real(8), dimension(1:3,NR_RHS), intent(out) :: Fval
!------------------------------------------------------------------------------
! solution variables
  real(8), dimension(3)   :: u
  real(8), dimension(3,3) :: gradu
  real(8), dimension(3,3) :: epsilon
  real(8), dimension(3,3) :: sigma
  real(8), dimension(3)   :: divsigma
! counters
  integer :: iload
!------------------------------------------------------------------------------
!
! initialize source terms
  Fval(1:3,1:NR_RHS) = 0.d0
!
  select case(NEXACT)
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                     |
!==============================================================================
  case(0)
!
!==============================================================================
!  KNOWN EXACT SOLUTION                                                       |
!==============================================================================
  case(1)
!   compute exact soution
    call elast_solution(X, u,gradu,epsilon,sigma,divsigma)

    do iload=1,NR_RHS
      Fval(1:3,iload) =-divsigma
    enddo

!==============================================================================
!  KNOWN EXACT SOLUTION , HOMOGENEOUS RHS                                     |
!==============================================================================
  case(2)
!
  endselect
!
endsubroutine getf
!
!
!
!------------------------------------------------------------------------------
!> @brief Neumann load
!!
!> @param[in]  Mdle  - element (middle node) number
!> @param[in]  Ibc   - boundary conditions
!> @param[in]  X     - physical coordinates
!> @param[in]  Rn    - outward unit vector
!> @param[out] Gval  - Neumann loads
!------------------------------------------------------------------------------
!
subroutine getg(Mdle,Ibc,X,Rn, Gval)
  use control    , only: NEXACT
  use assembly   , only: NR_RHS
  implicit none
!------------------------------------------------------------------------------
  integer,                     intent(in)  :: Mdle
  real(8), dimension(3),        intent(in)  :: X,Rn
  integer,                     intent(in)  :: Ibc
  real(8), dimension(3,NR_RHS), intent(out) :: Gval
!------------------------------------------------------------------------------
! solution variables
  real(8), dimension(3)   :: u
  real(8), dimension(3,3) :: gradu
  real(8), dimension(3,3) :: epsilon
  real(8), dimension(3,3) :: sigma
  real(8), dimension(3)   :: divsigma
! counters
  integer :: iload,j
!------------------------------------------------------------------------------
!
! initialize source terms
  Gval(1:3,1:NR_RHS) = 0.d0
!
  select case(NEXACT)
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                      |
!==============================================================================
  case(0)
!
!==============================================================================
!  KNOWN EXACT SOLUTION                                                       |
!==============================================================================
  case(1,2)
!
!   compute exact soution
    call elast_solution(X, u,gradu,epsilon,sigma,divsigma)
!
    select case(Ibc)
!   NEUMANN BC
    case(2)
!
!     loop over rhs's
      do iload=1,NR_RHS
!       g_i = sigma_ij n_j
        do j=1,3
          Gval(1:3,iload) = Gval(1:3,iload) + sigma(1:3,j)*Rn(j)
        enddo
      enddo
!   INTERFACE BC
    case(3)
      !  nothing here yet
    endselect
  endselect
!
endsubroutine getg
