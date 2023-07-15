!------------------------------------------------------------------------------
!> Purpose     Compute source term
!!
!> @param[in]  Mdle  - element (middle node) number
!> @param[in]  X     - physical coordinates
!!
!> @param[out] Fval  - rhs
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine getf(Mdle,X, Fval)
!
      use control   , only: NEXACT
      use assembly  , only: NR_RHS
!
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(in)  :: X(3)
      real(8), intent(out) :: Fval(1:3,NR_RHS)
!
!  ...solution variables
      real(8) :: u(3)
      real(8) :: gradu(3,3)
      real(8) :: epsilon(3,3)
      real(8) :: sigma(3,3)
      real(8) :: divsigma(3)
!
      integer :: iload
!
!------------------------------------------------------------------------------
!
!  ...initialize source terms
      Fval(1:3,1:NR_RHS) = 0.d0
!
      select case(NEXACT)
!  ...unknown exact solution
      case(0)
         continue
!
!  ...known exact solution
      case(1)
!     ...compute exact soution
         call elast_solution(X, u,gradu,epsilon,sigma,divsigma)
         do iload=1,NR_RHS
            Fval(1:3,iload) =-divsigma
         enddo
!
!  ...homogenous RHS
      case(2)
         continue
      end select
!
   end subroutine getf




!------------------------------------------------------------------------------
!> @brief      Neumann load
!!
!> @param[in]  Mdle  - element (middle node) number
!> @param[in]  Ibc   - boundary conditions
!> @param[in]  X     - physical coordinates
!> @param[in]  Rn    - outward unit vector
!> @param[out] Gval  - Neumann loads
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine getg(Mdle,Ibc,X,Rn, Gval)
!
      use control    , only: NEXACT
      use assembly   , only: NR_RHS
!
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(in)  :: X(3),Rn(3)
      integer, intent(in)  :: Ibc
      real(8), intent(out) :: Gval(3,NR_RHS)
!
!  ...solution variables
      real*8, dimension(3)   :: u
      real*8, dimension(3,3) :: gradu
      real*8, dimension(3,3) :: epsilon
      real*8, dimension(3,3) :: sigma
      real*8, dimension(3)   :: divsigma
!
      integer :: iload,j
!
!------------------------------------------------------------------------------
!
!  ...initialize source terms
      Gval(1:3,1:NR_RHS) = 0.d0
!
      select case(NEXACT)
!  ...unknown exact solution
      case(0)
         continue
!
!  ...known exact solution
      case(1,2)
!
!     ...compute exact soution
         call elast_solution(X, u,gradu,epsilon,sigma,divsigma)
!
         select case(Ibc)
!     ...Neumann BC
         case(2)
!
!        ...loop over rhs's
            do iload=1,NR_RHS
!              g_i = sigma_ij n_j
               do j=1,3
                  Gval(1:3,iload) = Gval(1:3,iload) + sigma(1:3,j)*Rn(j)
               enddo
            enddo
!
!     ...interface BC
         case(3)
!           TODO: implement
         end select
      endselect
!
   endsubroutine getg
