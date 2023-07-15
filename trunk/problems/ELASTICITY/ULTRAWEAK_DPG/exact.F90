!------------------------------------------------------------------------------
!> @brief      Compute exact solution at point
!!
!> @param[in]  X       - a point in physical space
!> @param[in]  Icase   - node case (specifies what variables are supported)
!!
!> @param[out] zvalH   - value of the H1 solution
!> @param[out] ZdvalH  - corresponding first derivatives
!> @param[out] Zd2valH - corresponding second derivatives
!> @param[out] ZvalE   - value of the H(curl) solution
!> @param[out] ZdvalE  - corresponding first derivatives
!> @param[out] Zd2valE - corresponding second derivatives
!> @param[out] ZvalV   - value of the H(div) solution
!> @param[out] ZdvalV  - corresponding first derivatives
!> @param[out] Zd2valV - corresponding second derivatives
!> @param[out] ZvalQ   - value of the H(div) solution
!> @param[out] ZdvalQ  - corresponding first derivatives
!> @param[out] Zd2valQ - corresponding second derivatives
!!
!> @date       July 2023
!------------------------------------------------------------------------------
!
   subroutine exact(X,Icase, ValH,DvalH,D2valH, &
                             ValE,DvalE,D2valE, &
                             ValV,DvalV,D2valV, &
                             ValQ,DvalQ,D2valQ)
!
      use data_structure3D
!
      implicit none
!
!------------------------------------------------------------------------------
      real(8), intent(in)  :: X(3)
      integer, intent(in)  :: Icase
      real(8), intent(out) ::   ValH(  MAXEQNH    )
      real(8), intent(out) ::  DvalH(  MAXEQNH,3  )
      real(8), intent(out) :: D2valH(  MAXEQNH,3,3)
      real(8), intent(out) ::   ValE(3,MAXEQNE    )
      real(8), intent(out) ::  DvalE(3,MAXEQNE,3  )
      real(8), intent(out) :: D2valE(3,MAXEQNE,3,3)
      real(8), intent(out) ::   ValV(3,MAXEQNV    )
      real(8), intent(out) ::  DvalV(3,MAXEQNV,3  )
      real(8), intent(out) :: D2valV(3,MAXEQNV,3,3)
      real(8), intent(out) ::   ValQ(  MAXEQNQ    )
      real(8), intent(out) ::  DvalQ(  MAXEQNQ,3  )
      real(8), intent(out) :: D2valQ(  MAXEQNQ,3,3)
!
      real(8) :: u(3)
      real(8) :: gradu(3,3)
      real(8) :: epsilon(3,3)
      real(8) :: sigma(3,3)
      real(8) :: divsigma(3)
!
!------------------------------------------------------------------------------
!
!  ...initialize exact solution
      ValH = 0.d0 ; DvalH = 0.d0 ; D2valH = 0.d0
      ValE = 0.d0 ; DvalE = 0.d0 ; D2valE = 0.d0
      ValV = 0.d0 ; DvalV = 0.d0 ; D2valV = 0.d0
      ValQ = 0.d0 ; DvalQ = 0.d0 ; D2valQ = 0.d0
!
      call elast_solution(X, u,gradu,epsilon,sigma,divsigma)
!
      ValH(:) = u(:)
      ValV(:,:) = sigma(:,:)  !  stress is symmetric so this works
      ValQ(1:3)  = u(:)
!
!          ( sigma1   sigma2  sigma4 )
!  sigma = ( sigma2   sigma3  sigma5 )
!          ( sigma4   sigma5  sigma6 )
!
      ValQ(4)  = sigma(1,1)
      ValQ(5)  = sigma(1,2)
      ValQ(6)  = sigma(2,2)
      ValQ(7)  = sigma(1,3)
      ValQ(8)  = sigma(2,3)
      ValQ(9)  = sigma(3,3)
!
   end subroutine exact
