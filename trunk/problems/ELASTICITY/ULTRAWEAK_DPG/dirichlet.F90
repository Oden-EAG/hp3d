!--------------------------------------------------------------------
!> Purpose : calculate dirichlet boundary condition
!! @param[in]  X      - physical coordinates of a point
!! @param[in]  Icase  - node case
!!
!! @param[out] ValH   - value of the H1 solution
!! @param[out] DvalH  - H1 corresponding first derivatives
!! @param[out] ValE   - value of the H(curl) solution
!! @param[out] DvalE  - H(curl) corresponding first derivatives
!! @param[out] ValV   - value of the H(div) solution
!! @param[out] DvalV  - H(div) corresponding first derivatives
!--------------------------------------------------------------------
   subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!
      use control    , only : NEXACT
      use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
!
      implicit none
!
      real(8), intent(in)  :: X(3)
      integer, intent(in)  :: Icase,Mdle
!
!  ...exact solution
      real*(8), intent(out) ::   ValH(  MAXEQNH    )
      real*(8), intent(out) ::  DvalH(  MAXEQNH,3  )
      real*(8)              :: d2valH(  MAXEQNH,3,3)
      real*(8), intent(out) ::   ValE(3,MAXEQNE    )
      real*(8), intent(out) ::  DvalE(3,MAXEQNE,3  )
      real*(8)              :: d2valE(3,MAXEQNE,3,3)
      real*(8), intent(out) ::   ValV(3,MAXEQNV    )
      real*(8), intent(out) ::  DvalV(3,MAXEQNV,3  )
      real*(8)              :: d2valV(3,MAXEQNV,3,3)
      real*(8)              ::   valQ(  MAXEQNQ    )
      real*(8)              ::  dvalQ(  MAXEQNQ,3  )
      real*(8)              :: d2valQ(  MAXEQNQ,3,3)
!
      integer :: iprint = 0
!
!--------------------------------------------------------------------
!
!  ...initialize
      ValH = ZERO; DvalH = ZERO
      ValE = ZERO; DvalE = ZERO
      ValV = ZERO; DvalV = ZERO
      select case(NEXACT)
!
!==============================================================================
!  UNKNOWN EXACT SOLUTION                                                     |
!==============================================================================
      case(0)
!
! TODO: Fix this
         if ((X(3).gt.0.99d0).and.(X(2).eq.0.d0)) then
            ValV(2,2) = -1.d0
         endif
!
!==============================================================================
!  KNOWN EXACT SOLUTION                                                       |
!==============================================================================
      case(1,2)
!     ...use the exact solution to determine Dirichlet data
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      case default
         write(*,*)'dirichlet: UNKNOWN EXACT SOLUTION FLAG', NEXACT
         stop 1
      end select
!
#if DEBUG_MODE
      if (iprint == 1) then
         write(*,1001) X(1:3), ValH(1:MAXEQNH)
 1001    format(' dirchlet: X, zvalH = ',3(e12.5,2x),2x,10(e12.5,2x))
      endif
#endif
!
   end subroutine dirichlet
