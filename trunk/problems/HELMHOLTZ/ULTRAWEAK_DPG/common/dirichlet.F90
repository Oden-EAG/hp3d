!---------------------------------------------------------------------
!> @brief      Computes Dirichlet data at a point
!!
!> @param[in]  X        - a point in physical space
!> @param[in]  Icase    - node case (specifies what variables are supported)
!> @param[out] ValH     - value of the H1 solution
!> @param[out] DvalH    - corresponding first derivatives
!> @param[out] DvalE    - value of the H(curl) solution
!> @param[out] DdvalE   - corresponding first derivatives
!> @param[out] DvalV    - value of the H(div) solution
!> @param[out] DdvalV   - corresponding first derivatives
!!
!> @date       July 2023
!----------------------------------------------------------------------
   subroutine dirichlet(Mdle,X,Icase, ValH,DvalH,ValE,DvalE,ValV,DvalV)
!    
      use control,            only: NEXACT
      use parameters,         only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO,ZONE
      use common_prob_data_UW
!   
      implicit none
!   
      real(8), intent(in)  :: X(3)
      integer, intent(in)  :: Icase,Mdle
!
!  ...exact solution
      complex(8) :: ValH  (  MAXEQNH    )
      complex(8) :: DvalH (  MAXEQNH,3  )
      complex(8) :: D2valH(  MAXEQNH,3,3)
      complex(8) :: ValE  (3,MAXEQNE    )
      complex(8) :: DvalE (3,MAXEQNE,3  )
      complex(8) :: D2valE(3,MAXEQNE,3,3)
      complex(8) :: ValV  (3,MAXEQNV    )
      complex(8) :: DvalV (3,MAXEQNV,3  )
      complex(8) :: D2valV(3,MAXEQNV,3,3)
      complex(8) :: ValQ  (  MAXEQNQ    )
      complex(8) :: DvalQ (  MAXEQNQ,3  )
      complex(8) :: D2valQ(  MAXEQNQ,3,3)
!
!  ...printing flag
      integer :: iprint = 0
!   
!--------------------------------------------------------------------
!
!  ...initialize
      ValH = ZERO; DvalH = ZERO
      ValE = ZERO; DvalE = ZERO
      ValV = ZERO; DvalV = ZERO
!   
      select case(NEXACT)
!      
!  ...unknown exact solution
      case(0)
!
         continue
!      
!  ...known exact solution
      case(1,2)
!     ...use the exact solution to determine Dirichlet data
         call exact(X,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                             ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
      case default
         write(*,*)'dirichlet: UNKNOWN EXACT SOLUTION FLAG', NEXACT
         stop 1
      end select
!
   end subroutine dirichlet
