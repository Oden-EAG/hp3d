!---------------------------------------------------------------------------------------
!> @brief     Display user-defined quantity
!!
!> @param[in] Mdle   - element (middle node) number
!> @param[in] Xi     - master element coordinates
!> @param[in] X      - physical coordinates
!> @param[in] Rn     - outward normal unit vector
!> @param[in] ZsolH  - H1    sol
!> @param[in] ZgradH - H1    grad
!> @param[in] ZsolE  - Hcurl sol
!> @param[in] ZcurlE - Hcurl curl
!> @param[in] ZsolV  - Hdiv  sol
!> @param[in] ZdivV  - Hdiv  div
!> @param[in] ZsolQ  - L2    sol
!!
!> @param[out] val   - quantity to display
!!
!> @date      July 2023
!---------------------------------------------------------------------------------------
!
   subroutine soldis(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
!
      use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
      use data_structure3D
!
      implicit none
!
      integer, intent(in)  :: Mdle
      real(8), intent(in)  :: Xi(3), X(3), Rn(3)
      real(8), intent(in)  :: ZsolH(  MAXEQNH), ZgradH(  MAXEQNH,3)
      real(8), intent(in)  :: ZsolE(3,MAXEQNE), ZcurlE(3,MAXEQNE  )
      real(8), intent(in)  :: ZsolV(3,MAXEQNV), ZdivV (  MAXEQNV  )
      real(8), intent(in)  :: ZsolQ(  MAXEQNQ)
      real(8), intent(out) :: Val
!
!  ...exact solution workspace
      integer :: icase
      real(8) ::   zvalH(  MAXEQNH    )
      real(8) ::  zdvalH(  MAXEQNH,3  )
      real(8) :: zd2valH(  MAXEQNH,3,3)
      real(8) ::   zvalE(3,MAXEQNE    )
      real(8) ::  zdvalE(3,MAXEQNE,3  )
      real(8) :: zd2valE(3,MAXEQNE,3,3)
      real(8) ::   zvalV(3,MAXEQNV    )
      real(8) ::  zdvalV(3,MAXEQNV,3  )
      real(8) :: zd2valV(3,MAXEQNV,3,3)
      real(8) ::   zvalQ(  MAXEQNQ    )
      real(8) ::  zdvalQ(  MAXEQNQ,3  )
      real(8) :: zd2valQ(  MAXEQNQ,3,3)
!
      real(8) :: tmpVal(3)
!
!---------------------------------------------------------------------------------------
!
      select case (IEXACT_DISP)
      case(1)
         icase=0
         call exact(X,icase, zvalH,zdvalH,zd2valH, &
                             zvalE,zdvalE,zd2valE, &
                             zvalV,zdvalV,zd2valV, &
                             zvalQ,zdvalQ,zd2valQ)
!     ...exact solution - displacement
         tmpVal = zvalQ(1:3)
      case(2)
!     ...approximate solution - displacement
         tmpVal = ZsolQ(1:3)
      case(3)
         icase=0
         call exact(X,icase, zvalH,zdvalH,zd2valH, &
                             zvalE,zdvalE,zd2valE, &
                             zvalV,zdvalV,zd2valV, &
                             zvalQ,zdvalQ,zd2valQ)

!     ...solution error - displacement
         tmpVal = zvalQ(1:3)-ZsolQ(1:3)
      end select

      select case (ICHOOSE_DISP)
      case(1,2,3)
         Val = tmpVal(ICHOOSE_DISP)
      case(4)
         Val = dsqrt(tmpVal(1)**2 + tmpVal(2)**2 + tmpVal(3)**2)
      end select
!
   end subroutine soldis




!---------------------------------------------------------------------------------------
!> @brief   Show the quantities to display
!---------------------------------------------------------------------------------------
   subroutine soldis_select
!
      use parameters, only : NSTD_OUT
      use common_prob_data, only : ICHOOSE_DISP, IEXACT_DISP
!
      implicit none
!
      integer :: iprev
!
!---------------------------------------------------------------------------------------
!
      write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
      read(*,*) iprev
      if (iprev.eq.1) then
         call disp_soldis(NSTD_OUT); return
      endif
!
      write(*,*) 'EXACT SOLUTION (1) or APPROXIMATE SOLUTION (2) or SIGNED ERROR (3)?'
      read(*,*) IEXACT_DISP
      write(*,200)
200   format('SET VARIABLE     1) u1    2) u2    3) u3    4) ||u||')
!
      read(*,*) ICHOOSE_DISP
      select case(ICHOOSE_DISP)
      case(1,2,3,4)
      case default
         write(*,200)
      end select
!
      call disp_soldis(NSTD_OUT)
!
   end subroutine soldis_select




!---------------------------------------------------------------------------------------
!> @breif   Soldis display interface
!---------------------------------------------------------------------------------------
   subroutine disp_soldis(Nstream)
      use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
!
      implicit none
!
      integer, intent(in) :: Nstream
!
!---------------------------------------------------------------------------------------
!
      write(Nstream,100)
      write(Nstream,310)
!
      select case (IEXACT_DISP)
      case(1); write(Nstream,101)
      case(2); write(Nstream,102)
      case(3); write(Nstream,103)
      end select
!
      select case (ICHOOSE_DISP)
      case(1); write(Nstream,301)
      case(2); write(Nstream,302)
      case(3); write(Nstream,303)
      case(4); write(Nstream,304)
      end select
!
      write(Nstream,310)
!
100   format('DISPLAY SETUP')
101   format('EXACT SOLUTION IS CHOSEN')
102   format('APPROX SOLUTION IS CHOSEN')
103   format('(EXACT SOLUTION - APPROX SOLUTION) IS CHOSEN')
301   format('u1 IS CHOSEN')
302   format('u2 IS CHOSEN')
303   format('u3 IS CHOSEN')
304   format('||u|| IS CHOSEN')
310   format('-----------------------')
!
   end subroutine disp_soldis
