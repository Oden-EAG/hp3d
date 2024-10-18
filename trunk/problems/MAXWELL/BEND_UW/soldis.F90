!
#include "typedefs.h"
!
!---------------------------------------------------------------------------------------
!> @brief     Display user-defined quantity
!!
!> @param[in]  Mdle   - element (middle node) number
!> @param[in]  Xi     - master element coordinates
!> @param[in]  X      - physical coordinates
!> @param[in]  Rn     - outward normal unit vector
!> @param[in]  ZsolH  - H1    sol
!> @param[in]  ZgradH - H1    grad
!> @param[in]  ZsolE  - Hcurl sol
!> @param[in]  ZcurlE - Hcurl curl
!> @param[in]  ZsolV  - Hdiv  sol
!> @param[in]  ZdivV  - Hdiv  div
!> @param[in]  ZsolQ  - L2    sol
!> @param[out] val   - quantity to display
!!
!> @date       July 2023
!---------------------------------------------------------------------------------------
   subroutine soldis(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
!
      use data_structure3D
      use commonParam
!
      implicit none
!
      integer,    intent(in)  :: Mdle
      real(8),    intent(in)  :: Xi(3), X(3), Rn(3)
      complex(8), intent(in)  :: ZsolH (  MAXEQNH  )
      complex(8), intent(in)  :: ZgradH(  MAXEQNH,3)
      complex(8), intent(in)  :: ZsolE (3,MAXEQNE  )
      complex(8), intent(in)  :: ZcurlE(3,MAXEQNE  )
      complex(8), intent(in)  :: ZsolV (3,MAXEQNV  )
      complex(8), intent(in)  :: ZdivV (  MAXEQNV  )
      complex(8), intent(in)  :: ZsolQ (  MAXEQNQ  )
      real(8),    intent(out) :: Val
!
!  ...work space for routine 'exact'
      complex(8) :: rntimesE(3)
      complex(8) :: zvalH(MAXEQNH)  , zdvalH(MAXEQNH,3)  , zd2valH(MAXEQNH,3,3)
      complex(8) :: zvalE(3,MAXEQNE), zdvalE(3,MAXEQNE,3), zd2valE(3,MAXEQNE,3,3)
      complex(8) :: zvalV(3,MAXEQNV), zdvalV(3,MAXEQNV,3), zd2valV(3,MAXEQNV,3,3)
      complex(8) :: zvalQ(MAXEQNQ)  , zdvalQ(MAXEQNQ,3)  , zd2valQ(MAXEQNQ,3,3)
!
      integer :: icase = 1, iprint
!
!!===========================
!!  ...redirect to system routine
!      call soldis_system(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
!      return
!!===========================
!
      iprint=0
      if (iprint.eq.1) then
         write(*,7001) Mdle,X(1:3)
 7001    format('soldis: Mdle = ',i5,'  X = ',3f8.3)
         write(*,7002) ZsolH(1)
 7002    format('ZsolH = ',e12.5)
      endif
!
      select case (IEXACT_DISP)
!  ...exact solution
      case(1)
         call exact(X, icase, &
                    zvalH,zdvalH,zd2valH, &
                    zvalE,zdvalE,zd2valE, &
                    zvalV,zdvalV,zd2valV, &
                    zvalQ,zdvalQ,zd2valQ)
         select case (ICHOOSE_COMP)
!
!     ...exact E (tangential) trace
         case(1)
            call zcross_product(rn,ZvalE(1:3,1), rntimesE)
            Val = dreal(rntimesE(1))
         case(2)
            call zcross_product(rn,ZvalE(1:3,1), rntimesE)
            Val = dreal(rntimesE(2))
         case(3)
            call zcross_product(rn,ZvalE(1:3,1), rntimesE)
            Val = dreal(rntimesE(3))
!
!     ...exact H (tangential) trace
         case(4)
            call zcross_product(rn,ZvalE(1:3,2), rntimesE)
            Val = dreal(rntimesE(1))
         case(5)
            call zcross_product(rn,ZvalE(1:3,2), rntimesE)
            Val = dreal(rntimesE(2))
         case(6)
            call zcross_product(rn,ZvalE(1:3,2), rntimesE)
            Val = dreal(rntimesE(3))
!
!     ...exact E field
         case(7)
            Val = dreal(zvalQ(1))
         case(8)
            Val = dreal(zvalQ(2))
         case(9)
            Val = dreal(zvalQ(3))
!
!     ...exact H field
         case(10)
            Val = dreal(zvalQ(4))
         case(11)
            Val = dreal(zvalQ(5))
         case(12)
            Val = dreal(zvalQ(6))
         end select
!
!  ...approximate solution
      case(0)
!
         select case (ICHOOSE_COMP)
!
!     ...exact E (tangential) trace
         case(1)
            call zcross_product(rn,zsolE(1:3,1), rntimesE)
            Val = dreal(rntimesE(1))
         case(2)
            call zcross_product(rn,zsolE(1:3,1), rntimesE)
            Val = dreal(rntimesE(2))
         case(3)
            call zcross_product(rn,zsolE(1:3,1), rntimesE)
            Val = dreal(rntimesE(3))
!
!     ...exact H (tangential) trace
         case(4)
            call zcross_product(rn,zsolE(1:3,2), rntimesE)
            Val = dreal(rntimesE(1))
         case(5)
            call zcross_product(rn,zsolE(1:3,2), rntimesE)
            Val = dreal(rntimesE(2))
         case(6)
            call zcross_product(rn,zsolE(1:3,2), rntimesE)
            Val = dreal(rntimesE(3))
!
!     ...exact E field
         case(7)
            Val = dreal(zsolQ(1))
         case(8)
            Val = dreal(zsolQ(2))
         case(9)
            Val = dreal(zsolQ(3))
!
!     ...exact H field
         case(10)
            Val = dreal(zsolQ(4))
         case(11)
            Val = dreal(zsolQ(5))
         case(12)
            Val = dreal(zsolQ(6))
!
         end select
      end select
      if (iprint.eq.1) then
         write(*,7003) IEXACT_DISP, ICHOOSE_COMP, Val
 7003    format('IEXACT_DISP, ICHOOSE_COMP, Val = ',2i2,e12.5)
         call pause
      endif
!
   end subroutine soldis
!
!
!---------------------------------------------------------------------------------------
!> @brief select quantities to display
!---------------------------------------------------------------------------------------
   subroutine soldis_select
!
      use control,    only: NEXACT    ! exact solution flag
      use parameters, only: NSTD_OUT  ! display file
      use commonParam
!
      implicit none
      integer :: iprev
!---------------------------------------------------------------------------------------
!
!!===========================
!!  ...redirect to system routine if you do not want to customize the routine
!      call soldis_select_system
!      return
!!===========================
!
      select case(ICHOOSE_COMP)
      case(1,2)
        write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
        read(*,*) iprev
        if (iprev.eq.1) then
          call disp_soldis(NSTD_OUT) ; return
        endif
      end select
!
      IEXACT_DISP=0
      if ((NEXACT.eq.1).or.(NEXACT.eq.2)) then
   10   write(*,*) 'DISPLAY APPROXIMATE OR EXACT SOLUTION (0/1) ?'
        read(*,*) IEXACT_DISP
        if ((IEXACT_DISP.ne.0).and.(IEXACT_DISP.ne.1)) goto 10
      endif
!
   20 write(*,*) 'SET VARIABLE: EEhat(1-3), HHhat(4-6), EHfld(7-12)'
      read(*,*) ICHOOSE_COMP
      if ((ICHOOSE_COMP.lt.1).and.(ICHOOSE_COMP.gt.12)) goto 20
!
      call disp_soldis(NSTD_OUT)
!
   end subroutine soldis_select
!
!
!---------------------------------------------------------------------------------------
!> @brief show the quantities selected to display
!---------------------------------------------------------------------------------------
   subroutine disp_soldis(Nstream)
!
      use commonParam, only: IEXACT_DISP
!
      implicit none
      integer, intent(in) :: Nstream
!
      write(Nstream,1000)
      write(Nstream,3100)
      select case (IEXACT_DISP)
      case(1); write(Nstream,1010)
      case(0); write(Nstream,1020)
      end select
!
      write(Nstream,3100)
!
 1000 format('DISPLAY SETUP')
 1010 format('DISPLAYING EXACT SOLUTION')
 1020 format('DISPLAYING APPROXIMATE SOLUTION')
 3100 format('-----------------------')
!
   end subroutine disp_soldis
