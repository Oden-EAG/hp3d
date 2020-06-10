!-------------------------------------------------------------------------------------
!> Purpose - Selection for quantities to be displayed by graphics
!!
!> @data Nov 14
!-------------------------------------------------------------------------------------
!
subroutine soldis_select
!
      implicit none
!
!     redirect to system routine
      call soldis_select_system
!      
!      
endsubroutine soldis_select
!
!
!
!-------------------------------------------------------------------------------------
!> Purpose - Compute quantities to be displayed by graphics
!!
!> @param[in]  Mdle   - middle node
!> @param[in]  Xi     - master element coordinates
!> @param[in]  X      - coordinates of a point
!> @param[in]  Rn     - outward normal unit vector
!> @param[in]  ZsolH  - value of H1      solution 
!> @param[in]  ZgradH - grad  of H1      solution
!> @param[in]  ZsolE  - value of H(curl) solution
!> @param[in]  ZcurlE - curl  of H(curl) solution
!> @param[in]  ZsolV  - value of H(div)  solution
!> @param[in]  ZdivV  - div   of H(div)  solution
!> @param[in]  ZsolQ  - value of L^2     solution
!> @param[out] Val    - quantity to display
!!
!> @data Nov 14
!-------------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine soldis(Mdle,Xi,X,Rn,SolH,GradH,SolE,CurlE,SolV,DivV,SolQ, Val)
!
      use data_structure3D , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ
!
      implicit none
      integer,                     intent(in)  :: Mdle
      real*8,dimension(3),         intent(in)  :: Xi,X,Rn
      VTYPE,dimension(  MAXEQNH  ),intent(in)  :: SolH
      VTYPE,dimension(  MAXEQNH,3),intent(in)  :: GradH
      VTYPE,dimension(3,MAXEQNE  ),intent(in)  :: SolE
      VTYPE,dimension(3,MAXEQNE  ),intent(in)  :: CurlE
      VTYPE,dimension(3,MAXEQNV  ),intent(in)  :: SolV
      VTYPE,dimension(  MAXEQNV  ),intent(in)  :: DivV
      VTYPE,dimension(  MAXEQNQ  ),intent(in)  :: SolQ
      real*8,                      intent(out) :: Val
!
!-------------------------------------------------------------------------------------
!
!     redirect to system routine
      call soldis_system(Mdle,Xi,X,Rn,SolH,GradH,SolE,CurlE,SolV,DivV,SolQ, Val)
!
!
endsubroutine soldis
