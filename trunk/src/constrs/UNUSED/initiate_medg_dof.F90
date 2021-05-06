#include "typedefs.h"
!-------------------------------------------------------------------------
!> Purpose : initiate H1 and H(curl), i.e., L2 dof for a son of 
!            a mid-edge node
!
!> @date May 20
!
!> @param[in]  Nson  = 1,2,3 (the son number)
!> @param[in]  Nord  - order for a mid-edge node
!> @param[in]  NvarH,NvarE - number of H1 and H(curl) components
!> @param[in]  Xnod  - geometry dof for the edge
!> @param[in]  DofH  - H1 dof for the edge
!> @param[in]  DofE  - H(curl) dof for the edge
!> @param[in]  NdofH,NdofE  - number of dof for the son
!
!> @param[out] Xnod_son    geometry dof for the sons
!> @param[out] DofH_son    - H1 dof for the sons
!> @param[out] DofE_son             - H(curl) dof for the sons
!-------------------------------------------------------------------------
!
      subroutine initiate_medg_dof(Nson,Nord,NvarH,NvarE, &
                                   Xnod,DofH,DofE, &
                                   Xnod_son,DofH_son,DofE_son, &
                                   NdofH,NdofE)
!
      use constraints
      implicit none
!
      integer,             intent(in)  :: Nson,Nord,NvarH,NvarE,NdofH,NdofE
      real*8, dimension(NDIMEN,MAXP+1),  intent(in)  :: Xnod
      VTYPE,  dimension(NvarH, MAXP+1),  intent(in)  :: DofH
      VTYPE,  dimension(NvarE, MAXP),    intent(in)  :: DofE
!
      real*8, dimension(NDIMEN,NdofH),   intent(out) :: Xnod_son
      VTYPE,  dimension(NvarH ,NdofH),   intent(out) :: DofH_son
      VTYPE,  dimension(NvarE, NdofE),   intent(out) :: DofE_son
!
!  ...locals
      integer ::  is
      real*8, dimension(MAXP+1,MAXP-1) :: coeff
!
      select case(Nson)
!  
!  ...mid-edge sons
      case(1,2)
        coeff(1:Nord-1,1:Nord-1) = RRRH(1,1:Nord-1,Nson,1:Nord-1)
        call Rtransform1(Xnod(1,3),NDIMEN, Nord-1,coeff, Xnod_son,Nord-1)
        call Vtransform1(DofH(1,3),NvarH,  Nord-1,coeff, DofH_son,Nord-1)
!
        coeff(1:Nord,1:Nord) = RRRE(1,1:Nord,is,1:Nord)
        call Vtransform1(DofE,    NvarE,   Nord,  coeff, DofE_son,Nord)
!
!  ...vertex son
      case(3)
        coeff(1:2,1) = .5d0
        coeff(3:Nord+1,1) = RRRH(1,1:Nord-1,3,1)
        call Rtransform1(Xnod,NDIMEN, Nord+1,coeff, Xnod_son,1)
        call Vtransform1(DofH,NvarH,  Nord+1,coeff, DofH_son,1)
      end select
!
!
      end subroutine initiate_medg_dof


      subroutine Rtransform1(Uu,M,Iu,Coeff, Vv,Iv)
      implicit none
!
      integer,                  intent(in)  :: M,Iu,Iv
      real*8, dimension(M,Iu),  intent(in)  :: Uu
      real*8, dimension(Iu,Iv), intent(in)  :: Coeff
      real*8, dimension(M,Iv),  intent(out) :: Vv
!
      integer :: i,j
!
      do i=1,Iu
        do j=1,Iv
          Vv(1:M,j) = Vv(1:M,j) &
                    + Uu(1:M,i)*Coeff(i,j)
        enddo
      enddo
!
      end subroutine Rtransform1

      subroutine Vtransform1(Uu,M,Iu,Coeff, Vv,Iv)
      use parameters, only: ZERO
      implicit none
!
      integer,                 intent(in)  :: M,Iu,Iv
      VTYPE, dimension(M,Iu),  intent(in)  :: Uu
      real*8,dimension(Iu,Iv), intent(in)  :: Coeff
      VTYPE, dimension(M,Iv),  intent(out) :: Vv
!
      integer :: i,j
!
      do i=1,Iu
        do j=1,Iv
          Vv(1:M,j) = Vv(1:M,j) &
                    + Uu(1:M,i)*Coeff(i,j)
        enddo
      enddo
!
      end subroutine Vtransform1

