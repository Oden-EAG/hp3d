#include "typedefs.h"
!-------------------------------------------------------------------------
!> Purpose : initiate H1, H(curl) and H(div), i.e. L2 dof for a son of 
!            a mid-face rectangular node
!
!> @date May 20
!
!> @param[in]  Kref    - refinement type
!> @param[in]  Nson    = 1,...,9 (the son number)
!> @param[in]  Norder  - polynomial order for the quad element nodes
!> @param[in]  Norient - orientation for the quad element nodes
!> @param[in]  NvarH,NvarE,NvarV - number of H1, H(curl) and H(div) components
!> @param[in]  Xnod    - geometry dof for the face
!> @param[in]  DofH    - H1 dof for the face
!> @param[in]  DofE    - H(curl) dof for the face
!> @param[in]  DofV    - H(div) (L2) dof for the face
!> @param[in]  NdofH,NdofE,NdofV - number of dof for the son
!
!> @param[out] Xnod_son - geometry dof for the son
!> @param[out] DofH_son - H1 dof for the son
!> @param[out] DofE_son - H(curl) dof for the son
!> @param[out] DofV_son - H(div) dof for the son
!-------------------------------------------------------------------------
!
      subroutine initiate_mdlq_dof(Kref,Nson,Norder,Norient,NvarH,NvarE,NvarV, &
                                   Xnod,DofH,DofE,DofV, &
                                   Xnod_son,DofH_son,DofE_son,DofV_son, &
                                   NdofH,NdofE,NdofV)
!
      use constraints
      implicit none
!
      integer, intent(in)  :: Kref,Nson,NvarH,NvarE,NvarV,NdofH,NdofE,NdofV
      integer,dimension(5),  intent(in)                :: Norder
      integer,dimension(4),  intent(in)                :: Norient
      real*8, dimension(NDIMEN,MAXquadH),  intent(in)  :: Xnod
      VTYPE,  dimension(NvarH, MAXquadH),  intent(in)  :: DofH
      VTYPE,  dimension(NvarE, MAXquadE),  intent(in)  :: DofE
      VTYPE,  dimension(NvarV, MAXquadV),  intent(in)  :: DofV
!
      real*8, dimension(NDIMEN,NdofH),  intent(out) :: Xnod_son
      VTYPE,  dimension(NvarH, NdofH),  intent(out) :: DofH_son
      VTYPE,  dimension(NvarE, NdofE),  intent(out) :: DofE_son
      VTYPE,  dimension(NvarV ,NdofV),  intent(out) :: DofV_son
!
!  ...locals
      integer                          :: nord1,nord2,i,j,naH,naE
      real*8, dimension(MAXP+1,MAXP-1) :: coeffH1,coeffH2
      real*8, dimension(MAXP,MAXP)     :: coeffE1,coeffE2
      real*8, dimension(MAXP+1)        :: rHv
!
      real*8, allocatable  :: xnod1(:,:,:)
      VTYPE,  allocatable  :: dofH1(:,:,:)
      VTYPE,  allocatable  :: dofE1(:,:,:)
      VTYPE,  allocatable  :: dofE2(:,:,:)
!
      rHv(1:2) = .5d0; rHv(3:MAXP+1) = RRRH(1,1:MAXP-1,3,1)
      naH=0; naE=0
      do i=1,4
        naH = naH + Norder(i)-1; naE = naE + Norder(i)
      enddo
      call decode(Norder(5), nord1,nord2)
      allocate(xnod1(NDIMEN,nord1+1,nord2+1), dofH1(NvarH,nord1+1,nord2+1), &
               dofE1(NvarE,nord1,nord2+1), dofE2(NvarE,nord1+1,nord2))
      call reset_mdlq_dof(Norder(1:4),Norient,NvarH,NvarE,nord1,nord2, &
                          Xnod,DofH,DofE, xnod1,dofH1,dofE1,dofE2)
!
      Xnod_son = 0.d0; DofH_son = ZERO; DofE_son = ZERO; DofV_son = ZERO
!
      select case(Kref)
!
!  ...isotropic refinement
      case(11)
!  
      select case(Nson)
!
!  ...mid-face sons
      case(1,2,3,4)
        select case(Nson)
        case(1); i=1;j=1
        case(2); i=2;j=1
        case(3); i=2;j=2
        case(4); i=1;j=2
        end select
        coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,i,1:nord1-1)
        coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,j,1:nord2-1)
        call Rtransform2(Xnod(1,naH+1),NDIMEN,nord1-1,nord2-1, &
                         coeffH1,coeffH2,  &
                         Xnod_son,nord1-1,nord2-1)
        call Vtransform2(DofH(1,naH+1),NvarH,nord1-1,nord2-1, &
                         coeffH1,coeffH2,  &
                         DofH_son,nord1-1,nord2-1)
!
        coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,i,1:nord1)
        coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,j,1:nord2)
        call Vtransform2(DofE(1,naE+1),NvarE,nord1,nord2-1, &
                         coeffE1,coeffH2,  &
                         DofE_son,nord1,nord2-1)
        call Vtransform2(DofE(1,naE+nord1*(nord2-1)+1),NvarE,nord1-1,nord2, &
                         coeffH1,coeffE2,  &
                         DofE_son(1,nord1*(nord2-1)+1),nord1-1,nord2)
!
        call Vtransform2(DofV,NvarV,nord1,nord2, &
                         coeffE1,coeffE2,  &
                         DofV_son,nord1,nord2)
!
!  ...fourth or second mid-edge son
      case(8,6)
        select case(Nson)
        case(8); i=1
        case(6); i=2
        end select
        coeffH1(1:2,1:nord1-1) = 0.d0
        coeffH1(3:nord1+1,1:nord1-1) = RRRH(1,1:nord1-1,i,1:nord1-1)
        call Rtransform2(Xnod1,NDIMEN,nord1+1,nord2+1, &
                         coeffH1,rHv,  &
                         Xnod_son,nord1-1,1)
!
        call Vtransform2(DofH1,NvarH,nord1+1,nord2+1, &
                         coeffH1,rHv,  &
                         DofH_son,nord1-1,1)
!
        coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,i,1:nord1)
        call Vtransform2(DofE1,NvarE,nord1,nord2+1, &
                         coeffE1,rHv,  &
                         DofE_son,nord1,1)
!
!  ...first or third mid-edge son
      case(5,7)
        select case(Nson)
        case(5); j=1
        case(7); j=2
        end select
        coeffH2(1:2,1:nord2-1) = 0.d0
        coeffH2(3:nord2+1,1:nord2-1) = RRRH(1,1:nord2-1,j,1:nord2-1)
        call Rtransform2(Xnod1,NDIMEN,nord1+1,nord2+1, &
                         rHv,coeffH2,  &
                         Xnod_son,1,nord2-1)
!
        call Vtransform2(DofH1,NvarH,nord1+1,nord2+1, &
                         rHv,coeffH2,  &
                         DofH_son,1,nord2-1)
!
        coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,j,1:nord2)
        call Vtransform2(DofE2,NvarE,nord1+1,nord2, &
                         rHv,coeffE2,  &
                         DofE_son,1,nord2)
!
!  ...vertex son
      case(9)
        call Rtransform2(Xnod1,NDIMEN,nord1+1,nord2+1, &
                         rHv,rHv,  &
                         Xnod_son,1,1)
!
        call Vtransform2(DofH1,NvarH,nord1+1,nord2+1, &
                         rHv,rHv,  &
                         DofH_son,1,1)
      end select
!
!  ...anisotropic refinement across the first axis
      case(10)
        select case(Nson)
!
!  .....mid-face sons
        case(1,2)
          coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,Nson,1:nord1-1)
          call setI(coeffH2,nord2-1)
          call Rtransform2(Xnod(1,naH+1),NDIMEN,nord1-1,nord2-1, &
                           coeffH1,coeffH2,  &
                           Xnod_son,nord1-1,nord2-1)
          call Vtransform2(DofH(1,naH+1),NvarH,nord1-1,nord2-1, &
                           coeffH1,coeffH2,  &
                           DofH_son,nord1-1,nord2-1)
!
          coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,Nson,1:nord1)
          call setI(coeffE2, nord2)
          call Vtransform2(DofE(1,naE+1),NvarE,nord1,nord2-1, &
                           coeffE1,coeffH2,  &
                           DofE_son,nord1,nord2-1)
          call Vtransform2(DofE(1,naE+nord1*(nord2-1)),NvarE,nord1-1,nord2, &
                           coeffH1,coeffE2,  &
                           DofE_son(1,nord1*(nord2-1)+1),nord1-1,nord2)
!
          call Vtransform2(DofV,NvarV,nord1,nord2, &
                           coeffE1,coeffE2,  &
                           DofV_son,nord1,nord2)
!
!  .....mid-edge son
        case(3)
          call setI0(coeffH2,nord2-1)
          call Rtransform2(Xnod1,NDIMEN,nord1+1,nord2+1, &
                           rHv,coeffH2,  &
                           Xnod_son,1,nord2-1)
!
          call Vtransform2(DofH1,NvarH,nord1+1,nord2+1, &
                           rHv,coeffH2,  &
                           DofH_son,1,nord2-1)
!
          call setI(coeffE2,nord2)
          call transform2(DofE1,NvarE,nord1+1,nord2, &
                          rHv,coeffE2,  &
                          DofE_son,1,nord2)
        end select
!
!  ...anisotropic refinement across the second axis
      case(01)
        select case(Nson)
!
!  .....mid-face sons
        case(1,2)
          call setI(coeffH1, nord1-1)
          coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,Nson,1:nord2-1)
          call Rtransform2(Xnod(1,naH+1),NDIMEN,nord1-1,nord2-1, &
                           coeffH1,coeffH2,  &
                           Xnod_son,nord1-1,nord2-1)
          call Vtransform2(DofH(1,naH+1),NvarH,nord1-1,nord2-1, &
                           coeffH1,coeffH2,  &
                           DofH_son,nord1-1,nord2-1)
!
          call setI(coeffE1,nord1)
          coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,Nson,1:nord2)
          call Vtransform2(DofE(1,naE+1),NvarE,nord1,nord2-1, &
                           coeffE1,coeffH2,  &
                           DofE_son,nord1,nord2-1)
          call Vtransform2(DofE(1,naE+nord1*(nord2-1)),NvarE,nord1-1,nord2, &
                           coeffH1,coeffE2,  &
                           DofE_son(1,nord1*(nord2-1)+1),nord1-1,nord2)
!
          call Vtransform2(DofV,NvarV,nord1,nord2, &
                           coeffE1,coeffE2,  &
                           DofV_son,nord1,nord2)
!
!  .....mid-edge son
        case(3)
          call setI0(coeffH1, nord1-1)
          call Rtransform2(Xnod1,NDIMEN,nord1+1,nord2+1, &
                           coeffH1,rhv,  &
                           Xnod_son,nord1-1,1)
          call Vtransform2(DofH1,nvarH,nord1+1,nord2+1, &
                           coeffH1,rhv,  &
                           DofH_son,nord1-1,1)
!
          call setI(coeffE1,nord1)
          call Vtransform2(DofE1,NvarE,nord1,nord2+1, &
                           coeffE1,rHv,  &
                           DofE_son,nord1,1)
        end select
!  
      end select
      deallocate(xnod1,dofH1,dofE1,dofE2)
!
!
      end subroutine initiate_mdlq_dof


      subroutine Rtransform2(Uu,M,Iu,Ju,Coeff1,Coeff2, Vv,Iv,Jv)
      use parameters, only: ZERO
      implicit none
!
      integer,                   intent(in)  :: M,Iu,Ju,Iv,Jv
      real*8, dimension(M,Iu,Ju), intent(in)  :: Uu
      real*8, dimension(Iu,Iv),   intent(in)  :: Coeff1
      real*8, dimension(Ju,Jv),   intent(in)  :: Coeff2
      real*8, dimension(M,Iv,Jv), intent(out) :: Vv
!
      integer :: i,j
      real*8, dimension(M,Ju,Iv)  :: aa
!
      aa = ZERO
      do i=1,Iu
        do j=1,Iv
          aa(1:M,1:Ju,j) = aa(1:M,1:Ju,j) &
                         + Uu(1:M,j,1:Ju)*Coeff1(i,j)
        enddo
      enddo
      do i=1,Ju
        do j=1,Jv
          Vv(1:M,1:Iv,j) = Vv(1:M,1:Iv,j) &
                         + aa(1:M,i,1:Iv)*Coeff2(i,j)
        enddo
      enddo
!
      end subroutine Rtransform2

      subroutine Vtransform2(Uu,M,Iu,Ju,Coeff1,Coeff2, Vv,Iv,Jv)
      use parameters, only: ZERO
      implicit none
!
      integer,                   intent(in)  :: M,Iu,Ju,Iv,Jv
      VTYPE, dimension(M,Iu,Ju), intent(in)  :: Uu
      real*8,dimension(Iu,Iv),   intent(in)  :: Coeff1
      real*8,dimension(Ju,Jv),   intent(in)  :: Coeff2
      VTYPE, dimension(M,Iv,Jv), intent(out) :: Vv
!
      integer :: i,j
      VTYPE, dimension(M,Ju,Iv)  :: aa
!
      aa = ZERO
      do i=1,Iu
        do j=1,Iv
          aa(1:M,1:Ju,j) = aa(1:M,1:Ju,j) &
                         + Uu(1:M,j,1:Ju)*Coeff1(i,j)
        enddo
      enddo
      do i=1,Ju
        do j=1,Jv
          Vv(1:M,1:Iv,j) = Vv(1:M,1:Iv,j) &
                         + aa(1:M,i,1:Iv)*Coeff2(i,j)
        enddo
      enddo
!
      end subroutine Vtransform2








      subroutine setI(A,Na)
      use parameters, only: ZERO
      implicit none
!
      integer,                 intent(in)  :: Na
      VTYPE, dimension(Na,Na), intent(out) :: A
!
      integer :: i
!
      A = ZERO
      do i=1,Na
        A(i,i) = 1.d0
      enddo
!
      end subroutine setI









      subroutine setI0(A,Na)
      use parameters, only: ZERO
      implicit none
!
      integer,                 intent(in)  :: Na
      VTYPE, dimension(Na+2,Na), intent(out) :: A
!
      integer :: i
!
      A = ZERO
      do i=1,Na
        A(2+i,i) = 1.d0
      enddo
!
      end subroutine setI0

