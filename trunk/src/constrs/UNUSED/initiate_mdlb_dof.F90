#include "typedefs.h"
!-------------------------------------------------------------------------
!> Purpose : initiate H1, H(curl) and H(div), and L2 dof for a son of 
!            the middle node of a brick element
!
!> @date May 20
!
!> @param[in]  Kref    - refinement type
!> @param[in]  Nson    - the son number
!> @param[in]  Norder  - polynomial order for the element
!> @param[in]  NvarH,NvarE,NvarV,NvarQ - number of H1, H(curl),H(div) and
!                        L2 components
!> @param[in]  Xnod    - geometry dof for the element
!> @param[in]  DofH    - H1 dof for the element
!> @param[in]  DofE    - H(curl) dof for the element
!> @param[in]  DofV    - H(div) dof for the element
!> @param[in]  DofQ    - L2  dof for the element
!
!> @param[out] Xnod_son - geometry dof for the son
!> @param[out] DofH_son - H1 dof for the son
!> @param[out] DofE_son - H(curl) dof for the son
!> @param[out] DofV_son - H(div) dof for the son
!> @param[out] DofQ_son - L2 dof for the son
!-------------------------------------------------------------------------
!
      subroutine initiate_mdlb_dof(Kref,Nson,Norder,NvarH,NvarE,NvarV,NvarQ &
                                   Xnod,DofH,DofE,DofV,DofQ, &
                                   Xnod_son,DofH_son,DofE_son,DofV_son,DofQ_son)
!
      use constraints
      implicit none
!
      integer,              intent(in)  :: Kref,Nson,NvarH,NvarE,NvarV,NvarQ
      integer,dimension(19),               intent(in)  :: Norder
      real*8, dimension(NDIMEN,MAXquadH),  intent(in)  :: Xnod
      VTYPE,  dimension(NvarH ,MAXquadH),  intent(in)  :: DofH
      VTYPE,  dimension(NvarE, MAXquadE),  intent(in)  :: DofE
      VTYPE,  dimension(NvarV, MAXquadV),  intent(in)  :: DofV
      VTYPE,  dimension(NvarQ, MAXquadQ),  intent(in)  :: DofQ
!
      real*8, dimension(NDIMEN,MAXmdlbH),  intent(out) :: Xnod_son
      VTYPE,  dimension(NvarH, MAXmdlbH),  intent(out) :: DofH_son
      VTYPE,  dimension(NvarE, MAXmdlbE),  intent(out) :: DofE_son
      VTYPE,  dimension(NvarV ,MAXmdlbV),  intent(out) :: DofV_son
      VTYPE,  dimension(NvarQ ,MAXmdlbQ),  intent(out) :: DofQ_son
!
!  ...locals
      character(len=4)                 :: nod_type
      integer                          :: nord1,nord2,nord3,i,j,k,ndofH,ndofE
      integer,dimension(19)            :: naH,naE,naV
      real*8, dimension(MAXP+1,MAXP-1) :: coeffH1,coeffH2,coeffH3
      real*8, dimension(MAXP,MAXP)     :: coeffE1,coeffE2,coeffE3
!
      Xnod_son = 0.d0
      DofH_son = ZERO; DofE_son = ZERO; DofV_son = ZERO; DofQ_son = ZERO
!
      ndofH=8; ndofE=0; ndofV=0
      do i=1,12
        naH(i) = ndofH; naE(i) = ndofE; naV(i) = ndofV
        ndofH = ndofH + Norder(i)-1; ndofE = ndofE + Norder(i)
      enddo
      do i=13,18
        call decode(Norder(i), nord1,nord2)
        naH(i) = ndofH; naE(i) = ndofE; naV(i) = ndofV
        ndofH = ndofH + (nord1-1)*(nord2-1)
        ndofE = ndofE + nord1*(nord2-1) + (nord1-1)*nord2
        ndofV = ndofV + nord1*nord2
      enddo
      naH(19) = ndofH; naE(19) = ndofE; naV(19) = ndofV
!
      call ddecode(Norder(19), nord1,nord2,nord3)
!
      select case(Kref)
!
!  ...isotropic refinement
      case(111)
!  
        select case(Nson)
!
!  .....middle nodes sons
        case(1,2,3,4,5,6,7,8)
          select case(Nson)
          case(1); i=1;j=1;k=1
          case(2); i=2;j=1;k=1
          case(3); i=2;j=2;k=1
          case(4); i=1;j=2;k=1
          case(5); i=1;j=1;k=2
          case(6); i=2;j=1;k=2
          case(7); i=2;j=2;k=2
          case(8); i=1;j=2;k=2
          end select
          coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,i,1:nord1-1)
          coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,j,1:nord2-1)
          coeffH3(1:nord3-1,1:nord3-1) = RRRH(1,1:nord3-1,j,1:nord3-1)
          coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,i,1:nord1)
          coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,j,1:nord2)
          coeffE3(1:nord3,1:nord3) = RRRE(1,1:nord3,j,1:nord3)
          node = 'mdlb'
!
!  .....mid-face sons parallel to axes 1,2
        case(13,14,15,16)
          select case(Nson)
          case(13); i=1;j=1
          case(14); i=2;j=1
          case(15); i=2;j=2
          case(16); i=1;j=2
          end select
          coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,i,1:nord1-1)
          coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,j,1:nord2-1)
          coeffH3(1:nord3-1,1) = RRRH(1,1:nord3-1,3,1)
          coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,i,1:nord1)
          coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,j,1:nord2)
          nod_type = 'mf12'
!
!  .....mid-face sons parallel to axes 2,3
        case(9,11,17,19)
          select case(Nson)
          case(9);  j=1;k=1
          case(11); j=2;k=1
          case(17); j=2;k=2
          case(19); j=1;k=2
          end select
          coeffH1(1:nord1-1,1) = RRRH(1,1:nord1-1,3,1)
          coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,j,1:nord2-1)
          coeffH3(1:nord3-1,1:nord3-1) = RRRH(1,1:nord3-1,k,1:nord3-1)
          coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,j,1:nord2)
          coeffE3(1:nord3,1:nord3) = RRRE(1,1:nord3,k,1:nord3)
          nod_type = 'mf23'
!
!  .....mid-face sons parallel to axes 1,3
        case(12,10,20,18)
          select case(Nson)
          case(12); i=1;k=1
          case(10); i=2;k=1
          case(20); i=2;k=2
          case(18); i=1;k=2
          end select
          coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,i,1:nord1-1)
          coeffH2(1:nord2-1,1) = RRRH(1,1:nord2-1,3,1)
          coeffH3(1:nord3-1,1:nord3-1) = RRRH(1,1:nord3-1,k,1:nord3-1)
          coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,j,1:nord1)
          coeffE3(1:nord3,1:nord3) = RRRE(1,1:nord3,i,1:nord3)
          nod_type = 'mf13'








        coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,i,1:nord1-1)
        coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,j,1:nord2-1)
        call transform2(Xnod(1,naH(5)+1),NDIMEN,nord1-1,nord2-1, &
                        coeffH1,coeffH2,  &
                        Xnod_son,nord1-1,nord2-1)
        call transform2(DofH(1,naH(5)+1),NvarH,nord1-1,nord2-1, &
                        coeffH1,coeffH2,  &
                        DofH_son,nord1-1,nord2-1)
        coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,i,1:nord1)
        coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,j,1:nord2)
        call transform2(DofE(1,naE(5)+1),NvarE,nord1,nord2-1, &
                        coeffE1,coeffH2,  &
                        DofE_son,nord1,nord2-1)
        call transform2(DofE(1,naE(5)+nord1*(nord2-1)),NvarE,nord1-1,nord2, &
                        coeffH1,coeffE2,  &
                        DofE_son(1,nord1*(nord2-1)+1),nord1-1,nord2)
        call transform2(DofV,NvarV,nord1,nord2, &
                        coeffE1,coeffE2,  &
                        DofV_son,nord1,nord2)
!
!  ...fourth or second mid-edge son
      case(8,6)
        select case(Nson)
        case(8); i=1
        case(6); i=2
        end select
        coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,i,1:nord1-1)
        coeffH2(1:nord2-1,1) = RRRH(1,1:nord2-1,3,1)
        call transform2(Xnod(1,naH(5)+1),NDIMEN,nord1-1,nord2-1, &
                        coeffH1,coeffH2,  &
                        Xnod_son,nord1-1,1)
        call transform2(Xnod(1,naH(1)+1),NDIMEN,Norder(1)-1,1, &
                        coeffH1,.5d0,  &
                        Xnod_son,Norder(1)-1,1)
        call transform2(Xnod(1,naH(3)+1),NDIMEN,Norder(3)-1,1, &
                        coeffH1,.5d0,  &
                        Xnod_son,Norder(3)-1,1)
!
        call transform2(DofH(1,naH(5)+1),NvarH,nord1-1,nord2-1, &
                        coeffH1,coeffH2,  &
                        DofH_son,nord1-1,1)
        call transform2(DofH(1,naH(1)+1),NvarH,Norder(1)-1,1, &
                        coeffH1,.5d0,  &
                        DofH_son,Norder(1)-1,1)
        call transform2(DofH(1,naH(3)+1),NvarH,Norder(3)-1,1, &
                        coeffH1,.5d0,  &
                        DofH_son,Norder(3)-1,1)
!
        coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,i,1:nord1)
        call transform2(DofE(1,naE(5)+1),NvarE,nord1,nord2-1, &
                        coeffE1,coeffH2,  &
                        DofE_son,nord1,1)
        call transform2(DofE(1,naE(1)+1),NvarE,Norder(1),1, &
                        coeffE1,.5d0,  &
                        DofE_son,Norder(1),1)
        call transform2(DofE(1,naE(3)+1),NvarE,Norder(3),1, &
                        coeffE1,.5d0,  &
                        DofE_son,Norder(3),1)
!
!  ...first or third mid-edge son
      case(5,7)
        select case(Nson)
        case(5); j=1
        case(7); j=2
        end select
        coeffH1(1:nord1-1,1) = RRRH(1,1:nord1-1,3,1)
        coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,j,1:nord2-1)
        call transform2(Xnod(1,naH(5)+1),NDIMEN,nord1-1,nord2-1, &
                        coeffH1,coeffH2,  &
                        Xnod_son,1,nord2-1)
        call transform2(Xnod(1,naH(2)+1),NDIMEN,1,Norder(2)-1, &
                        .5d0,coeffH2,  &
                        Xnod_son,1,Norder(2)-1)
        call transform2(Xnod(1,naH(4)+1),NDIMEN,1,Norder(4)-1, &
                        .5d0,coeffH2,  &
                        Xnod_son,1,Norder(4)-1)
!
        call transform2(DofH(1,naH(5)+1),NvarH,nord1-1,nord2-1, &
                        coeffH1,coeffH2,  &
                        DofH_son,1,nord2-1)
        call transform2(DofH(1,naH(2)+1),NvarH,1,Norder(2)-1, &
                        .5d0,coeffH2,  &
                        DofH_son,1,Norder(2)-1)
        call transform2(DofH(1,naH(4)+1),NvarH,1,Norder(4)-1, &
                        .5d0,coeffH2,  &
                        DofH_son,1,Norder(4)-1)
!
        coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,j,1:nord2)
        call transform2(DofE(1,naE(5)+nord1*(nord2-1)+1),NvarE,nord1-1,nord2, &
                        coeffE1,coeffH2,  &
                        DofE_son,1,nord2)
        call transform2(DofE(1,naE(2)+1),NvarE,1,Norder(2), &
                        .5d0,coeffE2,  &
                        DofE_son,1,Norder(2))
        call transform2(DofE(1,naE(4)+1),NvarE,1,Norder(2), &
                        .5d0,coeffE2,  &
                        DofE_son,1,Norder(4))
!
!  ...vertex son
      case(9)
        coeffH1(1:nord1-1,1) = RRRH(1,1:nord1-1,3,1)
        coeffH2(1:nord2-1,1) = RRRH(1,1:nord2-1,3,1)
        call transform2(Xnod(1,naH(5)+1),NDIMEN,nord1-1,nord2-1, &
                        coeffH1,coeffH2,  &
                        Xnod_son,1,1)
        call transform2(Xnod(1,naH(1)+1),NDIMEN,Norder(1)-1,1, &
                        coeffH1,.5d0,  &
                        Xnod_son,1,1)
        call transform2(Xnod(1,naH(3)+1),NDIMEN,Norder(3)-1,1, &
                        coeffH1,.5d0,  &
                        Xnod_son,1,1)
        call transform2(Xnod(1,naH(2)+1),NDIMEN,1,Norder(2)-1, &
                        .5d0,coeffH2,  &
                        Xnod_son,1,1)
        call transform2(Xnod(1,naH(4)+1),NDIMEN,1,Norder(4)-1, &
                        .5d0,coeffH2,  &
                        Xnod_son,1,1)
        do j=1,4
          Xnod_son(1:NDIMEN,1) = Xnod_son(1:NDIMEN,1) &
                               + Xnod(1:NDIMEN,j)*.25d0
        enddo
!
        call transform2(DofH(1,naH(5)+1),NvarH,nord1-1,nord2-1, &
                        coeffH1,coeffH2,  &
                        DofH_son,1,1)
        call transform2(DofH(1,naH(1)+1),NvarH,Norder(1)-1,1, &
                        coeffH1,.5d0,  &
                        DofH_son,1,1)
        call transform2(DofH(1,naH(3)+1),NvarH,Norder(3)-1,1, &
                        coeffH1,.5d0,  &
                        DofH_son,1,1)
        call transform2(DofH(1,naH(2)+1),NvarH,1,Norder(2)-1, &
                        .5d0,coeffH2,  &
                        DofH_son,1,1)
        call transform2(DofH(1,naH(4)+1),NvarH,1,Norder(4)-1, &
                        .5d0,coeffH2,  &
                        DofH_son,1,1)
        do j=1,4
          DofH_son(1:NvarH,1) = DofH_son(1:NvarH,1) &
                                + DofH(1:NvarH,j)*.25d0
        enddo
      end select
!
!  ...anisotropic refinement across the first axis
      case(10)
        select case(Nson)
!
!  .....mid-face sons
        case(1,2)
          coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,Nson,1:nord1-1)
          coeffH2(1:nord2-1,1:nord2-1) = 0.d0
          do i=1,nord2-1
            coeffH2(i,i) = 1.d0
          enddo
          call transform2(Xnod(1,naH(5)+1),NDIMEN,nord1-1,nord2-1, &
                          coeffH1,coeffH2,  &
                          Xnod_son,nord1-1,nord2-1)
          call transform2(DofH(1,naH(5)+1),NvarH,nord1-1,nord2-1, &
                          coeffH1,coeffH2,  &
                          DofH_son,nord1-1,nord2-1)
          coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,Nson,1:nord1)
          coeffE2(1:nord2,1:nord2) = 0.d0
          do i=1,nord2
            coeffE2(i,i) = 1.d0
          enddo
          call transform2(DofE(1,naE(5)+1),NvarE,nord1,nord2-1, &
                          coeffE1,coeffH2,  &
                          DofE_son,nord1,nord2-1)
          call transform2(DofE(1,naE(5)+nord1*(nord2-1)),NvarE,nord1-1,nord2, &
                          coeffH1,coeffE2,  &
                          DofE_son(1,nord1*(nord2-1)+1),nord1-1,nord2)
          call transform2(DofV,NvarV,nord1,nord2, &
                          coeffE1,coeffE2,  &
                          DofV_son,nord1,nord2)
!
!  .....mid-edge son
        case(3)
          coeffH1(1:nord1-1,1) = RRRH(1,1:nord1-1,3,1)
          coeffH2(1:nord2-1,1:nord2-1) = 0.d0
          do i=1,nord2-1
            coeffH2(i,i) = 1.d0
          enddo
          call transform2(Xnod(1,naH(5)+1),NDIMEN,nord1-1,nord2-1, &
                          coeffH1,coeffH2,  &
                          Xnod_son,1,nord2-1)
          call transform2(Xnod(1,naH(2)+1),NDIMEN,1,Norder(2)-1, &
                          .5d0,coeffH2,  &
                          Xnod_son,1,Norder(2)-1)
          call transform2(Xnod(1,naH(4)+1),NDIMEN,1,Norder(4)-1, &
                          .5d0,coeffH2,  &
                          Xnod_son,1,Norder(4)-1)
!
          call transform2(DofH(1,naH(5)+1),NvarH,nord1-1,nord2-1, &
                          coeffH1,coeffH2,  &
                          DofH_son,1,nord2-1)
          call transform2(DofH(1,naH(2)+1),NvarH,1,Norder(2)-1, &
                          .5d0,coeffH2,  &
                          DofH_son,1,Norder(2)-1)
          call transform2(DofH(1,naH(4)+1),NvarH,1,Norder(4)-1, &
                          .5d0,coeffH2,  &
                          DofH_son,1,Norder(4)-1)
!
          coeffE2(1:nord2,1:nord2) = 0.d0
          do i=1,nord2
            coeffE2(i,i) = 1.d0
          enddo
          call transform2(DofE(1,naE(5)+nord1*(nord2-1)+1),NvarE,nord1-1,nord2, &
                          coeffE1,coeffH2,  &
                          DofE_son,1,nord2)
          call transform2(DofE(1,naE(2)+1),NvarE,1,Norder(2), &
                          .5d0,coeffE2,  &
                          DofE_son,1,Norder(2))
          call transform2(DofE(1,naE(4)+1),NvarE,1,Norder(2), &
                          .5d0,coeffE2,  &
                          DofE_son,1,Norder(4))
        end select
!
!  ...anisotropic refinement across the second axis
      case(01)
        select case(Nson)
!
!  .....mid-face sons
        case(1,2)
          coeffH1(1:nord1-1,1:nord1-1) = 0.d0
          do i=1,nord1-1
            coeffH1(i,i) = 1.d0
          enddo
          coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,Nson,1:nord2-1)
          call transform2(Xnod(1,naH(5)+1),NDIMEN,nord1-1,nord2-1, &
                          coeffH1,coeffH2,  &
                          Xnod_son,nord1-1,nord2-1)
          call transform2(DofH(1,naH(5)+1),NvarH,nord1-1,nord2-1, &
                          coeffH1,coeffH2,  &
                          DofH_son,nord1-1,nord2-1)
          coeffE1(1:nord1,1:nord1) = 0.d0
          do i=1,nord1
            coeffE1(i,i) = 1.d0
          enddo
          coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,Nson,1:nord2)
          call transform2(DofE(1,naE(5)+1),NvarE,nord1,nord2-1, &
                          coeffE1,coeffH2,  &
                          DofE_son,nord1,nord2-1)
          call transform2(DofE(1,naE(5)+nord1*(nord2-1)),NvarE,nord1-1,nord2, &
                          coeffH1,coeffE2,  &
                          DofE_son(1,nord1*(nord2-1)+1),nord1-1,nord2)
          call transform2(DofV,NvarV,nord1,nord2, &
                          coeffE1,coeffE2,  &
                          DofV_son,nord1,nord2)
!
!  .....mid-edge son
        case(3)
          coeffH1(1:nord1-1,1:nord1-1) = 0.d0
          do i=1,nord1-1
            coeffH1(i,i) = 1.d0
          enddo
          coeffH2(1:nord2-1,1) = RRRH(1,1:nord2-1,3,1)
          call transform2(Xnod(1,naH(5)+1),NDIMEN,nord1-1,nord2-1, &
                          coeffH1,coeffH2,  &
                          Xnod_son,nord1-1,1)
          call transform2(Xnod(1,naH(1)+1),NDIMEN,Norder(1)-1,1, &
                          coeffH1,.5d0,  &
                          Xnod_son,Norder(1)-1,1)
          call transform2(Xnod(1,naH(3)+1),NDIMEN,Norder(3)-1,1, &
                          coeffH1,.5d0,  &
                          Xnod_son,Norder(3)-1,1)
!
          call transform2(DofH(1,naH(5)+1),NvarH,nord1-1,nord2-1, &
                          coeffH1,coeffH2,  &
                          DofH_son,nord1-1,1)
          call transform2(DofH(1,naH(1)+1),NvarH,Norder(1)-1,1, &
                          coeffH1,.5d0,  &
                          DofH_son,Norder(1)-1,1)
          call transform2(DofH(1,naH(3)+1),NvarH,Norder(3)-1,1, &
                          coeffH1,.5d0,  &
                          DofH_son,Norder(3)-1,1)
!
          coeffE1(1:nord1,1:nord1) = 0.d0
          do i=1,nord1
            coeffE1(i,i) = 1.d0
          enddo
          call transform2(DofE(1,naE(5)+1),NvarE,nord1,nord2-1, &
                          coeffE1,coeffH2,  &
                          DofE_son,nord1,1)
          call transform2(DofE(1,naE(1)+1),NvarE,Norder(1),1, &
                          coeffE1,.5d0,  &
                          DofE_son,Norder(1),1)
          call transform2(DofE(1,naE(3)+1),NvarE,Norder(3),1, &
                          coeffE1,.5d0,  &
                          DofE_son,Norder(3),1)
        end select
!
!  ...anisotropic refinement across the first axis
      case(100)
        select case(Nson)
!
!  .....middle sons
        case(1,2)
          coeffH1(1:nord1-1,1:nord1-1) = RRRH(1,1:nord1-1,Nson,1:nord1-1)
          call setI(coeffH2,nord2-1)
          call setI(coeffH3,nord3-1)
          coeffE1(1:nord1,1:nord1) = RRRE(1,1:nord1,Nson,1:nord1)
          call setI(coeffE2,nord2)
          call setI(coeffE3,nord3)
          nod_type = 'mdlb'
!
!  .....mid-face son
        case(3)
          coeffH1(1:nord1-1,1) = RRRH(1,1:nord1-1,3,1)
          call setI(coeffH2,nord2-1)
          call setI(coeffH3,nord3-1)
          call setI(coeffE2,nord2)
          call setI(coeffE3,nord3)
          nod_type = 'mf23'
        end select
      end select
!
!  ...anisotropic refinement across the second axis
      case(010)
        select case(Nson)
!
!  .....middle sons
        case(1,2)
          call setI(coeffH1,nord1-1)
          coeffH2(1:nord2-1,1:nord2-1) = RRRH(1,1:nord2-1,Nson,1:nord2-1)
          call setI(coeffH3,nord3-1)
          call setI(coeffE1,nord1)
          coeffE2(1:nord2,1:nord2) = RRRE(1,1:nord2,Nson,1:nord2)
          call setI(coeffE3,nord3)
          nod_type = 'mdlb'
!
!  .....mid-face son
        case(3)
          call setI(coeffH1,nord1-1)
          call setI(coeffH3,nord3-1)
          coeffH2(1:nord2-1,1) = RRRH(1,1:nord2-1,3,1)
          call setI(coeffE1,nord1)
          call setI(coeffE3,nord3)
          nod_type = 'mf13'
        end select
      end select
!
!  ...anisotropic refinement across the third axis
      case(001)
        select case(Nson)
!
!  .....middle sons
        case(1,2)
          call setI(coeffH1,nord1-1)
          call setI(coeffH2,nord2-1)
          coeffH3(1:nord3-1,1:nord3-1) = RRRH(1,1:nord3-1,Nson,1:nord3-1)
          call setI(coeffE1,nord1)
          call setI(coeffE2,nord2)
          coeffE3(1:nord3,1:nord3) = RRRE(1,1:nord3,Nson,1:nord3)
          nod_type = 'mdlb'
!
!  .....mid-face son
        case(3)
          call setI(coeffH1,nord1-1)
          call setI(coeffH2,nord2-1)
          coeffH3(1:nord3-1,1) = RRRH(1,1:nord3-1,3,1)
          call setI(coeffE1,nord1)
          call setI(coeffE2,nord2)
          nod_type = 'mf12'
        end select
      end select
!
!  ...initiate dof for the node
      select case(nod_type}
!
!  ...a middle node
      case('mdlb')
        call transform3(Xnod(1,naH(19)+1),NDIMEN,nord1-1,nord2-1,nord3-1 &
                        coeffH1,coeffH2,coeffH3,  &
                        Xnod_son,nord1-1,nord2-1,nord3-1)
        call transform3(DofH(1,naH(19)+1),NvarH, nord1-1,nord2-1,nord3-1 &
                        coeffH1,coeffH2,coeffH3,  &
                        DofH_son,nord1-1,nord2-1,nord3-1)
!
        call transform3(DofE(1,naE(19)+1),NvarE, nord1,nord2-1,nord3-1 &
                        coeffE1,coeffH2,coeffH3,  &
                        DofE_son,nord1,nord2-1,nord3-1)
        ndofE = nord1*(nord2-1)*(nord3-1) 
        call transform3(DofE(1,naE(19)+ndofE+1),NvarE, nord1-1,nord2,nord3-1 &
                        coeffH1,coeffE2,coeffH3,  &
                        DofE_son(1,ndofE+1),nord1-1,nord2,nord3-1)
        ndofE = ndofE+(nord1-1)*nord2*(nord3-1) 
        call transform3(DofE(1,naE(19)+ndofE+1),NvarE, nord1-1,nord2-1,nord3 &
                        coeffH1,coeffH2,coeffE3,  &
                        DofE_son(1,ndofE+1),nord1-1,nord2-1,nord3)
!
        call transform3(DofV(1,naV(19)+1),NvarV, nord1-1,nord2,nord3 &
                        coeffH1,coeffE2,coeffE3,  &
                        DofV_son,nord1-1,nord2,nord3)
        ndofV = (nord1-1)*nord2*nord3 
        call transform3(DofV(1,naV(19)+ndofV+1),NvarV, nord1,nord2-1,nord3 &
                        coeffE1,coeffH2,coeffE3,  &
                        DofV_son(1,ndofV+1),nord1,nord2-1,nord3)
        ndofV = ndofV+nord1*(nord2-1)*nord3 
        call transform3(DofV(1,naV(19)+ndofV+1),NvarV, nord1,nord2,nord3-1 &
                        coeffE1,coeffE2,coeffH3,  &
                        DofV_son(1,ndofV+1),nord1,nord2,nord3-1)
!
        call transform3(DofQ,NvarQ, nord1,nord2,nord3 &
                        coeffE1,coeffE2,coeffE3,  &
                        DofQ_son,nord1,nord2,nord3)
!
!  ...a mid-face node parallel to axes 1,2
      case('mf12')
!
!  .....contributions from the middle node
        call transform3(Xnod(1,naH(19)+1),NDIMEN,nord1-1,nord2-1,nord3-1, &
                        coeffH1,coeffH2,coeffH3,  &
                        Xnod_son,nord1-1,nord2-1,1)
        call transform3(DofH(1,naH(19)+1),NvarH, nord1-1,nord2-1,nord3-1, &
                        coeffH1,coeffH2,coeffH3,  &
                        DofH_son,nord1-1,nord2-1,1)
!
        call transform3(DofE(1,naE(19)+1),NvarE, nord1,nord2-1,nord3-1, &
                        coeffE1,coeffH2,coeffH3,  &
                        DofE_son,nord1,nord2-1,1)
        ndofE = nord1*(nord2-1)*(nord3-1)
        call transform3(DofE(1,naE(19)+ndofE+1),NvarE, nord1-1,nord2,nord3-1, &
                        coeffH1,coeffE2,coeffH3,  &
                        DofE_son(1,nord1*(nord2-1)+1),nord1-1,nord2,1)
!
        call transform3(DofV(1,naV(19)+1),NvarV, nord1,nord2,nord3-1, &
                        coeffE1,coeffE2,coeffH3,  &
                        DofV_son,nord1,nord2,1)
!
!  .....contributions from the bottom mid-face node
        call decode(Norder(13), nordf1,nordf2)
        call transform3(Xnod(1,naH(13)+1),NDIMEN,nordf1-1,nordf2-1,1, &
                        coeffH1,coeffH2,.5d0,  &
                        Xnod_son,nord1-1,nord2-1,1)
        call transform3(DofH(1,naH(13)+1),NvarH, nordf1-1,nordf2-1,1, &
                        coeffH1,coeffH2,.5d0,  &
                        DofH_son,nord1-1,nord2-1,1)
!
        call transform3(DofE(1,naE(13)+1),NvarE, nordf1,nordf2-1,1, &
                        coeffE1,coeffH2,.5d0,  &
                        DofE_son,nord1,nord2-1,1)
        ndofE = nord1*(nord2-1)
        call transform3(DofE(1,naE(13)+ndofE+1),NvarE, nordf1-1,nordf2,1, &
                        coeffH1,coeffE2,.5d0,  &
                        DofE_son(1,ndofE+1),nord1-1,nord2,1)
!
        call transform3(DofV(1,naV(13)+1),NvarV, nordf1,nordf2,1, &
                        coeffE1,coeffE2,.5d0,  &
                        DofV_son,nord1,nord2,1)
!
!  .....contributions from the top mid-face node
        call decode(Norder(14), nordf1,nordf2)
        call transform3(Xnod(1,naH(14)+1),NDIMEN,nordf1-1,nordf2-1,1, &
                        coeffH1,coeffH2,.5d0,  &
                        Xnod_son,nord1-1,nord2-1,1)
        call transform3(DofH(1,naH(14)+1),NvarH, nordf1-1,nordf2-1,1, &
                        coeffH1,coeffH2,.5d0,  &
                        DofH_son,nord1-1,nord2-1,1)
!
        call transform3(DofE(1,naE(14)+1),NvarE, nordf1,nordf2-1,1, &
                        coeffE1,coeffH2,.5d0,  &
                        DofE_son,nord1,nord2-1,1)
        ndofE = nord1*(nord2-1)
        call transform3(DofE(1,naE(14)+ndofE+1),NvarE, nordf1-1,nordf2,1, &
                        coeffH1,coeffE2,.5d0,  &
                        DofE_son(1,ndofE+1),nord1-1,nord2,1)
!
        call transform3(DofV(1,naV(14)+1),NvarV, nordf1,nordf2,1, &
                        coeffE1,coeffE2,.5d0,  &
                        DofV_son,nord1,nord2,1)
!
!  ...a mid-face node parallel to axes 2,3
      case('mf23')
!
!  .....contributions from the middle node
        call transform3(Xnod(1,naH(19)+1),NDIMEN,nord1-1,nord2-1,nord3-1 &
                        coeffH1,coeffH2,coeffH3,  &
                        Xnod_son,1,nord2-1,nord3-1)
        call transform3(DofH(1,naH(19)+1),NvarH, nord1-1,nord2-1,nord3-1 &
                        coeffH1,coeffH2,coeffH3,  &
                        DofH_son,1,nord2-1,nord3-1)
!
        ndofE = nord1*(nord2-1)*(nord3-1)
        call transform3(DofE(1,naE(19)+ndofE+1),NvarE, nord1,nord2,nord3-1 &
                        coeffH1,coeffE2,coeffH3,  &
                        DofE_son,1,nord2,nord3-1)
        ndofE = ndofE+(nord1-1)*nord2*(nord3-1)
        call transform3(DofE(1,naE(19)+ndofE+1),NvarE, nord1-1,nord2-1,nord3 &
                        coeffH1,coeffH2,coeffE3,  &
                        DofE_son(1,nord2*(nord3-1)+1),1,nord2-1,nord3)
!
        ndofV = nord1*nord2*(nord3-1)
        call transform3(DofV(1,naV(19)+ndofV+1),NvarV, nord1-1,nord2,nord3 &
                        coeffH1,coeffE2,coeffE3,  &
                        DofV_son,1,nord2,nord3)
!
!  .....contributions from the left mid-face node
        call decode(Norder(18), nordf2,nordf3)
        call transform3(Xnod(1,naH(18)+1),NDIMEN,1,nordf2-1,nordf3-1, &
                        .5d0,coeffH2,coeffH3,  &
                        Xnod_son,1,nord2-1,nord3-1)
        call transform3(DofH(1,naH(18)+1),NvarH, 1,nordf2-1,nordf3-1, &
                        .5d0,coeffH2,coeffH3,  &
                        DofH_son,1,nord2-1,nord3-1)
!
        call transform3(DofE(1,naE(18)+1),NvarE, 1,nordf2,nordf3-1, &
                        .5d0,coeffE2,coeffH3,  &
                        DofE_son,1,nord2,nord3-1)
        ndofE = nord1*(nord2-1)
        call transform3(DofE(1,naE(18)+ndofE+1),NvarE, 1,nordf2-1,nordf3, &
                        .5d0,coeffH2,coeffE3,  &
                        DofE_son(1,ndofE+1),1,nord2-1,nord3)
!
        call transform3(DofV(1,naV(18)+1),NvarV, 1,nordf2,nordf3, &
                        .5d0,coeffE2,coeffE3,  &
                        DofV_son,1,nord2,nord3)
!
!  .....contributions from the right mid-face node
        call decode(Norder(16), nordf2,nordf3)
        call transform3(Xnod(1,naH(16)+1),NDIMEN,1,nordf2-1,nordf3-1, &
                        .5d0,coeffH2,coeffH3,  &
                        Xnod_son,1,nord2-1,nord3-1)
        call transform3(DofH(1,naH(16)+1),NvarH, 1,nordf2-1,nordf3-1, &
                        .5d0,coeffH2,coeffH3,  &
                        DofH_son,1,nord2-1,nord3-1)
!
        call transform3(DofE(1,naE(16)+1),NvarE, 1,nordf2,nordf3-1, &
                        .5d0,coeffE2,coeffH3,  &
                        DofE_son,1,nord2,nord3-1)
        ndofE = nord1*(nord2-1)
        call transform3(DofE(1,naE(16)+ndofE+1),NvarE, 1,nordf2-1,nordf3, &
                        .5d0,coeffH2,coeffE3,  &
                        DofE_son(1,ndofE+1),1,nord2-1,nord3)
!
        call transform3(DofV(1,naV(16)+1),NvarV, 1,nordf2,nordf3, &
                        .5d0,coeffE2,coeffE3,  &
                        DofV_son,1,nord2,nord3)
!
!  ...a mid-face node parallel to axes 1,3
      case('mf13')
!
!  .....contributions from the middle node
        call transform3(Xnod(1,naH(19)+1),NDIMEN,nord1-1,nord2-1,nord3-1 &
                        coeffH1,coeffH2,coeffH3,  &
                        Xnod_son,nord1-1,1,nord3-1)
        call transform3(DofH(1,naH(19)+1),NvarH, nord1-1,nord2-1,nord3-1 &
                        coeffH1,coeffH2,coeffH3,  &
                        DofH_son,nord1-1,1,nord3-1)
!
        ndofE = 0
        call transform3(DofE(1,naE(19)+ndofE+1),NvarE, nord1,nord2,nord3-1 &
                        coeffE1,coeffH2,coeffH3,  &
                        DofE_son,nord1,1,nord3-1)
        ndofE = ndofE+nord1*(nord2-1)*(nord3-1)+(nord1-1)*nord2*(nord3-1)
        call transform3(DofE(1,naE(19)+ndofE+1),NvarE, nord1-1,nord2-1,nord3 &
                        coeffE1,coeffH2,coeffE3,  &
                        DofE_son(1,nord1*(nord3-1)+1),nord1-1,1,nord3)
!
        ndofV = nord1*nord2*(nord3-1)+(nord1-1)*nord2*nord3
        call transform3(DofV(1,naV(19)+ndofV+1),NvarV, nord1,nord2-1,nord3 &
                        coeffE1,coeffH2,coeffE3,  &
                        DofV_son,nord1,1,nord3)
!
!  .....contributions from the front mid-face node
        call decode(Norder(18), nordf2,nordf3)
        call transform3(Xnod(1,naH(18)+1),NDIMEN,1,nordf2-1,nordf3-1, &
                        .5d0,coeffH2,coeffH3,  &
                        Xnod_son,1,nord2-1,nord3-1)
        call transform3(DofH(1,naH(18)+1),NvarH, 1,nordf2-1,nordf3-1, &
                        .5d0,coeffH2,coeffH3,  &
                        DofH_son,1,nord2-1,nord3-1)
!
        call transform3(DofE(1,naE(18)+1),NvarE, 1,nordf2,nordf3-1, &
                        .5d0,coeffE2,coeffH3,  &
                        DofE_son,1,nord2,nord3-1)
        ndofE = nord1*(nord2-1)
        call transform3(DofE(1,naE(18)+ndofE+1),NvarE, 1,nordf2-1,nordf3, &
                        .5d0,coeffH2,coeffE3,  &
                        DofE_son(1,ndofE+1),1,nord2-1,nord3)
!
        call transform3(DofV(1,naV(18)+1),NvarV, 1,nordf2,nordf3, &
                        .5d0,coeffE2,coeffE3,  &
                        DofV_son,1,nord2,nord3)
!
!  .....contributions from the back mid-face node
        call decode(Norder(16), nordf2,nordf3)
        call transform3(Xnod(1,naH(16)+1),NDIMEN,1,nordf2-1,nordf3-1, &
                        .5d0,coeffH2,coeffH3,  &
                        Xnod_son,1,nord2-1,nord3-1)
        call transform3(DofH(1,naH(16)+1),NvarH, 1,nordf2-1,nordf3-1, &
                        .5d0,coeffH2,coeffH3,  &
                        DofH_son,1,nord2-1,nord3-1)
!
        call transform3(DofE(1,naE(16)+1),NvarE, 1,nordf2,nordf3-1, &
                        .5d0,coeffE2,coeffH3,  &
                        DofE_son,1,nord2,nord3-1)
        ndofE = nord1*(nord2-1)
        call transform3(DofE(1,naE(16)+ndofE+1),NvarE, 1,nordf2-1,nordf3, &
                        .5d0,coeffH2,coeffE3,  &
                        DofE_son(1,ndofE+1),1,nord2-1,nord3)
!
        call transform3(DofV(1,naV(16)+1),NvarV, 1,nordf2,nordf3, &
                        .5d0,coeffE2,coeffE3,  &
                        DofV_son,1,nord2,nord3)
!
      end select

!
!
      end subroutine initiate_mdlq_dof


      subroutine transform3(Uu,M,Iu,Ju,Ku,Coeff1,Coeff2,Coeff3, Vv,Iv,Jv,Kv)
      use parameters, only: ZERO
      implicit none
!
      integer,                      intent(in)  :: M,Iu,Ju,Ku,Iv,Jv,Kv
      VTYPE, dimension(M,Iu,Ju,Ku), intent(in)  :: Uu
      VTYPE, dimension(Iu,Iv),      intent(in)  :: Coeff1
      VTYPE, dimension(Ju,Jv),      intent(in)  :: Coeff2
      VTYPE, dimension(Ku,Kv),      intent(in)  :: Coeff3
      VTYPE, dimension(M,Iv,Jv,Kv), intent(out) :: Vv
!
      integer :: i,i1,j,j1,k,k1
      VTYPE, dimension(M,Iv,Iu,Ku)  :: aa
      VTYPE, dimension(M,Iv,IV,Ku)  :: bb
!
      aa = ZERO
      do i=1,Iu
        do i1=1,Iv
          aa(1:M,i1,1:Ju,1:Ku) = aa(1:M,i1,1:Ju,1:Ku) &
        enddo                  + Uu(1:M,i,1:Ju,1:Ku)*Coeff1(i,i1)
      enddo
      bb = ZERO
      do j=1,Ju
        do j1=1,Jv
          bb(1:M,1:Iv,j1,1:Ku) = bb(1:M,1:Iv,j1,1:Ku) &
                               + aa(1:M,1:Iv,j,1:Ku)*Coeff2(j,j1)
        enddo
      enddo
      do k=1,Ku
        do k1=1,Kv
          Vv(1:M,1:Iv,1:Jv,k1) = Vv(1:M,1:Iv,1:Jv,k1) &
                               + bb(1:M,1:Iv,1:Jv,k)*Coeff3(k,k1)
        enddo
      enddo
!
      end subroutine transform3


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
      VTYPE, dimension(Na+1,Na-1), intent(out) :: A
!
      integer :: i
!
      A = ZERO
      do i=1,Na-1
        A(2+i,i) = 1.d0
      enddo
!
      end subroutine setI0


      rh(1:2) = .5d0
      rh(3:MAXP+1) = RRRH(1,1:MAXP-1,3,1)
      

