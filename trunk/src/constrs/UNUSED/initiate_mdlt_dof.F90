#include "typedefs.h"
!-------------------------------------------------------------------------
!> Purpose : initiate H1, H(curl) and H(div), i.e. L2 dof for the sons of 
!            a mid-face triangular node
!
!> @date May 20
!
!> @param[in]  Nson    = 1,...,7 (the son number)
!> @param[in]  Norder  - polynomial order for the face
!> @param[in]  NvarH,NvarE,NvarV - number of H1, H(curl) and H(div) components
!> @param[in]  Xnod    - geometry dof for the face
!> @param[in]  DofH    - H1 dof for the face
!> @param[in]  DofE    - H(curl) dof for the face
!> @param[in]  DofV    - H(div) dof for the face
!> @param[in]  NdofH,NdofE,NdofV - number of dof for the son
!
!> @param[out] Xnod_son - geometry dof for the son
!> @param[out] DofH_son - H1 dof for the son
!> @param[out] DofE_son - H(curl) dof for the son
!> @param[out] DofV_son - H(div) dof for the son
!-------------------------------------------------------------------------
!
      subroutine initiate_mdlt_dof(Nson,Norder,NvarH,NvarE,NvarV, &
                                   Xnod,DofH,DofE,DofV, &
                                   Xnod_son,DofH_son,DofE_son,DofV_son, &
                                   NdofH,NdofE,NdofV)
!
      use constraints
      implicit none
!
      integer, intent(in)  :: Nson,NvarH,NvarE,NvarV,NdofH,NdofE,NdofV    
      integer,dimension(4),                intent(in)  :: Norder
      real*8, dimension(NDIMEN,MAXtriaH),  intent(in)  :: Xnod
      VTYPE,  dimension(NvarH, MAXtriaH),  intent(in)  :: DofH
      VTYPE,  dimension(NvarE, MAXtriaE),  intent(in)  :: DofE
      VTYPE,  dimension(NvarV, MAXtriaV),  intent(in)  :: DofV
!
      real*8, dimension(NDIMEN,NdofH),  intent(out) :: Xnod_son
      VTYPE,  dimension(NvarH, NdofH),  intent(out) :: DofH_son
      VTYPE,  dimension(NvarE, NdofE),  intent(out) :: DofE_son
      VTYPE,  dimension(NvarV, NdofV),  intent(out) :: DofV_son
!
!  ...locals
      integer :: nord,i,nrdofH,nrdofE, &
                 ndofHf,ndofEf,ndofVf,nvoid
      integer,dimension(4)                 :: naH,naE
      real*8, dimension(MAXmdltH,MAXmdltH) :: coeffH
      real*8, dimension(MAXmdltE,MAXmdltE) :: coeffE
      real*8, dimension(MAXmdltV,MAXmdltV) :: coeffV
!
!
      Xnod_son = 0.d0; DofH_son = ZERO; DofE_son = ZERO; DofV_son = ZERO
!
      nord = Norder(4)
!
      nrdofH=3; nrdofE=0
      do i=1,4
        naH(i) = nrdofH; naE(i) = nrdofE;
        nrdofH = nrdofH + Norder(i)-1; nrdofE = nrdofE + Norder(i)
      enddo
!  
!  ...mid-face parent node
      call ndof_nod('mdlt',Norder(4), ndofHf,ndofEf,ndofVf,nvoid)
      coeffH(1:ndofH,1:ndofH) = RRTH(1,1:ndofHf,Nson,1:ndofH)
      call Rtransform1(Xnod(1,naH(4)+1),NDIMEN,ndofHf,coeffH, Xnod_son,ndofH)
      call Vtransform1(DofH(1,naH(4)+1),NvarH, ndofHf,coeffH, DofH_son,ndofH)
      coeffE(1:ndofE,1:ndofE) = RRTE(1,1:ndofEf,Nson,1:ndofE)
      call Rtransform1(DofE(1,naE(4)+1),NvarE, ndofEf,coeffE, DofE_son,ndofE)
      coeffV(1:ndofV,1:ndofV) = RRTQ(1:ndofVf,Nson,1:ndofV)
      call Vtransform1(DofV,            NvarV, ndofVf,coeffV, DofV_son,ndofV)
!
!  ...mid-edge parent nodes
      do i=1,3
        call ndof_nod('medg',Norder(i), ndofHf,ndofEf,nvoid,nvoid)
        coeffH(1:ndofH,1:ndofH) = RRTH(1+i,1:ndofHf,Nson,1:ndofH)
        call transform1(Xnod(1,naH(i)+1),NDIMEN,ndofHf,coeffH, Xnod_son,ndofH)
        call transform1(DofH(1,naH(i)+1),NvarH, ndofHf,coeffH, DofH_son,ndofH)
        coeffE(1:ndofE,1:ndofE) = RRTE(1+i,1:ndofEf,Nson,1:ndofE)
        call transform1(DofE(1,naE(i)+1),NvarE, ndofEf,coeffE, DofE_son,ndofE)
      enddo
!
!
      end subroutine initiate_mdlt_dof
