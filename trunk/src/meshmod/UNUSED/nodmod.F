
#include "typedefs.h"
c
c---------------------------------------------------------------------
c
c   routine name       - nodmod
c
c---------------------------------------------------------------------
c
c   latest revision    - May 20
c
c   purpose            - routine modifies order of approximation
c                        for a node
c
c   arguments :
c
c     in:
c             Nod      - node number
c             Nordn    - new order of approximation
c
c-----------------------------------------------------------------------
c
      subroutine nodmod(Nod,Nordn)
c 
      use data_structure3D
      implicit none
      common /ccopy_dofG/ iprint_copy_dofG
      common /ccopy_dofH/ iprint_copy_dofH
      common /ccopy_dofE/ iprint_copy_dofE
      common /ccopy_dofV/ iprint_copy_dofV
      common /ccopy_dofQ/ iprint_copy_dofQ
      integer :: iprint_copy_dofG,iprint_copy_dofH,iprint_copy_dofE,
     .           iprint_copy_dofV,iprint_copy_dofQ
c
      integer                              :: Nod,Nordn
      real*8, allocatable   ::  xnod(:,:)
      VTYPE,  allocatable   :: zdofH(:,:)
      VTYPE,  allocatable   :: zdofE(:,:)
      VTYPE,  allocatable   :: zdofV(:,:)
      VTYPE,  allocatable   :: zdofQ(:,:)
c
      character(len=4) :: ntype
      integer          :: iprint,icase,nordo,ndofH,ndofE,ndofV,ndofQ,
     .                    ndofHo,ndofEo,ndofVo,ndofQo,
     .                    nvarH,nvarE,nvarV,nvarQ
c
      select case(Nod)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
      iprint_copy_dofG=iprint
      iprint_copy_dofH=iprint; iprint_copy_dofE=iprint
      iprint_copy_dofV=iprint; iprint_copy_dofQ=iprint
c
c  ...node case
      icase = NODES(Nod)%case
      ntype = NODES(Nod)%type
      nordo = NODES(Nod)%order
      if (iprint.eq.1) then
        write(*,7010) Nod,ntype,nordo,icase
 7010   format('nodmod: Nod,type,order,icase = ',i6,2x,a5,2x,2i3)
        write(*,7013) NREQNH(icase),NREQNE(icase),
     .                NREQNV(icase),NREQNQ(icase)
 7013   format('nodmod: NREQNH,NREQNE,NREQNV,NREQNQ(icase) = ',4i2)
      endif
      if (Nordn.eq.nordo) return
c
c  ...determine the current number of dof
      call ndof_nod(ntype,nordo,
     .              ndofH,ndofE,ndofV,ndofQ)
c
c  ...save current gdofs in local array and update number of gdofs
      if (Associated(NODES(Nod)%dof%coord)) then
        allocate(xnod(NDIMEN,ndofH))
        xnod(1:NDIMEN,1:ndofH) = NODES(Nod)%dof%coord(1:NDIMEN,1:ndofH)
        deallocate(NODES(Nod)%dof%coord)
      endif
c      
c  ...save current dof in local arrays
      if (Associated(NODES(Nod)%dof%zdofH)) then
        nvarH = NREQNH(icase)*NRCOMS
        allocate(zdofH(nvarH,ndofH))
        zdofH(1:nvarH,1:ndofH) = NODES(Nod)%dof%zdofH(1:nvarH,1:ndofH)
        deallocate(NODES(Nod)%dof%zdofH)
        NRDOFSH = NRDOFSH - ndofH*NREQNH(icase)
      endif
      if (Associated(NODES(Nod)%dof%zdofE)) then
        nvarE = NREQNE(icase)*NRCOMS
        allocate(zdofE(nvarE,ndofE))
        zdofE(1:nvarE,1:ndofE) = NODES(Nod)%dof%zdofE(1:nvarE,1:ndofE)
        deallocate(NODES(Nod)%dof%zdofE)
        NRDOFSE = NRDOFSE - ndofE*NREQNE(icase)
      endif
      if (Associated(NODES(Nod)%dof%zdofV)) then
        nvarV = NREQNV(icase)*NRCOMS
        allocate(zdofV(nvarV,ndofV))
        zdofV(1:nvarV,1:ndofV) = NODES(Nod)%dof%zdofV(1:nvarV,1:ndofV)
        deallocate(NODES(Nod)%dof%zdofV)
        NRDOFSV = NRDOFSV - ndofV*NREQNV(icase)
      endif
      if (Associated(NODES(Nod)%dof%zdofQ)) then
        nvarQ = NREQNQ(icase)*NRCOMS
        allocate(zdofQ(nvarQ,ndofQ))
        zdofQ(1:nvarQ,1:ndofQ) = NODES(Nod)%dof%zdofQ(1:nvarQ,1:ndofQ)
        deallocate(NODES(Nod)%dof%zdofQ)
        NRDOFSQ = NRDOFSQ - ndofQ*NREQNQ(icase)
      endif
c
c  ...save number of dof corresponding to the old order
      ndofHo = ndofH
      ndofEo = ndofE
      ndofVo = ndofV
      ndofQo = ndofQ
c
c  ...modify the order of approximation 
      NODES(Nod)%order = Nordn
c
c  ...calculate the new number of dof for the node
      call ndof_nod(ntype,Nordn,
     .              ndofH,ndofE,ndofV,ndofQ)
      if (iprint.eq.1) then
        write(*,7015) ndofH,ndofE,ndofV,ndofQ
 7015   format('nodmod: ndofH,ndofE,ndofV,ndofQ = ',4i5)
      endif
c
c  ...allocate memory for geometry and solution dof
c     and copy old dof from the local arrays
      if (ndofH.gt.0) then
        allocate(NODES(Nod)%dof%coord(NDIMEN,ndofH))
        if (ndofHo.gt.0) then
          call copy_dofG(ntype,nordo,Nordn,ndofHo,ndofH,
     .                   xnod,NODES(Nod)%dof%coord)
          deallocate(xnod)
        else
          NODES(Nod)%dof%coord=0.d0
        endif
      endif
      if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
        nvarH = NREQNH(icase)*NRCOMS
        allocate( NODES(Nod)%dof%zdofH(nvarH, ndofH))
        NRDOFSH = NRDOFSH + ndofH*NREQNH(icase)
        if (ndofHo.gt.0) then
          call copy_dofH(ntype,nordo,Nordn,ndofHo,ndofH,
     .                   nvarH,zdofH,NODES(Nod)%dof%zdofH)
          deallocate(zdofH)
        else
          NODES(Nod)%dof%zdofH = ZERO
        endif
      endif
      if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
        nvarE = NREQNE(icase)*NRCOMS
        allocate( NODES(Nod)%dof%zdofE(nvarE, ndofE))
        NRDOFSE = NRDOFSE + ndofE*NREQNE(icase)
        if (ndofEo.gt.0) then
          call copy_dofE(ntype,nordo,Nordn,ndofEo,ndofE,
     .                   nvarE,zdofE,NODES(Nod)%dof%zdofE)
          deallocate(zdofE)
        else
          NODES(Nod)%dof%zdofE = ZERO
        endif
      endif
      if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
        nvarV = NREQNV(icase)*NRCOMS
        allocate( NODES(Nod)%dof%zdofV(nvarV, ndofV))
        NRDOFSV = NRDOFSV + ndofV*NREQNV(icase)
        if (ndofVo.gt.0) then
          call copy_dofV(ntype,nordo,Nordn,ndofVo,ndofV,
     .                   nvarV,zdofV,NODES(Nod)%dof%zdofV)
          deallocate(zdofV)
        else
          NODES(Nod)%dof%zdofV = ZERO
        endif
      endif
      if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
        nvarQ = NREQNQ(icase)*NRCOMS
        allocate( NODES(Nod)%dof%zdofQ(nvarQ, ndofQ))
        NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(icase)
        if (ndofQo.gt.0) then
          call copy_dofQ(ntype,nordo,Nordn,ndofQo,ndofQ,
     .                   nvarQ,zdofQ,NODES(Nod)%dof%zdofQ)
          deallocate(zdofQ)
        else
          NODES(Nod)%dof%zdofQ = ZERO
        endif
      endif
c
c
      if (iprint.eq.1) then
        write(*,7020) Nod
 7020   format('nodmod: Nod = ',i6,' HAS BEEN UPDATED ')
        write(*,7030) nordo,Nordn
 7030   format('        OLD ORDER = ',i3,' NEW ORDER = ',i3)
        call pause
      endif
c
c
      end subroutine nodmod
c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofG
c
c---------------------------------------------------------------------
c
c   latest revision    - May 20
c
c   purpose            - routine copies G1 dof for a node of order
c                        Nordo to a node of order Nordn
c
c   arguments :
c
c     in:
c        Ntype         - node type
c        Nordo,Nordn   - old and new orders of the node
c        Xnodo,Xnodn - old and new dof
c        NdofGo,NdofGn - old and new # dof
c
c-----------------------------------------------------------------------
c
      subroutine copy_dofG(Ntype,Nordo,Nordn,NdofGo,NdofGn,
     .                     Xnodo,Xnodn)
c
      use parameters, only: NDIMEN
      implicit none
      common /ccopy_dofG/ iprint
      integer :: iprint
c
      character(len=4)            :: Ntype
      integer                     :: Nordo,Nordn,NdofGo,NdofGn
      real*8, dimension(NDIMEN,NdofGo) :: Xnodo
      real*8, dimension(NDIMEN,NdofGn) :: Xnodn
c
      integer :: nord1o,nord2o,nord3o,
     .           nord1n,nord2n,nord3n
c
      if (iprint.eq.1) then
        write(*,7010) Ntype,Nordo,Nordn,NdofGo,NdofGn
 7010   format(' copy_dofG: Ntype,Nordo,Nordn,NdofGo,NdofGn = ',
     .                      a4,2i4,2i6)
      endif
c
      select case(Ntype)
c
c  ...edge, triangle, tet, pyramid
      case('medg','mdlt','mdln','mdld')
        call copy_1array(NDIMEN,Xnodo,NdofGo, Xnodn,NdofGn)
c
c  ...quad
      case('mdlq')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
        call copy_2array(NDIMEN,Xnodo,nord1o-1,nord2o-1, 
     .                          Xnodn,nord1n-1,nord2n-1)
c
c  ...prism
      case('mdlp')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
        call copy_2array(NDIMEN,Xnodo,(nord1o-2)*(nord1o-1)/2,nord2o-1, 
     .                          Xnodn,(nord1n-2)*(nord1n-1)/2,nord2n-1)
c
c  ...hexa
      case('mdlb')
        call ddecode(Nordo, nord1o,nord2o,nord3o)
        call ddecode(Nordn, nord1n,nord2n,nord3n)
        call copy_3array(NDIMEN,Xnodo,nord1o-1,nord2o-1,nord3o-1,
     .                          Xnodn,nord1n-1,nord2n-1,nord3n-1)
c
      end select
      if (iprint.eq.1) write(*,*) 'copy_dofG: DONE'
c
c
      end subroutine copy_dofG
c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofH
c
c---------------------------------------------------------------------
c
c   latest revision    - May 20
c
c   purpose            - routine copies H1 dof for a node of order
c                        Nordo to a node of order Nordn
c
c   arguments :
c
c     in:
c        Ntype         - node type
c        Nordo,Nordn   - old and new orders of the node
c        NdofHo,NdofHn - old and new # dof
c        NvarH         - first dimension of ZdofHo,ZdofHn
c        ZdofHo,ZdofHn - old and new dof
c
c-----------------------------------------------------------------------
c
      subroutine copy_dofH(Ntype,Nordo,Nordn,NdofHo,NdofHn,
     .                     NvarH,ZdofHo,ZdofHn)
c
      use parameters, only: MAXEQNH
      implicit none
      common /ccopy_dofH/ iprint
      integer :: iprint
c
      character(len=4)                :: Ntype
      integer                         :: Nordo,Nordn,NdofHo,NdofHn,NvarH
      VTYPE, dimension (NvarH,NdofHo) :: ZdofHo
      VTYPE, dimension (NvarH,NdofHn) :: ZdofHn
c
      integer :: nord1o,nord2o,nord3o,
     .           nord1n,nord2n,nord3n
c
      if (iprint.eq.1) then
        write(*,7010) Ntype,Nordo,Nordn,NdofHo,NdofHn
 7010   format(' copy_dofH: Ntype,Nordo,Nordn,NdofHo,NdofHn = ',
     .                      a4,2i4,2i6)
      endif
c
      select case(Ntype)
c
c  ...edge, triangle, tet, pyramid
      case('medg','mdlt','mdln','mdld')
        call copy_1array(NvarH,ZdofHo,NdofHo, ZdofHn,NdofHn)
c
c  ...quad
      case('mdlq')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
        call copy_2array(NvarH,ZdofHo,nord1o-1,nord2o-1, 
     .                         ZdofHn,nord1n-1,nord2n-1)
c
c  ...prism
      case('mdlp')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
        call copy_2array(NvarH,
     .                   ZdofHo,(nord1o-2)*(nord1o-1)/2,nord2o-1, 
     .                   ZdofHn,(nord1n-2)*(nord1n-1)/2,nord2n-1)
c
c  ...hexa
      case('mdlb')
        call ddecode(Nordo, nord1o,nord2o,nord3o)
        call ddecode(Nordn, nord1n,nord2n,nord3n)
        call copy_3array(NvarH,ZdofHo,nord1o-1,nord2o-1,nord3o-1,
     .                         ZdofHn,nord1n-1,nord2n-1,nord3n-1)
c
      end select
      if (iprint.eq.1) write(*,*) 'copy_dofH: DONE'
c
c
      end subroutine copy_dofH
c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofE
c
c---------------------------------------------------------------------
c
c   latest revision    - May 20
c
c   purpose            - routine copies H(curl) dof for a node of order
c                        Nordo to a node of order Nordn
c
c   arguments :
c
c     in:
c        Ntype         - node type
c        Nordo,Nordn   - old and new orders of the node
c        NdofEo,NdofEn - old and new # dof
c        NvarE         - first dimension of ZdofEo,ZdofEn
c        ZdofEo,ZdofEn - old and new dof
c
c-----------------------------------------------------------------------
c
      subroutine copy_dofE(Ntype,Nordo,Nordn,NdofEo,NdofEn,
     .                     NvarE,ZdofEo,ZdofEn)
c
      use parameters, only: MAXEQNE
      implicit none
      common /ccopy_dofE/ iprint
      integer :: iprint
c
      character(len=4)               :: Ntype
      integer                        :: Nordo,Nordn,NdofEo,NdofEn,NvarE
      VTYPE, dimension(NvarE,NdofEo) :: ZdofEo
      VTYPE, dimension(NvarE,NdofEn) :: ZdofEn
c
      integer :: nord1o,nord2o,nord3o,
     .           nord1n,nord2n,nord3n,
     .           ibego,ibegn
c
      if (iprint.eq.1) then
        write(*,7010) Ntype,Nordo,Nordn,NdofEo,NdofEn
 7010   format(' copy_dofE: Ntype,Nordo,Nordn,NdofEo,NdofEn = ',
     .                      a4,2i4,2i6)
      endif
c
      select case(Ntype)
c
c  ...segment
      case('medg')
        call copy_1array(NvarE,ZdofEo,NdofEo, ZdofEn,NdofEn)
c
c  ...triangle
      case('mdlt')
c
c  .....the 2 families are intertwined to form a hierarchy
        call copy_1array(NvarE,ZdofEo,NdofEo, ZdofEn,NdofEn)
c
c  ...quadrilateral
      case('mdlq')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
        ibego=1; ibegn=1
c
c  .....family 1
        call copy_2array(NvarE,ZdofEo(1,ibego),nord1o,nord2o-1, 
     .                         ZdofEn(1,ibegn),nord1n,nord2n-1)
        ibego=ibego+nord1o*(nord2o-1)
        ibegn=ibegn+nord1n*(nord2n-1)
c
c  .....family 2
        call copy_2array(NvarE,ZdofEo(1,ibego),nord1o-1,nord2o, 
     .                         ZdofEn(1,ibegn),nord1n-1,nord2n)
        ibego=ibego+(nord1o-1)*nord2o - 1
        ibegn=ibegn+(nord1n-1)*nord2n - 1
c
        if ((ibegn.ne.NdofEn).or.(ibego.ne.NdofEo)) then
          write(*,*) 'copy_dofE: INCONSISTENCY 1'
          stop 1
        endif
c
c  ...tetrahedron
      case('mdln')
c
c  .....the 3 families are intertwined to form a hierarchy
        call copy_1array(NvarE,ZdofEo,NdofEo, ZdofEn,NdofEn)
c
c  ...prism
      case('mdlp')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
        ibego=1; ibegn=1
c
c  .....first 2 families are (horizontally) intertwined 
        call copy_2array(NvarE,
     .                   ZdofEo(1,ibego),nord1o*(nord1o-1),nord2o-1, 
     .                   ZdofEn(1,ibegn),nord1n*(nord1n-1),nord2n-1)
        ibego=ibego+nord1o*(nord1o-1)*(nord2o-1)
        ibegn=ibegn+nord1n*(nord1n-1)*(nord2n-1)
c
c  .....3rd family
        call copy_2array(NvarE,
     .                   ZdofEo(1,ibego),(nord1o-2)*(nord1o-1)/2,nord2o, 
     .                   ZdofEn(1,ibegn),(nord1n-2)*(nord1n-1)/2,nord2n)
        ibego=ibego+(nord1o-2)*(nord1o-1)/2*nord2o - 1
        ibegn=ibegn+(nord1n-2)*(nord1n-1)/2*nord2n - 1
c
        if ((ibegn.ne.NdofEn).or.(ibego.ne.NdofEo)) then
          write(*,*) 'ibegn,NdofEn,ibego,NdofEo = ',
     .                ibegn,NdofEn,ibego,NdofEo
          write(*,*) 'copy_dofE: INCONSISTENCY 2'
          stop 1
        endif
c
c  ...hexahedron
      case('mdlb')
        call ddecode(Nordo, nord1o,nord2o,nord3o)
        call ddecode(Nordn, nord1n,nord2n,nord3n)
        ibego=1; ibegn=1
c
c  .....family 1
        call copy_3array(NvarE,
     .                   ZdofEo(1,ibego),nord1o,nord2o-1,nord3o-1,
     .                   ZdofEn(1,ibegn),nord1n,nord2n-1,nord3n-1)
        ibego=ibego+nord1o*(nord2o-1)*(nord3o-1)
        ibegn=ibegn+nord1n*(nord2n-1)*(nord3n-1)
c
c  .....family 2
        call copy_3array(NvarE,
     .                   ZdofEo(1,ibego),nord1o-1,nord2o,nord3o-1,
     .                   ZdofEn(1,ibegn),nord1n-1,nord2n,nord3n-1)
        ibego=ibego+(nord1o-1)*nord2o*(nord3o-1)
        ibegn=ibegn+(nord1n-1)*nord2n*(nord3n-1)
c
c  .....family 3
        call copy_3array(NvarE,
     .                   ZdofEo(1,ibego),nord1o-1,nord2o-1,nord3o,
     .                   ZdofEn(1,ibegn),nord1n-1,nord2n-1,nord3n)
        ibego=ibego+(nord1o-1)*(nord2o-1)*nord3o - 1
        ibegn=ibegn+(nord1n-1)*(nord2n-1)*nord3n - 1
c
        if ((ibegn.ne.NdofEn).or.(ibego.ne.NdofEo)) then
          write(*,*) 'ibegn,NdofEn,ibego,NdofEo = ',
     .                ibegn,NdofEn,ibego,NdofEo
          write(*,*) 'copy_dofE: INCONSISTENCY 3'
          stop 1
        endif
c
c  ...pyramid
      case('mdld')
        ibego=1; ibegn=1
c
c  .....FAMILY 1 (gradients of H1 bubbles) 
        call copy_3array(NvarE,
     .                   ZdofEo(1,ibego),Nordo-1,Nordo-1,Nordo-1,
     .                   ZdofEn(1,ibegn),Nordn-1,Nordn-1,Nordn-1)
        ibego=ibego+(Nordo-1)**3
        ibegn=ibegn+(Nordn-1)**3
c
c  .....FAMILY 2 AND 3 (induced from quad face functions)
        call copy_3array(NvarE,ZdofEo(1,ibego),Nordo,Nordo-1,Nordo-1,
     .                         ZdofEn(1,ibegn),Nordn,Nordn-1,Nordn-1)
        ibego=ibego+Nordo*(Nordo-1)**2
        ibegn=ibegn+Nordn*(Nordn-1)**2
        call copy_3array(NvarE,ZdofEo(1,ibego),Nordo-1,Nordo,Nordo-1,
     .                         ZdofEn(1,ibegn),Nordn-1,Nordn,Nordn-1)
        ibego=ibego+Nordo*(Nordo-1)**2
        ibegn=ibegn+Nordn*(Nordn-1)**2
c
c  .....FAMILY 4 
        call copy_2array(NvarE,
     .                   ZdofEo(1,ibego),Nordo-1,Nordo-1, 
     .                   ZdofEn(1,ibegn),Nordn-1,Nordn-1)
        ibego=ibego+(Nordo-1)**2 - 1
        ibegn=ibegn+(Nordn-1)**2 - 1
c
        if ((ibegn.ne.NdofEn).or.(ibego.ne.NdofEo)) then
          write(*,*) 'copy_dofE: INCONSISTENCY 4'
          stop 1
        endif
c
      end select
      if (iprint.eq.1) write(*,*) 'copy_dofE: DONE'
c
c    
      end subroutine copy_dofE
c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofV
c
c---------------------------------------------------------------------
c
c   latest revision    - May 20
c
c   purpose            - routine copies H(div) dof for a node of order
c                        Nordo to a node of order Nordn
c
c   arguments :
c
c     in:
c        Ntype         - node type
c        Nordo,Nordn   - old and new orders of the node
c        NdofEo,NdofEn - old and new # dof
c        NvarV         - first dimension of ZdofVo,ZdofVn
c        ZdofEo,ZdofEn - old and new dof
c
c-----------------------------------------------------------------------
c
      subroutine copy_dofV(Ntype,Nordo,Nordn,NdofVo,NdofVn,
     .                     NvarV,ZdofVo,ZdofVn)
c
      use parameters, only: MAXEQNV
      implicit none
      common /ccopy_dofV/ iprint
      integer :: iprint
c
      character(len=4)               :: Ntype
      integer                        :: Nordo,Nordn,NdofVo,NdofVn,NvarV
      VTYPE, dimension(NvarV,NdofVo) :: ZdofVo
      VTYPE, dimension(NvarV,NdofVn) :: ZdofVn
c
      integer :: nord1o,nord2o,nord3o,
     .           nord1n,nord2n,nord3n,
     .           ibego,ibegn
c
      if (iprint.eq.1) then
        write(*,7010) Ntype,Nordo,Nordn,NdofVo,NdofVn
 7010   format(' copy_dofV: Ntype,Nordo,Nordn,NdofVo,NdofVn = ',
     .                      a4,2i4,2i6)
      endif
c
      select case(Ntype)
c
c  ...triangle
      case('mdlt')
c
c  .....a hierarchical family
        call copy_1array(NvarV,ZdofVo,NdofVo, ZdofVn,NdofVn)
c
c  ...quadrilateral
      case('mdlq')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
c
c  .....tensor product of hierarchical 1D families
        call copy_2array(NvarV,ZdofVo,nord1o,nord2o, 
     .                           ZdofVn,nord1n,nord2n)
c
c  ...tetrahedron
      case('mdln')
c
c  .....3 intertwined hierarchical families form a hierarchical family
        call copy_1array(NvarV,ZdofVo,NdofVo, ZdofVn,NdofVn)
c
c  ...prism
      case('mdlp')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
        ibego=1; ibegn=1
c
c  .....first 2 families are (horizontally) intertwined 
        call copy_2array(NvarV,
     .                   ZdofVo(1,ibego),(nord1o-1)*nord1o,nord2o, 
     .                   ZdofVn(1,ibegn),(nord1n-1)*nord1n,nord2n)
        ibego=ibego+(nord1o-1)*nord1o*nord2o
        ibegn=ibegn+(nord1n-1)*nord1n*nord2n
c
c  .....3rd family
        call copy_2array(NvarV,
     .                   ZdofVo(1,ibego),nord1o*(nord1o+1)/2,nord2o-1, 
     .                   ZdofVn(1,ibegn),nord1n*(nord1n+1)/2,nord2n-1)
        ibego=ibego+nord1o*(nord1o+1)/2*(nord2o-1) - 1
        ibegn=ibegn+nord1n*(nord1n+1)/2*(nord2n-1) - 1
c
        if ((ibegn.ne.NdofVn).or.(ibego.ne.NdofVo)) then
          write(*,*) 'ibegn,NdofVn,ibego,NdofVo = ',
     .                 ibegn,NdofVn,ibego,NdofVo
          write(*,*) 'copy_dofV: INCONSISTENCY 1'
          stop 1
        endif
c
c  ...hexahedron
      case('mdlb')
        call ddecode(Nordo, nord1o,nord2o,nord3o)
        call ddecode(Nordn, nord1n,nord2n,nord3n)
        ibego=1; ibegn=1
c
c  .....family 1
        call copy_3array(NvarV,ZdofVo(1,ibego),nord1o,nord2o,nord3o-1,
     .                         ZdofVn(1,ibegn),nord1n,nord2n,nord3n-1)
        ibego=ibego+nord1o*nord2o*(nord3o-1)
        ibegn=ibegn+nord1n*nord2n*(nord3n-1)
c
c  .....family 2
        call copy_3array(NvarV,ZdofVo(1,ibego),nord2o,nord3o,nord1o-1,
     .                         ZdofVn(1,ibegn),nord2n,nord3n,nord1n-1)
        ibego=ibego+nord2o*nord3o*(nord1o-1)
        ibegn=ibegn+nord2n*nord3n*(nord1n-1)
c
c  .....family 3
        call copy_3array(NvarV,ZdofVo(1,ibego),nord3o,nord1o,nord2o-1,
     .                         ZdofVn(1,ibegn),nord3n,nord1n,nord2n-1)
        ibego=ibego+nord3o*nord1o*(nord2o-1) - 1
        ibegn=ibegn+nord3n*nord1n*(nord2n-1) - 1
c
        if ((ibegn.ne.NdofVn).or.(ibego.ne.NdofVo)) then
          write(*,*) 'ibegn,NdofVn,ibego,NdofVo = ',
     .                ibegn,NdofVn,ibego,NdofVo
          write(*,*) 'copy_dofV: INCONSISTENCY 2'
          stop 1
        endif
c
c  ...pyramid
      case('mdld')
        ibego=1; ibegn=1
c
c  .....FAMILY 1 AND 2 (curl of families 2 and 3 from H(curl))
        call copy_3array(NvarV,ZdofVo(1,ibego),Nordo,Nordo-1,Nordo-1,
     .                         ZdofVn(1,ibegn),Nordn,Nordn-1,Nordn-1)
        ibego=ibego+Nordo*(Nordo-1)**2
        ibegn=ibegn+Nordn*(Nordn-1)**2
        call copy_3array(NvarV,ZdofVo(1,ibego),Nordo,Nordo-1,Nordo-1,
     .                         ZdofVn(1,ibegn),Nordn,Nordn-1,Nordn-1)
        ibego=ibego+Nordo*(Nordo-1)**2
        ibegn=ibegn+Nordn*(Nordn-1)**2
c
c  .....FAMILY 3 (curl of family 4 from H(curl))
        call copy_2array(NvarV,ZdofVo(1,ibego),Nordo-1,Nordo-1, 
     .                         ZdofVn(1,ibegn),Nordn-1,Nordn-1)
        ibego=ibego+(Nordo-1)**2
        ibegn=ibegn+(Nordn-1)**2
c
c  .....FAMILY 4 (induced from quad face functions)
        call copy_3array(NvarV,ZdofVo(1,ibego),Nordo,Nordo,Nordo-1,
     .                         ZdofVn(1,ibegn),Nordn,Nordn,Nordn-1)
        ibego=ibego+Nordo**2*(Nordo-1)
        ibegn=ibegn+Nordn**2*(Nordn-1)
c
c  .....FAMILY 5
        call copy_2array(NvarV,ZdofVo(1,ibego),Nordo-1,Nordo-1, 
     .                         ZdofVn(1,ibegn),Nordn-1,Nordn-1)
        ibego=ibego+(Nordo-1)**2 
        ibegn=ibegn+(Nordn-1)**2
c
c  .....FAMILY 6 AND 7
        call copy_1array(NvarV,ZdofVo(1,ibego),Nordo-1, 
     .                         ZdofVn(1,ibegn),Nordn-1)
        ibego=ibego+Nordo-1
        ibegn=ibegn+Nordn-1
        call copy_1array(NvarV,ZdofVo(1,ibego),Nordo-1, 
     .                         ZdofVn(1,ibegn),Nordn-1)
        ibego=ibego+Nordo-1 - 1
        ibegn=ibegn+Nordn-1 - 1
c
        if ((ibegn.ne.NdofVn).or.(ibego.ne.NdofVo)) then
          write(*,*) 'ibegn,NdofVn,ibego,NdofVo = ',
     .                ibegn,NdofVn,ibego,NdofVo
          write(*,*) 'copy_dofV: INCONSISTENCY 3'
          stop 1
        endif
c
      end select
      if (iprint.eq.1) write(*,*) 'copy_dofV: DONE'
c
c    
      end subroutine copy_dofV
c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofQ
c
c---------------------------------------------------------------------
c
c   latest revision    - May 20
c
c   purpose            - routine copies L2 dof for a node of order
c                        Nordo to a node of order Nordn
c
c   arguments :
c
c     in:
c        Ntype         - node type
c        Nordo,Nordn   - old and new orders of the node
c        ZdofQo,ZdofQn - old and new dof
c        NvarQ         - first dimension of ZdofQo,ZdofQn
c        NdofQo,NdofQn - old and new # dof
c
c-----------------------------------------------------------------------
c
      subroutine copy_dofQ(Ntype,Nordo,Nordn,NdofQo,NdofQn,
     .                     NvarQ,ZdofQo,ZdofQn)
c
      use parameters, only: MAXEQNQ
      implicit none
      common /ccopy_dofQ/ iprint
      integer :: iprint
c      
      character(len=4)                :: Ntype
      integer                         :: Nordo,Nordn,NdofQo,NdofQn,NvarQ
      VTYPE, dimension (NvarQ,NdofQo) :: ZdofQo
      VTYPE, dimension (NvarQ,NdofQn) :: ZdofQn
c
      integer :: nord1o,nord2o,nord3o,
     .           nord1n,nord2n,nord3n
c
      if (iprint.eq.1) then
        write(*,7010) Ntype,Nordo,Nordn,NdofQo,NdofQn
 7010   format(' copy_dofQ: Ntype,Nordo,Nordn,NdofQo,NdofQn = ',
     .                      a4,2i4,2i6)
      endif
c
      select case(Ntype)
c
c  ...tet
      case('mdln')
c
c  .....a hierarchical family
        call copy_1array(NvarQ,ZdofQo,NdofQo, ZdofQn,NdofQn)
c
c  ...prism
      case('mdlp')
        call decode(Nordo, nord1o,nord2o)
        call decode(Nordn, nord1n,nord2n)
        call copy_2array(NvarQ,ZdofQo,nord1o*(nord1o+1)/2,nord2o, 
     .                         ZdofQn,nord1n*(nord1n+1)/2,nord2n)
c
c  ...hexa
      case('mdlb')
c
c  .....tensor products of 1D hierarchical families
        call ddecode(Nordo, nord1o,nord2o,nord3o)
        call ddecode(Nordn, nord1n,nord2n,nord3n)
        call copy_3array(NvarQ,ZdofQo,nord1o,nord2o,nord3o,
     .                         ZdofQn,nord1n,nord2n,nord3n)
c
c  ...pyramid
      case('mdld')
c
c  .....tensor products of 1D hierarchical families
        call copy_3array(NvarQ,ZdofQo,Nordo,Nordo,Nordo,
     .                           ZdofQn,Nordn,Nordn,Nordn)
c
      end select
      if (iprint.eq.1) write(*,*) 'copy_dofQ: DONE'
c
c
      end subroutine copy_dofQ
c
c
      subroutine copy_1array(M,A,Ia, B,Ib)
c
      use parameters, only: ZERO
      implicit none
      integer                 :: M,Ia,Ib,i
      VTYPE, dimension(M,Ia) :: A
      VTYPE, dimension(M,Ib) :: B
c
      B = ZERO
      i=min(Ia,Ib)
      B(1:M,1:i) =  A(1:M,1:i)
c
      end subroutine copy_1array
c
c
      subroutine copy_2array(M,A,Ia,Ja, B,Ib,Jb)
c
      use parameters, only: ZERO
      implicit none
      integer                 :: M,Ia,Ja,Ib,Jb,i,j
      VTYPE, dimension(M,Ia,Ja) :: A
      VTYPE, dimension(M,Ib,Jb) :: B
c
      B = ZERO
      i=min(Ia,Ib)
      j=min(Ja,Jb)
      B(1:M,1:i,1:j) =  A(1:M,1:i,1:j)
c
      end subroutine copy_2array
c
c
      subroutine copy_3array(M,A,Ia,Ja,Ka, B,Ib,Jb,Kb)
c
      use parameters, only: ZERO
      implicit none
      integer                 :: M,Ia,Ja,Ka,Ib,Jb,Kb,i,j,k
      VTYPE, dimension(M,Ia,Ja,Ka) :: A
      VTYPE, dimension(M,Ib,Jb,Kb) :: B
c
      B = ZERO
      i=min(Ia,Ib)
      j=min(Ja,Jb)
      k=min(Ka,Kb)
      B(1:M,1:i,1:j,1:k) =  A(1:M,1:i,1:j,1:k)
c
      end subroutine copy_3array
c










