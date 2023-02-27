!
#include "typedefs.h"
!
!---------------------------------------------------------------------
!
!   routine name       - nodmod
!
!---------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine modifies the order of approximation
!                        for a node
!
!   arguments:
!
!     in:
!             Nod      - node number
!             Nordn    - new order of approximation
!
!-----------------------------------------------------------------------
!
      subroutine nodmod(Nod,Nordn)
!
      use data_structure3D
      use mpi_param, only: RANK
      use par_mesh , only: DISTRIBUTED
!
      implicit none
!
#if DEBUG_MODE
      common /ccopy_dofG/ iprint_copy_dofG
      common /ccopy_dofH/ iprint_copy_dofH
      common /ccopy_dofE/ iprint_copy_dofE
      common /ccopy_dofV/ iprint_copy_dofV
      common /ccopy_dofQ/ iprint_copy_dofQ
      integer :: iprint_copy_dofG,iprint_copy_dofH,iprint_copy_dofE, &
                 iprint_copy_dofV,iprint_copy_dofQ
#endif
!
      integer, intent(in)   :: Nod,Nordn
!
      real(8), allocatable ::  xnod(:,:)
      VTYPE  , allocatable :: zdofH(:,:)
      VTYPE  , allocatable :: zdofE(:,:)
      VTYPE  , allocatable :: zdofV(:,:)
      VTYPE  , allocatable :: zdofQ(:,:)
!
!  ...set whether current dofxs should be copied to modified node
      logical, parameter :: COPY_DOFS = .false.
!
      logical :: act_dof
      integer :: ntype,iprint,icase,nordo,subd, &
                 nvarH ,nvarE ,nvarV ,nvarQ,    &
                 ndofH ,ndofE ,ndofV ,ndofQ,    &
                 ndofHo,ndofEo,ndofVo,ndofQo
!
#if DEBUG_MODE
      iprint=0
      iprint_copy_dofG=iprint
      iprint_copy_dofH=iprint; iprint_copy_dofE=iprint
      iprint_copy_dofV=iprint; iprint_copy_dofQ=iprint
#endif
!
!  ...node case
      icase = NODES(Nod)%case
      ntype = NODES(Nod)%ntype
      nordo = NODES(Nod)%order
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7010) Nod,S_Type(ntype),nordo,icase
 7010    format('nodmod: Nod,type,order,icase = ',i6,2x,a5,2x,2i3)
         write(*,7013) NREQNH(icase),NREQNE(icase),NREQNV(icase),NREQNQ(icase)
 7013    format('nodmod: NREQNH,NREQNE,NREQNV,NREQNQ(icase) = ',4i2)
      endif
#endif
!
!  ...exit if the order is the same
      if (Nordn.eq.nordo) return
!
!  ...determine the current number of dof
      call ndof_nod(ntype,nordo, ndofH,ndofE,ndofV,ndofQ)
!
!  ...determine the number of components supported on this node
      nvarH = NREQNH(icase)*NRCOMS
      nvarE = NREQNE(icase)*NRCOMS
      nvarV = NREQNV(icase)*NRCOMS
      nvarQ = NREQNQ(icase)*NRCOMS
!
!  ...save current gdofs in local array and update number of gdofs
      if (ndofH.gt.0) then
         if (associated(NODES(Nod)%dof)) then
            if (associated(NODES(Nod)%dof%coord)) then
               if (COPY_DOFS) then
                  allocate(xnod(NDIMEN,ndofH))
                  xnod(1:NDIMEN,1:ndofH) = NODES(Nod)%dof%coord(1:NDIMEN,1:ndofH)
               endif
               deallocate(NODES(Nod)%dof%coord)
            endif
         endif
      endif
!
!  ...save current dof in local arrays
      if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
         if (associated(NODES(Nod)%dof)) then
            if (associated(NODES(Nod)%dof%zdofH)) then
               if (COPY_DOFS) then
                  allocate(zdofH(nvarH,ndofH))
                  zdofH(1:nvarH,1:ndofH) = NODES(Nod)%dof%zdofH(1:nvarH,1:ndofH)
               endif
               deallocate(NODES(Nod)%dof%zdofH)
            endif
         endif
         NRDOFSH = NRDOFSH - ndofH*NREQNH(icase)
      endif
      if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
         if (associated(NODES(Nod)%dof)) then
            if (associated(NODES(Nod)%dof%zdofE)) then
               if (COPY_DOFS) then
                  allocate(zdofE(nvarE,ndofE))
                  zdofE(1:nvarE,1:ndofE) = NODES(Nod)%dof%zdofE(1:nvarE,1:ndofE)
               endif
               deallocate(NODES(Nod)%dof%zdofE)
            endif
         endif
         NRDOFSE = NRDOFSE - ndofE*NREQNE(icase)
      endif
      if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
         if (associated(NODES(Nod)%dof)) then
            if (associated(NODES(Nod)%dof%zdofV)) then
               if (COPY_DOFS) then
                  allocate(zdofV(nvarV,ndofV))
                  zdofV(1:nvarV,1:ndofV) = NODES(Nod)%dof%zdofV(1:nvarV,1:ndofV)
               endif
               deallocate(NODES(Nod)%dof%zdofV)
            endif
         endif
         NRDOFSV = NRDOFSV - ndofV*NREQNV(icase)
      endif
      if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
         if (associated(NODES(Nod)%dof)) then
            if (associated(NODES(Nod)%dof%zdofQ)) then
               if (COPY_DOFS) then
                  allocate(zdofQ(nvarQ,ndofQ))
                  zdofQ(1:nvarQ,1:ndofQ) = NODES(Nod)%dof%zdofQ(1:nvarQ,1:ndofQ)
               endif
               deallocate(NODES(Nod)%dof%zdofQ)
            endif
         endif
         NRDOFSQ = NRDOFSQ - ndofQ*NREQNQ(icase)
      endif
      if (associated(NODES(Nod)%dof)) then
         deallocate(NODES(Nod)%dof)
      endif
!
!  ...save number of dof corresponding to the old order
      ndofHo = ndofH
      ndofEo = ndofE
      ndofVo = ndofV
      ndofQo = ndofQ
!
!  ...modify the order of approximation
      NODES(Nod)%order = Nordn
!
!  ...calculate the new number of dof for the node
      call ndof_nod(ntype,Nordn, ndofH,ndofE,ndofV,ndofQ)
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7015) ndofH,ndofE,ndofV,ndofQ
 7015    format('nodmod: ndofH,ndofE,ndofV,ndofQ = ',4i5)
      endif
#endif
!
!  ...find out whether node is inside my subdomain
      act_dof = .true.
      if (DISTRIBUTED) then
         call get_subd(Nod, subd)
         if (RANK .ne. subd) act_dof = .false.
      endif
!
      if (act_dof) then
         allocate(NODES(Nod)%dof)
         nullify(NODES(Nod)%dof%coord)
         nullify(NODES(Nod)%dof%zdofH)
         nullify(NODES(Nod)%dof%zdofE)
         nullify(NODES(Nod)%dof%zdofV)
         nullify(NODES(Nod)%dof%zdofQ)
      endif
!
!  ...allocate memory for geometry and solution dof
!     and copy old dof from the local arrays
      if (ndofH.gt.0) then
         if (act_dof) then
            allocate(NODES(Nod)%dof%coord(NDIMEN,ndofH))
            NODES(Nod)%dof%coord = 0.d0
            if (allocated(xnod)) then
               call copy_dofG(ntype,nordo,Nordn,ndofHo,ndofH, &
                              xnod,NODES(Nod)%dof%coord)
            endif
         endif
      endif
      if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
         if (act_dof) then
            allocate(NODES(Nod)%dof%zdofH(nvarH,ndofH))
            NODES(Nod)%dof%zdofH = ZERO
            if (allocated(zdofH)) then
               call copy_dofH(ntype,nordo,Nordn,ndofHo,ndofH, &
                              nvarH,zdofH,NODES(Nod)%dof%zdofH)
            endif
         endif
         NRDOFSH = NRDOFSH + ndofH*NREQNH(icase)
      endif
      if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
         if (act_dof) then
            allocate(NODES(Nod)%dof%zdofE(nvarE,ndofE))
            NODES(Nod)%dof%zdofE = ZERO
            if (allocated(zdofE)) then
               call copy_dofE(ntype,nordo,Nordn,ndofEo,ndofE, &
                              nvarE,zdofE,NODES(Nod)%dof%zdofE)
            endif
         endif
         NRDOFSE = NRDOFSE + ndofE*NREQNE(icase)
      endif
      if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
         if (act_dof) then
            allocate(NODES(Nod)%dof%zdofV(nvarV,ndofV))
            NODES(Nod)%dof%zdofV = ZERO
            if (allocated(zdofV)) then
               call copy_dofV(ntype,nordo,Nordn,ndofVo,ndofV, &
                              nvarV,zdofV,NODES(Nod)%dof%zdofV)
            endif
         endif
         NRDOFSV = NRDOFSV + ndofV*NREQNV(icase)
      endif
      if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
         if (act_dof) then
            allocate( NODES(Nod)%dof%zdofQ(nvarQ,ndofQ))
            NODES(Nod)%dof%zdofQ = ZERO
            if (allocated(zdofQ)) then
               call copy_dofQ(ntype,nordo,Nordn,ndofQo,ndofQ, &
                              nvarQ,zdofQ,NODES(Nod)%dof%zdofQ)
            endif
         endif
         NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(icase)
      endif
!
      if (allocated(xnod )) deallocate(xnod )
      if (allocated(zdofH)) deallocate(zdofH)
      if (allocated(zdofE)) deallocate(zdofE)
      if (allocated(zdofV)) deallocate(zdofV)
      if (allocated(zdofQ)) deallocate(zdofQ)
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7020) Nod
 7020    format('nodmod: Nod = ',i6,' HAS BEEN UPDATED ')
         write(*,7030) nordo,Nordn
 7030    format('        OLD ORDER = ',i3,', NEW ORDER = ',i3)
      endif
#endif
!
!
      end subroutine nodmod
!
!---------------------------------------------------------------------
!
!   routine name       - copy_dofG
!
!---------------------------------------------------------------------
!
!   latest revision    - June 2020
!
!   purpose            - routine copies G1 dof for a node of order
!                        Nordo to a node of order Nordn
!
!   arguments:
!
!     in:
!        Ntype         - node type
!        Nordo,Nordn   - old and new orders of the node
!        Xnodo,Xnodn   - old and new dof
!        NdofGo,NdofGn - old and new # dof
!
!-----------------------------------------------------------------------
!
      subroutine copy_dofG(Ntype,Nordo,Nordn,NdofGo,NdofGn, &
                           Xnodo,Xnodn)
!
      use node_types
      use parameters, only: NDIMEN
      implicit none
!
#if DEBUG_MODE
      common /ccopy_dofG/ iprint
      integer :: iprint
#endif
!
      integer :: Ntype,Nordo,Nordn,NdofGo,NdofGn
      real(8) :: Xnodo(NDIMEN,NdofGo)
      real(8) :: Xnodn(NDIMEN,NdofGn)
!
      integer :: nord1o,nord2o,nord3o, &
                 nord1n,nord2n,nord3n
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7010) S_Type(Ntype),Nordo,Nordn,NdofGo,NdofGn
 7010    format(' copy_dofG: Ntype,Nordo,Nordn,NdofGo,NdofGn = ', &
                             a4,2i4,2i6)
      endif
#endif
!
      select case(Ntype)
!
!  ...edge, triangle, tet, pyramid
      case(MEDG,MDLT,MDLN,MDLD)
         call copy_1array_r(NDIMEN,Xnodo,NdofGo, Xnodn,NdofGn)
!
!  ...quad
      case(MDLQ)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
         call copy_2array_r(NDIMEN,Xnodo,nord1o-1,nord2o-1, &
                                   Xnodn,nord1n-1,nord2n-1)
!
!  ...prism
      case(MDLP)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
         call copy_2array_r(NDIMEN,Xnodo,(nord1o-2)*(nord1o-1)/2,nord2o-1, &
                                   Xnodn,(nord1n-2)*(nord1n-1)/2,nord2n-1)
!
!  ...hexa
      case(MDLB)
         call ddecode(Nordo, nord1o,nord2o,nord3o)
         call ddecode(Nordn, nord1n,nord2n,nord3n)
         call copy_3array_r(NDIMEN,Xnodo,nord1o-1,nord2o-1,nord3o-1, &
                                   Xnodn,nord1n-1,nord2n-1,nord3n-1)
!
      end select
!
#if DEBUG_MODE
      if (iprint.eq.1) write(*,*) 'copy_dofG: DONE'
#endif
!
!
      end subroutine copy_dofG
!
!---------------------------------------------------------------------
!
!   routine name       - copy_dofH
!
!---------------------------------------------------------------------
!
!   latest revision    - June 2020
!
!   purpose            - routine copies H1 dof for a node of order
!                        Nordo to a node of order Nordn
!
!   arguments :
!
!     in:
!        Ntype         - node type
!        Nordo,Nordn   - old and new orders of the node
!        NdofHo,NdofHn - old and new # dof
!        NvarH         - first dimension of ZdofHo,ZdofHn
!        ZdofHo,ZdofHn - old and new dof
!
!-----------------------------------------------------------------------
!
      subroutine copy_dofH(Ntype,Nordo,Nordn,NdofHo,NdofHn, &
                           NvarH,ZdofHo,ZdofHn)
!
      use node_types
      implicit none
!
#if DEBUG_MODE
      common /ccopy_dofH/ iprint
      integer :: iprint
#endif
!
      integer :: Ntype,Nordo,Nordn,NdofHo,NdofHn,NvarH
      VTYPE   :: ZdofHo(NvarH,NdofHo)
      VTYPE   :: ZdofHn(NvarH,NdofHn)
!
      integer :: nord1o,nord2o,nord3o, &
                 nord1n,nord2n,nord3n
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7010) S_Type(Ntype),Nordo,Nordn,NdofHo,NdofHn
 7010    format(' copy_dofH: Ntype,Nordo,Nordn,NdofHo,NdofHn = ', &
                             a4,2i4,2i6)
      endif
#endif
!
      select case(Ntype)
!
!  ...edge, triangle, tet, pyramid
      case(MEDG,MDLT,MDLN,MDLD)
         call copy_1array(NvarH,ZdofHo,NdofHo, ZdofHn,NdofHn)
!
!  ...quad
      case(MDLQ)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
         call copy_2array(NvarH,ZdofHo,nord1o-1,nord2o-1, &
                                ZdofHn,nord1n-1,nord2n-1)
!
!  ...prism
      case(MDLP)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
         call copy_2array(NvarH,                                   &
                          ZdofHo,(nord1o-2)*(nord1o-1)/2,nord2o-1, &
                          ZdofHn,(nord1n-2)*(nord1n-1)/2,nord2n-1)
!
!  ...hexa
      case(MDLB)
         call ddecode(Nordo, nord1o,nord2o,nord3o)
         call ddecode(Nordn, nord1n,nord2n,nord3n)
         call copy_3array(NvarH,ZdofHo,nord1o-1,nord2o-1,nord3o-1, &
                                ZdofHn,nord1n-1,nord2n-1,nord3n-1)
!
      end select
!
#if DEBUG_MODE
      if (iprint.eq.1) write(*,*) 'copy_dofH: DONE'
#endif
!
!
      end subroutine copy_dofH
!
!---------------------------------------------------------------------
!
!   routine name       - copy_dofE
!
!---------------------------------------------------------------------
!
!   latest revision    - June 2020
!
!   purpose            - routine copies H(curl) dof for a node of order
!                        Nordo to a node of order Nordn
!
!   arguments:
!
!     in:
!        Ntype         - node type
!        Nordo,Nordn   - old and new orders of the node
!        NdofEo,NdofEn - old and new # dof
!        NvarE         - first dimension of ZdofEo,ZdofEn
!        ZdofEo,ZdofEn - old and new dof
!
!-----------------------------------------------------------------------
!
      subroutine copy_dofE(Ntype,Nordo,Nordn,NdofEo,NdofEn, &
                           NvarE,ZdofEo,ZdofEn)
!
      use node_types
      implicit none
!
#if DEBUG_MODE
      common /ccopy_dofE/ iprint
      integer :: iprint
#endif
!
      integer :: Ntype
      integer :: Nordo,Nordn,NdofEo,NdofEn,NvarE
      VTYPE   :: ZdofEo(NvarE,NdofEo)
      VTYPE   :: ZdofEn(NvarE,NdofEn)
!
      integer :: nord1o,nord2o,nord3o, &
                 nord1n,nord2n,nord3n, &
                 ibego,ibegn
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7010) Ntype,Nordo,Nordn,NdofEo,NdofEn
 7010    format(' copy_dofE: Ntype,Nordo,Nordn,NdofEo,NdofEn = ', &
                             a4,2i4,2i6)
      endif
#endif
!
      select case(Ntype)
!
!  ...segment
      case(MEDG)
         call copy_1array(NvarE,ZdofEo,NdofEo, ZdofEn,NdofEn)
!
!  ...triangle
      case(MDLT)
!
!     ...the 2 families are intertwined to form a hierarchy
         call copy_1array(NvarE,ZdofEo,NdofEo, ZdofEn,NdofEn)
!
!  ...quadrilateral
      case(MDLQ)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
         ibego=1; ibegn=1
!
!     ...family 1
         call copy_2array(NvarE,ZdofEo(1,ibego),nord1o,nord2o-1, &
                                ZdofEn(1,ibegn),nord1n,nord2n-1)
         ibego=ibego+nord1o*(nord2o-1)
         ibegn=ibegn+nord1n*(nord2n-1)
!
!     ...family 2
         call copy_2array(NvarE,ZdofEo(1,ibego),nord1o-1,nord2o, &
                                ZdofEn(1,ibegn),nord1n-1,nord2n)
         ibego=ibego+(nord1o-1)*nord2o - 1
         ibegn=ibegn+(nord1n-1)*nord2n - 1
!
         if ((ibegn.ne.NdofEn).or.(ibego.ne.NdofEo)) then
            write(*,*) 'copy_dofE: INCONSISTENCY 1'
            stop
         endif
!
!  ...tetrahedron
      case(MDLN)
!
!     ...the 3 families are intertwined to form a hierarchy
         call copy_1array(NvarE,ZdofEo,NdofEo, ZdofEn,NdofEn)
!
!  ...prism
      case(MDLP)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
         ibego=1; ibegn=1
!
!     ...first 2 families are (horizontally) intertwined
         call copy_2array(NvarE,                                      &
                          ZdofEo(1,ibego),nord1o*(nord1o-1),nord2o-1, &
                          ZdofEn(1,ibegn),nord1n*(nord1n-1),nord2n-1)
         ibego=ibego+nord1o*(nord1o-1)*(nord2o-1)
         ibegn=ibegn+nord1n*(nord1n-1)*(nord2n-1)
!
!     ...3rd family
         call copy_2array(NvarE,                                          &
                          ZdofEo(1,ibego),(nord1o-2)*(nord1o-1)/2,nord2o, &
                          ZdofEn(1,ibegn),(nord1n-2)*(nord1n-1)/2,nord2n)
!-----------------------------------------------------------------------
         ibego=ibego+(nord1o-2)*(nord1o-1)/2*nord2o - 1
         ibegn=ibegn+(nord1n-2)*(nord1n-1)/2*nord2n - 1
!
         if ((ibegn.ne.NdofEn).or.(ibego.ne.NdofEo)) then
            write(*,*) 'ibegn,NdofEn,ibego,NdofEo = ', &
                        ibegn,NdofEn,ibego,NdofEo
            write(*,*) 'copy_dofE: INCONSISTENCY 2'
            stop
         endif
!
!  ...hexahedron
      case(MDLB)
         call ddecode(Nordo, nord1o,nord2o,nord3o)
         call ddecode(Nordn, nord1n,nord2n,nord3n)
         ibego=1; ibegn=1
!
!     ...family 1
         call copy_3array(NvarE,                                    &
                          ZdofEo(1,ibego),nord1o,nord2o-1,nord3o-1, &
                          ZdofEn(1,ibegn),nord1n,nord2n-1,nord3n-1)
         ibego=ibego+nord1o*(nord2o-1)*(nord3o-1)
         ibegn=ibegn+nord1n*(nord2n-1)*(nord3n-1)
!
!     ...family 2
         call copy_3array(NvarE,                                    &
                          ZdofEo(1,ibego),nord1o-1,nord2o,nord3o-1, &
                          ZdofEn(1,ibegn),nord1n-1,nord2n,nord3n-1)
         ibego=ibego+(nord1o-1)*nord2o*(nord3o-1)
         ibegn=ibegn+(nord1n-1)*nord2n*(nord3n-1)
!
!     ...family 3
         call copy_3array(NvarE,                                    &
                          ZdofEo(1,ibego),nord1o-1,nord2o-1,nord3o, &
                          ZdofEn(1,ibegn),nord1n-1,nord2n-1,nord3n)
         ibego=ibego+(nord1o-1)*(nord2o-1)*nord3o - 1
         ibegn=ibegn+(nord1n-1)*(nord2n-1)*nord3n - 1
!
         if ((ibegn.ne.NdofEn).or.(ibego.ne.NdofEo)) then
            write(*,*) 'ibegn,NdofEn,ibego,NdofEo = ', &
                        ibegn,NdofEn,ibego,NdofEo
            write(*,*) 'copy_dofE: INCONSISTENCY 3'
            stop
         endif
!
!     ...pyramid
         case(MDLD)
         ibego=1; ibegn=1
!
!     ...FAMILY 1 (gradients of H1 bubbles)
         call copy_3array(NvarE,                                   &
                          ZdofEo(1,ibego),Nordo-1,Nordo-1,Nordo-1, &
                          ZdofEn(1,ibegn),Nordn-1,Nordn-1,Nordn-1)
         ibego=ibego+(Nordo-1)**3
         ibegn=ibegn+(Nordn-1)**3
!
!     ...FAMILY 2 AND 3 (induced from quad face functions)
         call copy_3array(NvarE,ZdofEo(1,ibego),Nordo,Nordo-1,Nordo-1, &
                                ZdofEn(1,ibegn),Nordn,Nordn-1,Nordn-1)
         ibego=ibego+Nordo*(Nordo-1)**2
         ibegn=ibegn+Nordn*(Nordn-1)**2
         call copy_3array(NvarE,ZdofEo(1,ibego),Nordo-1,Nordo,Nordo-1, &
                                ZdofEn(1,ibegn),Nordn-1,Nordn,Nordn-1)
         ibego=ibego+Nordo*(Nordo-1)**2
         ibegn=ibegn+Nordn*(Nordn-1)**2
!
!     ...FAMILY 4
         call copy_2array(NvarE,                           &
                          ZdofEo(1,ibego),Nordo-1,Nordo-1, &
                          ZdofEn(1,ibegn),Nordn-1,Nordn-1)
         ibego=ibego+(Nordo-1)**2 - 1
         ibegn=ibegn+(Nordn-1)**2 - 1
!
         if ((ibegn.ne.NdofEn).or.(ibego.ne.NdofEo)) then
            write(*,*) 'copy_dofE: INCONSISTENCY 4'
            stop
         endif
!
      end select
!
#if DEBUG_MODE
      if (iprint.eq.1) write(*,*) 'copy_dofE: DONE'
#endif
!
!
      end subroutine copy_dofE
!
!---------------------------------------------------------------------
!
!   routine name       - copy_dofV
!
!---------------------------------------------------------------------
!
!   latest revision    - June 2020
!
!   purpose            - routine copies H(div) dof for a node of order
!                        Nordo to a node of order Nordn
!
!   arguments:
!
!     in:
!        Ntype         - node type
!        Nordo,Nordn   - old and new orders of the node
!        NdofEo,NdofEn - old and new # dof
!        NvarV         - first dimension of ZdofVo,ZdofVn
!        ZdofEo,ZdofEn - old and new dof
!
!-----------------------------------------------------------------------
!
      subroutine copy_dofV(Ntype,Nordo,Nordn,NdofVo,NdofVn, &
                           NvarV,ZdofVo,ZdofVn)
!
      use node_types
      implicit none
!
#if DEBUG_MODE
      common /ccopy_dofV/ iprint
      integer :: iprint
#endif
!
      integer :: Ntype,Nordo,Nordn,NdofVo,NdofVn,NvarV
      VTYPE   :: ZdofVo(NvarV,NdofVo)
      VTYPE   :: ZdofVn(NvarV,NdofVn)
!
      integer :: nord1o,nord2o,nord3o, &
                 nord1n,nord2n,nord3n, &
                 ibego,ibegn
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7010) S_Type(Ntype),Nordo,Nordn,NdofVo,NdofVn
 7010    format(' copy_dofV: Ntype,Nordo,Nordn,NdofVo,NdofVn = ', &
                             a4,2i4,2i6)
      endif
#endif
!
      select case(Ntype)
!
!  ...triangle
      case(MDLT)
!
!     ...a hierarchical family
         call copy_1array(NvarV,ZdofVo,NdofVo, ZdofVn,NdofVn)
!
!  ...quadrilateral
      case(MDLQ)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
!
!     ...tensor product of hierarchical 1D families
         call copy_2array(NvarV,ZdofVo,nord1o,nord2o, &
                                ZdofVn,nord1n,nord2n)
!
!  ...tetrahedron
      case(MDLN)
!
!     ...3 intertwined hierarchical families form a hierarchical family
         call copy_1array(NvarV,ZdofVo,NdofVo, ZdofVn,NdofVn)
!
!  ...prism
      case(MDLP)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
         ibego=1; ibegn=1
!
!     ...first 2 families are (horizontally) intertwined
         call copy_2array(NvarV,                                    &
                          ZdofVo(1,ibego),(nord1o-1)*nord1o,nord2o, &
                          ZdofVn(1,ibegn),(nord1n-1)*nord1n,nord2n)
         ibego=ibego+(nord1o-1)*nord1o*nord2o
         ibegn=ibegn+(nord1n-1)*nord1n*nord2n
!
!     ...3rd family
         call copy_2array(NvarV,                                        &
                          ZdofVo(1,ibego),nord1o*(nord1o+1)/2,nord2o-1, &
                          ZdofVn(1,ibegn),nord1n*(nord1n+1)/2,nord2n-1)
         ibego=ibego+nord1o*(nord1o+1)/2*(nord2o-1) - 1
         ibegn=ibegn+nord1n*(nord1n+1)/2*(nord2n-1) - 1
!
         if ((ibegn.ne.NdofVn).or.(ibego.ne.NdofVo)) then
            write(*,*) 'ibegn,NdofVn,ibego,NdofVo = ', &
                        ibegn,NdofVn,ibego,NdofVo
            write(*,*) 'copy_dofV: INCONSISTENCY 1'
            stop
         endif
!
!  ...hexahedron
      case(MDLB)
         call ddecode(Nordo, nord1o,nord2o,nord3o)
         call ddecode(Nordn, nord1n,nord2n,nord3n)
         ibego=1; ibegn=1
!
!     ...family 1
         call copy_3array(NvarV,ZdofVo(1,ibego),nord1o,nord2o,nord3o-1, &
                                ZdofVn(1,ibegn),nord1n,nord2n,nord3n-1)
         ibego=ibego+nord1o*nord2o*(nord3o-1)
         ibegn=ibegn+nord1n*nord2n*(nord3n-1)
!
!     ...family 2
         call copy_3array(NvarV,ZdofVo(1,ibego),nord2o,nord3o,nord1o-1, &
                                ZdofVn(1,ibegn),nord2n,nord3n,nord1n-1)
         ibego=ibego+nord2o*nord3o*(nord1o-1)
         ibegn=ibegn+nord2n*nord3n*(nord1n-1)
!
!     ...family 3
         call copy_3array(NvarV,ZdofVo(1,ibego),nord3o,nord1o,nord2o-1, &
                                ZdofVn(1,ibegn),nord3n,nord1n,nord2n-1)
         ibego=ibego+nord3o*nord1o*(nord2o-1) - 1
         ibegn=ibegn+nord3n*nord1n*(nord2n-1) - 1
!
         if ((ibegn.ne.NdofVn).or.(ibego.ne.NdofVo)) then
            write(*,*) 'ibegn,NdofVn,ibego,NdofVo = ', &
                        ibegn,NdofVn,ibego,NdofVo
            write(*,*) 'copy_dofV: INCONSISTENCY 2'
            stop
         endif
!
!  ...pyramid
      case(MDLD)
         ibego=1; ibegn=1
!
!     ...FAMILY 1 AND 2 (curl of families 2 and 3 from H(curl))
         call copy_3array(NvarV,ZdofVo(1,ibego),Nordo,Nordo-1,Nordo-1, &
                                ZdofVn(1,ibegn),Nordn,Nordn-1,Nordn-1)
         ibego=ibego+Nordo*(Nordo-1)**2
         ibegn=ibegn+Nordn*(Nordn-1)**2
         call copy_3array(NvarV,ZdofVo(1,ibego),Nordo,Nordo-1,Nordo-1, &
                                ZdofVn(1,ibegn),Nordn,Nordn-1,Nordn-1)
         ibego=ibego+Nordo*(Nordo-1)**2
         ibegn=ibegn+Nordn*(Nordn-1)**2
!
!     ...FAMILY 3 (curl of family 4 from H(curl))
         call copy_2array(NvarV,ZdofVo(1,ibego),Nordo-1,Nordo-1, &
                                ZdofVn(1,ibegn),Nordn-1,Nordn-1)
         ibego=ibego+(Nordo-1)**2
         ibegn=ibegn+(Nordn-1)**2
!
!     ...FAMILY 4 (induced from quad face functions)
         call copy_3array(NvarV,ZdofVo(1,ibego),Nordo,Nordo,Nordo-1, &
                                ZdofVn(1,ibegn),Nordn,Nordn,Nordn-1)
         ibego=ibego+Nordo**2*(Nordo-1)
         ibegn=ibegn+Nordn**2*(Nordn-1)
!
!     ...FAMILY 5
         call copy_2array(NvarV,ZdofVo(1,ibego),Nordo-1,Nordo-1, &
                                ZdofVn(1,ibegn),Nordn-1,Nordn-1)
         ibego=ibego+(Nordo-1)**2
         ibegn=ibegn+(Nordn-1)**2
!
!     ...FAMILY 6 AND 7
         call copy_1array(NvarV,ZdofVo(1,ibego),Nordo-1, &
                                ZdofVn(1,ibegn),Nordn-1)
         ibego=ibego+Nordo-1
         ibegn=ibegn+Nordn-1
         call copy_1array(NvarV,ZdofVo(1,ibego),Nordo-1, &
                                ZdofVn(1,ibegn),Nordn-1)
         ibego=ibego+Nordo-1 - 1
         ibegn=ibegn+Nordn-1 - 1
!
         if ((ibegn.ne.NdofVn).or.(ibego.ne.NdofVo)) then
            write(*,*) 'ibegn,NdofVn,ibego,NdofVo = ', &
                        ibegn,NdofVn,ibego,NdofVo
            write(*,*) 'copy_dofV: INCONSISTENCY 3'
            stop
         endif
!
      end select
!
#if DEBUG_MODE
      if (iprint.eq.1) write(*,*) 'copy_dofV: DONE'
#endif
!
!
      end subroutine copy_dofV
!
!---------------------------------------------------------------------
!
!   routine name       - copy_dofQ
!
!---------------------------------------------------------------------
!
!   latest revision    - June 2020
!
!   purpose            - routine copies L2 dof for a node of order
!                        Nordo to a node of order Nordn
!
!   arguments:
!
!     in:
!        Ntype         - node type
!        Nordo,Nordn   - old and new orders of the node
!        ZdofQo,ZdofQn - old and new dof
!        NvarQ         - first dimension of ZdofQo,ZdofQn
!        NdofQo,NdofQn - old and new # dof
!
!-----------------------------------------------------------------------
!
      subroutine copy_dofQ(Ntype,Nordo,Nordn,NdofQo,NdofQn, &
                           NvarQ,ZdofQo,ZdofQn)
!
      use node_types
      implicit none
!
#if DEBUG_MODE
      common /ccopy_dofQ/ iprint
      integer :: iprint
#endif
!
      integer :: Ntype,Nordo,Nordn,NdofQo,NdofQn,NvarQ
      VTYPE   :: ZdofQo(NvarQ,NdofQo)
      VTYPE   :: ZdofQn(NvarQ,NdofQn)
!
      integer :: nord1o,nord2o,nord3o, &
                 nord1n,nord2n,nord3n
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7010) S_Type(Ntype),Nordo,Nordn,NdofQo,NdofQn
 7010    format(' copy_dofQ: Ntype,Nordo,Nordn,NdofQo,NdofQn = ', &
                             a4,2i4,2i6)
      endif
#endif
!
      select case(Ntype)
!
!  ...tet
      case(MDLN)
!
!     ...a hierarchical family
         call copy_1array(NvarQ,ZdofQo,NdofQo, ZdofQn,NdofQn)
!
!  ...prism
      case(MDLP)
         call decode(Nordo, nord1o,nord2o)
         call decode(Nordn, nord1n,nord2n)
         call copy_2array(NvarQ,ZdofQo,nord1o*(nord1o+1)/2,nord2o, &
                                ZdofQn,nord1n*(nord1n+1)/2,nord2n)
!
!  ...hexa
      case(MDLB)
!
!     ...tensor products of 1D hierarchical families
         call ddecode(Nordo, nord1o,nord2o,nord3o)
         call ddecode(Nordn, nord1n,nord2n,nord3n)
         call copy_3array(NvarQ,ZdofQo,nord1o,nord2o,nord3o, &
                                ZdofQn,nord1n,nord2n,nord3n)
!
!  ...pyramid
      case(MDLD)
!
!     ...tensor products of 1D hierarchical families
         call copy_3array(NvarQ,ZdofQo,Nordo,Nordo,Nordo, &
                                ZdofQn,Nordn,Nordn,Nordn)
!
      end select
!
#if DEBUG_MODE
      if (iprint.eq.1) write(*,*) 'copy_dofQ: DONE'
#endif
!
!
      end subroutine copy_dofQ
!
!-----------------------------------------------------------------------
!
      subroutine copy_1array(M,A,Ia, B,Ib)
!
      use parameters, only: ZERO
      implicit none
      integer :: M,Ia,Ib,i
      VTYPE   :: A(M,Ia)
      VTYPE   :: B(M,Ib)
!
      B = ZERO
      i = min(Ia,Ib)
      B(1:M,1:i) = A(1:M,1:i)
!
      end subroutine copy_1array
!
!-----------------------------------------------------------------------
!
      subroutine copy_2array(M,A,Ia,Ja, B,Ib,Jb)
!
      use parameters, only: ZERO
      implicit none
      integer :: M,Ia,Ja,Ib,Jb,i,j
      VTYPE   :: A(M,Ia,Ja)
      VTYPE   :: B(M,Ib,Jb)
!
      B = ZERO
      i = min(Ia,Ib)
      j = min(Ja,Jb)
      B(1:M,1:i,1:j) = A(1:M,1:i,1:j)
!
      end subroutine copy_2array
!
!-----------------------------------------------------------------------
!
      subroutine copy_3array(M,A,Ia,Ja,Ka, B,Ib,Jb,Kb)
!
      use parameters, only: ZERO
      implicit none
      integer :: M,Ia,Ja,Ka,Ib,Jb,Kb,i,j,k
      VTYPE   :: A(M,Ia,Ja,Ka)
      VTYPE   :: B(M,Ib,Jb,Kb)
!
      B = ZERO
      i = min(Ia,Ib)
      j = min(Ja,Jb)
      k = min(Ka,Kb)
      B(1:M,1:i,1:j,1:k) = A(1:M,1:i,1:j,1:k)
!
      end subroutine copy_3array

!-----------------------------------------------------------------------
!
      subroutine copy_1array_r(M,A,Ia, B,Ib)
!
      implicit none
      integer :: M,Ia,Ib,i
      real(8) :: A(M,Ia)
      real(8) :: B(M,Ib)
!
      B = 0.d0
      i = min(Ia,Ib)
      B(1:M,1:i) = A(1:M,1:i)
!
      end subroutine copy_1array_r
!
!-----------------------------------------------------------------------
!
      subroutine copy_2array_r(M,A,Ia,Ja, B,Ib,Jb)
!
      implicit none
      integer :: M,Ia,Ja,Ib,Jb,i,j
      real(8) :: A(M,Ia,Ja)
      real(8) :: B(M,Ib,Jb)
!
      B = 0.d0
      i = min(Ia,Ib)
      j = min(Ja,Jb)
      B(1:M,1:i,1:j) = A(1:M,1:i,1:j)
!
      end subroutine copy_2array_r
!
!-----------------------------------------------------------------------
!
      subroutine copy_3array_r(M,A,Ia,Ja,Ka, B,Ib,Jb,Kb)
!
      implicit none
      integer :: M,Ia,Ja,Ka,Ib,Jb,Kb,i,j,k
      real(8) :: A(M,Ia,Ja,Ka)
      real(8) :: B(M,Ib,Jb,Kb)
!
      B = 0.d0
      i = min(Ia,Ib)
      j = min(Ja,Jb)
      k = min(Ka,Kb)
      B(1:M,1:i,1:j,1:k) = A(1:M,1:i,1:j,1:k)
!
      end subroutine copy_3array_r

