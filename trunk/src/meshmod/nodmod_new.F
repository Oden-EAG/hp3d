
C THIS ROUTINE IS STILL INCOMPLETE. IT WILL EVENTUALLY SUBSTITUTE THE
C CURRENT VERSION OF NODMOD.
C THE PURPOSE OF THIS UPDATED VERSION IS TO UPDATE CORRECTLY THE DOFS
C IN THE DATASTRUCTURE AFTER A NODE IS P-REFINED


c---------------------------------------------------------------------
c
c   routine name       - nodmod_new
c
c---------------------------------------------------------------------
c
c   latest revision    - Aug 17
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
c
#include "syscom.blk"
c
      dimension xnod(NDIMEN,MAXmdlbH)
      dimension zdofH(MAXEQNH,MAXmdlbH),zdofE(MAXEQNE,MAXmdlbE),
     .          zdofV(MAXEQNV,MAXmdlbV),zdofQ(MAXEQNQ,MAXmdlbQ)
      character(len=4) :: ntype
c
      select case(Nod)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
c
c  ...node case
      icase = NODES(Nod)%case
      ntype = NODES(Nod)%type
      nordo = NODES(Nod)%order
      if (iprint.eq.1) then
        write(*,7001) Nod,ntype,nordo,icase
 7001   format('nodmod: Nod,type,order,icase = ',i6,2x,a5,2x,2i3)
      endif
      if (Nordn.eq.nordo) return
c
c  ...determine the current number of dof
      call ndof_nod(ntype,nordo,
     .              ndofH,ndofE,ndofV,ndofQ)
c
c  ...save current gdofs in local array and update number of gdofs
      if (Associated(NODES(Nod)%coord)) then
        xnod(1:NDIMEN,1:ndofH) = NODES(Nod)%coord(1:NDIMEN,1:ndofH)
        deallocate(NODES(Nod)%coord)
      endif
c
c  ...save current dof in local arrays
      if (Associated(NODES(Nod)%zdofH)) then
        nvar = NREQNH(icase)*NRCOMS
        zdofH(1:nvar,1:ndofH) = NODES(Nod)%zdofH(1:nvar,1:ndofH)
        deallocate(NODES(Nod)%zdofH)
        NRDOFSH = NRDOFSH - ndofH*NREQNH(icase)
      endif
      if (Associated(NODES(Nod)%zdofE)) then
        nvar = NREQNE(icase)*NRCOMS
        zdofE(1:nvar,1:ndofE) = NODES(Nod)%zdofE(1:nvar,1:ndofE)
        deallocate(NODES(Nod)%zdofE)
        NRDOFSE = NRDOFSE - ndofE*NREQNE(icase)
      endif
      if (Associated(NODES(Nod)%zdofV)) then
        nvar = NREQNV(icase)*NRCOMS
        zdofV(1:nvar,1:ndofV) = NODES(Nod)%zdofV(1:nvar,1:ndofV)
        deallocate(NODES(Nod)%zdofV)
        NRDOFSV = NRDOFSV - ndofV*NREQNV(icase)
      endif
      if (Associated(NODES(Nod)%zdofQ)) then
        nvar = NREQNQ(icase)*NRCOMS
        zdofQ(1:nvar,1:ndofQ) = NODES(Nod)%zdofQ(1:nvar,1:ndofQ)
        deallocate(NODES(Nod)%zdofQ)
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
c
c  ...allocate memory for geometry and solution dof
c     and copy old dof from the local arrays
      if (ndofH.gt.0) then
        allocate(NODES(Nod)%coord(NDIMEN,ndofH))
        call copy_dofG(ntype,nordo,Nordn,ndofHo,ndofH,
     .                 xnod,NODES(Nod)%coord)
      endif
      if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
        nvar = NREQNH(icase)*NRCOMS
        allocate( NODES(Nod)%zdofH(nvar, ndofH))
        NRDOFSH = NRDOFSH + ndofH*NREQNH(icase)
        call copy_dofH(ntype,nordo,Nordn,ndofHo,ndofH,
     .                 zdofH,NODES(Nod)%zdofH)
      endif
      if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
        nvar = NREQNE(icase)*NRCOMS
        allocate( NODES(Nod)%zdofE(nvar, ndofE))
        NRDOFSE = NRDOFSE + ndofE*NREQNE(icase)
        NODES(Nod)%zdofE(1:nvar, 1:ndofE) = ZERO
      endif
      if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
        nvar = NREQNV(icase)*NRCOMS
        allocate( NODES(Nod)%zdofV(nvar, ndofV))
        NRDOFSV = NRDOFSV + ndofV*NREQNV(icase)
        NODES(Nod)%zdofV(1:nvar, 1:ndofV) = ZERO
      endif
      if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
        nvar = NREQNQ(icase)*NRCOMS
        allocate( NODES(Nod)%zdofQ(nvar, ndofQ))
        NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(icase)
        call copy_dofQ(ntype,nordo,Nordn,ndofQo,ndofQ,
     .                 zdofQ,NODES(Nod)%zdofQ)
      endif
c
c
      if (iprint.eq.1) then
        write(*,7011) Nod
 7011   format('nodmod: Nod = ',i6,' HAS BEEN UPDATED ')
        write(*,7012) nordo,Nordn
 7012   format('        OLD ORDER = ',i3,' NEW ORDER = ',i3)
        call pause
      endif
c
c
      end
c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofG
c
c---------------------------------------------------------------------
c
c   latest revision    - Aug 17
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
      use data_structure3D
#include "syscom.blk"
      character(len=4) :: Ntype
c
      dimension Xnodo(3,NdofGo),Xnodn(3,NdofGn)
c
      Xnodn = 0.d0
      select case(Ntype)
c
c  ...edge, triangle, tet, pyramid
      case('medg','mdlt','mdln','mdld')
        ndofG = min(NdofGn,NdofGo)
        Xnodn(1:3,1:NdofG) = Xnodo(1:3,1:NdofG)
c
c  ...quad, prism
      case('mdlq','mdlp')
        call decode(Nordo, nordho,nordvo)
        call decode(Nordn, nordhn,nordvn)
        select case(Ntype)
        case('mdlq')
          ndofGhn = nordhn-1; ndofGho = nordho-1
        case('mdlp')
          ndofGhn = (nordhn-2)*(nordhn-1)/2
          ndofGho = (nordho-2)*(nordho-1)/2
        end select
        ndofGh = min(ndofGhn,ndofGho)
        ndofGv = min(nordvo-1,nordvn-1)
        ibegn=0; ibego=0
        do j=1,ndofGv
          Xnodn(1:3,ibegn+1:ibegn+ndofGh)
     .  = Xnodo(1:3,ibego+1:ibego+ndofGo)
c
c  .......update counters
          ibegn = ibegn + ndofGn
          ibego = ibego + ndofGo
        enddo
c
c  ...hexa
      case('mdlb')
        call decode(Nordo, nordho,nord3o)
        call decode(nordho, nord1o,nord2o)
        call decode(Nordn, nordhn,nord3n)
        call decode(nordhn, nord1n,nord2n)
        ibegn=0; ibego=0
        do k=1,min(nord3o,nord3n)-1
          do j=1,min(nord2o,nord2n)-1
            ndofG1 = min(nord1o,nord1n)-1
            Xnodn(1:3,ibegn+1:ibegn+ndofG1)
     .    = Xnodo(1:3,ibego+1:ibego+ndofG1)
c
c  .........update counters
            ibegn = ibegn + nord1n-1
            ibego = ibego + nord1o-1
          enddo
        enddo
      end select
c
c
      end subroutine
c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofH
c
c---------------------------------------------------------------------
c
c   latest revision    - Aug 17
c
c   purpose            - routine copies H1 dof for a node of order
c                        Nordo to a node of order Nordn
c
c   arguments :
c
c     in:
c        Ntype         - node type
c        Nordo,Nordn   - old and new orders of the node
c        ZdofHo,ZdofHn - old and new dof
c        NdofHo,NdofHn - old and new # dof
c
c-----------------------------------------------------------------------
c
      subroutine copy_dofH(Ntype,Nordo,Nordn,NdofHo,NdofHn,
     .                     ZdofHo,ZdofHn)
c
      use data_structure3D
#include "syscom.blk"
      character(len=4) :: Ntype
c
      dimension ZdofHo(MAXEQNH,NdofHo),ZdofHn(MAXEQNH,NdofHn)
c
      ZdofHn = ZERO
      select case(Ntype)
c
c  ...edge, triangle, tet, pyramid
      case('medg','mdlt','mdln','mdld')
        ndofH = min(NdofHn,NdofHo)
        ZdofHn(1:MAXEQNH,1:NdofH) = ZdofHo(1:MAXEQNH,1:NdofH)
c
c  ...quad, prism
      case('mdlq','mdlp')
        call decode(Nordo, nordho,nordvo)
        call decode(Nordn, nordhn,nordvn)
        select case(Ntype)
        case('mdlq')
          ndofHhn = nordhn-1; ndofHho = nordho-1
        case('mdlp')
          ndofHhn = (nordhn-2)*(nordhn-1)/2
          ndofHho = (nordho-2)*(nordho-1)/2
        end select
        ndofHh = min(ndofHhn,ndofHho)
        ndofHv = min(nordvo-1,nordvn-1)
        ibegn=0; ibego=0
        do j=1,ndofHv
          ZdofHn(1:MAXEQNH,ibegn+1:ibegn+ndofHh)
     .  = ZdofHo(1:MAXEQNH,ibego+1:ibego+ndofHh)
c
c  .......update counters
          ibegn = ibegn + ndofHn
          ibego = ibego + ndofHo
        enddo
c
c  ...hexa
      case('mdlb')
        call decode(Nordo, nordho,nord3o)
        call decode(nordho, nord1o,nord2o)
        call decode(Nordn, nordhn,nord3n)
        call decode(nordhn, nord1n,nord2n)
        ibegn=0; ibego=0
        do k=1,min(nord3o,nord3n)-1
          do j=1,min(nord2o,nord2n)-1
            ndofH1 = min(nord1o,nord1n)-1
            ZdofHn(1:MAXEQNH,ibegn+1:ibegn+ndofH1)
     .    = ZdofHo(1:MAXEQNH,ibego+1:ibego+ndofH1)
c
c  .........update counters
            ibegn = ibegn + nord1n-1
            ibego = ibego + nord1o-1
          enddo
        enddo
      end select
c
c
      end subroutine
c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofE
c
c---------------------------------------------------------------------
c
c   latest revision    - Aug 17
c
c   purpose            - routine copies H(curl) dof for a node of order
c                        Nordo to a node of order Nordn
c
c   arguments :
c
c     in:
c        Ntype         - node type
c        Nordo,Nordn   - old and new orders of the node
c        ZdofEo,ZdofEn - old and new dof
c        NdofEo,NdofEn - old and new # dof
c
c-----------------------------------------------------------------------
c
      subroutine copy_dofE(Ntype,Nordo,Nordn,NdofEo,NdofEn,
     .                     ZdofEo,ZdofEn)
c
      use data_structure3D
#include "syscom.blk"
      character(len=4) :: Ntype
c
      dimension ZdofEo(MAXEQNE,NdofEo),ZdofEn(MAXEQNE,NdofEn)
c
      ZdofEn = ZERO
      select case(Ntype)
      case('medg')
        ndofE = min(NdofEo,NdofEn)
        ZdofEn(1:MAXEQNE,1:NdofE) = ZdofEo(1:MAXEQNE,1:NdofE)
      case('mdlt')
      end select
c
c
      end subroutine























c
c---------------------------------------------------------------------
c
c   routine name       - copy_dofQ
c
c---------------------------------------------------------------------
c
c   latest revision    - Aug 17
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
c        NdofQo,NdofQn - old and new # dof
c
c-----------------------------------------------------------------------
c
      subroutine copy_dofQ(Ntype,Nordo,Nordn,NdofQo,NdofQn,
     .                     ZdofQo,ZdofQn)
c
      use data_structure3D
#include "syscom.blk"
      character(len=4) :: Ntype
c
      dimension ZdofQo(MAXEQNQ,NdofQo),ZdofQn(MAXEQNQ,NdofQn)
c
      ZdofQn = ZERO
      select case(Ntype)
c
c  ...tet, pyramid
      case('mdln','mdld')
        ndofQ = min(NdofQn,NdofQo)
        ZdofQn(1:MAXEQNQ,1:NdofQ) = ZdofQo(1:MAXEQNQ,1:NdofQ)
c
c  ...prism
      case('mdlp')
        call decode(Nordo, nordho,nordvo)
        call decode(Nordn, nordhn,nordvn)
        ndofQhn = (nordhn+1)*nordhn/2
        ndofQho = (nordho+1)*nordho/2
        ndofQh = min(ndofQhn,ndofQho)
        nordv = min(nordvo,nordvn)
        ibegn=0; ibego=0
        do j=1,nordv
          ZdofQn(1:MAXEQNQ,ibegn+1:ibegn+ndofQh)
     .  = ZdofQo(1:MAXEQNQ,ibego+1:ibego+ndofQh)
c
c  .......update counters
          ibegn = ibegn + ndofQn
          ibego = ibego + ndofQo
        enddo
c
c  ...hexa
      case('mdlb')
        call decode(Nordo, nordho,nord3o)
        call decode(nordho, nord1o,nord2o)
        call decode(Nordn, nordhn,nord3n)
        call decode(nordhn, nord1n,nord2n)
        ibegn=0; ibego=0
        do k=1,min(nord3o,nord3n)
          do j=1,min(nord2o,nord2n)
            ndofQ1 = min(nord1o,nord1n)
            ZdofQn(1:MAXEQNQ,ibegn+1:ibegn+ndofQ1)
     .    = ZdofQo(1:MAXEQNQ,ibego+1:ibego+ndofQ1)
c
c  .........update counters
            ibegn = ibegn + nord1n
            ibego = ibego + nord1o
          enddo
        enddo
      end select
c
c
      end subroutine
















