c---------------------------------------------------------------------
c
c   routine name       - nodmod
c
c---------------------------------------------------------------------
c
c   latest revision    - July 2019
c
c   purpose            - routine modifies order of approximation
c                        for a node.
c                  note: if mesh is distributed, routine will
c                        only allocate dofs if the node is within
c                        subdomain (NODES(Nod)%subd == rank).
c                remark: routine does not currently update the dofs
c                        correctly after a node is p-refined
c                        (see also: nodmod_new.F)
c
c                  TODO: change to F90 implementation
c                        fix dof update (see remark)
c
c   arguments :
c
c     in:
c             Nod      - node number
c             Newp     - new order of approximation
c
c-----------------------------------------------------------------------
c
      subroutine nodmod(Nod,Newp)
c
      use data_structure3D
      use element_data
      use MPI_param, only: RANK
      use par_mesh , only: DISTRIBUTED
c
#include "syscom.blk"
c
      dimension xnod(NDIMEN,MAXmdlbH)
      dimension zdofH(MAXEQNH,MAXmdlbH),zdofE(MAXEQNE,MAXmdlbE),
     .          zdofV(MAXEQNV,MAXmdlbV),zdofQ(MAXEQNQ,MAXmdlbQ)
c
      integer :: act_dof, subd
c
      select case(Nod)
      case(2348)
        iprint=0
      case default
        iprint=0
      end select
c
c  ...exit if the order is the same
      if (Newp.eq.NODES(Nod)%order) return
c
c  ...node case
      icase = NODES(Nod)%case
      if (iprint.eq.1) then
        write(*,7001) Nod,NODES(Nod)%type,NODES(Nod)%order,icase
 7001   format('nodmod: Nod,type,order,icase = ',i6,2x,a5,2x,2i3)
      endif
c
c  ...determine the current number of dof
      call ndof_nod(NODES(Nod)%type,NODES(Nod)%order,
     .              ndofH,ndofE,ndofV,ndofQ)
c
c  ...save current gdofs in local array and update number of gdofs
      if (ndofH.gt.0) then
        if (Associated(NODES(Nod)%coord)) then
          xnod(1:NDIMEN,1:ndofH) = NODES(Nod)%coord(1:NDIMEN,1:ndofH)
          deallocate(NODES(Nod)%coord)
        endif
      endif
c      
c  ...save current dof in local arrays
      if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
        if (Associated(NODES(Nod)%zdofH)) then
          nvar = NREQNH(icase)*NRCOMS
          zdofH(1:nvar,1:ndofH) = NODES(Nod)%zdofH(1:nvar,1:ndofH)
          deallocate(NODES(Nod)%zdofH)
        endif
        NRDOFSH = NRDOFSH - ndofH*NREQNH(icase)
      endif
      if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
        if (Associated(NODES(Nod)%zdofE)) then
          nvar = NREQNE(icase)*NRCOMS
          zdofE(1:nvar,1:ndofE) = NODES(Nod)%zdofE(1:nvar,1:ndofE)
          deallocate(NODES(Nod)%zdofE)
        endif
        NRDOFSE = NRDOFSE - ndofE*NREQNE(icase)
      endif
      if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
        if (Associated(NODES(Nod)%zdofV)) then
          nvar = NREQNV(icase)*NRCOMS
          zdofV(1:nvar,1:ndofV) = NODES(Nod)%zdofV(1:nvar,1:ndofV)
          deallocate(NODES(Nod)%zdofV)
        endif
        NRDOFSV = NRDOFSV - ndofV*NREQNV(icase)
      endif
      if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
        if (Associated(NODES(Nod)%zdofQ)) then
          nvar = NREQNQ(icase)*NRCOMS
          zdofQ(1:nvar,1:ndofQ) = NODES(Nod)%zdofQ(1:nvar,1:ndofQ)
          deallocate(NODES(Nod)%zdofQ)
        endif
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
      npold = NODES(Nod)%order
      NODES(Nod)%order = Newp
c
c  ...calculate the new number of dof for the node
      call ndof_nod(NODES(Nod)%type,NODES(Nod)%order,
     .              ndofH,ndofE,ndofV,ndofQ)
c
c  ...find out whether node is inside my subdomain
      act_dof = 1
      if (DISTRIBUTED .eq. 1) then
         call get_subd(Nod, subd)
         if (RANK .ne. subd) act_dof = 0
      endif
c
c  ...allocate memory for geometry and solution dof
c     and copy old dof from the local arrays
      if (ndofH.gt.0) then
        allocate(NODES(Nod)%coord(NDIMEN,ndofH))
        NODES(Nod)%coord(1:NDIMEN, 1:ndofH) = 0.d0
        ndofHm = min(ndofH,ndofHo)
        NODES(Nod)%coord(1:NDIMEN, 1:ndofHm) = 
     .              xnod(1:NDIMEN, 1:ndofHm)
      endif
      if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
        if (act_dof.eq.1) then
          nvar = NREQNH(icase)*NRCOMS
          allocate( NODES(Nod)%zdofH(nvar, ndofH))
          NODES(Nod)%zdofH(1:nvar, 1:ndofH) = ZERO
          ndofHm = min(ndofH,ndofHo)
          NODES(Nod)%zdofH(1:nvar, 1:ndofHm) =
     .               zdofH(1:nvar, 1:ndofHm)
        endif
        NRDOFSH = NRDOFSH + ndofH*NREQNH(icase)
      endif
      if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
        if (act_dof.eq.1) then
          nvar = NREQNE(icase)*NRCOMS
          allocate( NODES(Nod)%zdofE(nvar, ndofE))
          NODES(Nod)%zdofE(1:nvar, 1:ndofE) = ZERO
          ndofEm = min(ndofE,ndofEo)
          NODES(Nod)%zdofE(1:nvar, 1:ndofEm) =
     .               zdofE(1:nvar, 1:ndofEm)
        endif
        NRDOFSE = NRDOFSE + ndofE*NREQNE(icase)
      endif
      if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
        if (act_dof.eq.1) then
          nvar = NREQNV(icase)*NRCOMS
          allocate( NODES(Nod)%zdofV(nvar, ndofV))
          NODES(Nod)%zdofV(1:nvar, 1:ndofV) = ZERO
          ndofVm = min(ndofV,ndofVo)
          NODES(Nod)%zdofV(1:nvar, 1:ndofVm) =
     .               zdofV(1:nvar, 1:ndofVm)
        endif
        NRDOFSV = NRDOFSV + ndofV*NREQNV(icase)
      endif
      if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
        if (act_dof.eq.1) then
          nvar = NREQNQ(icase)*NRCOMS
          allocate( NODES(Nod)%zdofQ(nvar, ndofQ))
          NODES(Nod)%zdofQ(1:nvar, 1:ndofQ) = ZERO
          ndofQm = min(ndofQ,ndofQo)
          NODES(Nod)%zdofQ(1:nvar, 1:ndofQm) =
     .               zdofQ(1:nvar, 1:ndofQm)
        endif
        NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(icase)
      endif
c
c
      if (iprint.eq.1) then
        write(*,7011) Nod
 7011   format('nodmod: Nod = ',i6,' HAS BEEN UPDATED ')
        write(*,7012) npold,Newp
 7012   format('        OLD ORDER = ',i3,' NEW ORDER = ',i3)
        call pause
      endif
c
c
      end


