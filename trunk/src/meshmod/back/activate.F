c---------------------------------------------------------------------
c
c   routine name       - activate
c
c---------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - activate and initialize the dofs 
c
c   arguments :
c
c     in:
c               Nod    - nod
c
c-----------------------------------------------------------------------
c
      subroutine activate(Nod)
c
      use data_structure3D
#include "syscom.blk"
c
      select case(Nod)
      case(331)
        iprint=0
      case default
        iprint=0
      end select
      if (NODES(Nod)%act.eq.1) then
        write(*, 7077) Nod
 7077   format('activate: ERROR NOD = ',i10,' IS ALREADY ACTIVATED')
        return
      endif
c
c  ...get case number
      icase = NODES(Nod)%case
c
c  ...determine number of dof
      call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
      if (iprint.eq.1) then
        write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
 7001   format('activate: Nod=',i10,'ndofH,ndofE,ndofV,ndofQ= ',4i4)
      endif
c
c  ...allocate dofs
      if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
        nvar = NREQNH(icase)*NRCOMS
        allocate(NODES(Nod)%zdofH(nvar, ndofH))
        NODES(Nod)%zdofH = ZERO
        NRDOFSH = NRDOFSH + ndofH*NREQNH(icase)
      endif
      if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
        nvar = NREQNE(icase)*NRCOMS
        allocate(NODES(Nod)%zdofE(nvar, ndofE))
        NODES(Nod)%zdofE = ZERO
        NRDOFSE = NRDOFSE + ndofE*NREQNE(icase)
      endif
      if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
        nvar = NREQNV(icase)*NRCOMS
        allocate(NODES(Nod)%zdofV(nvar, ndofV))
        NODES(Nod)%zdofV = ZERO
        NRDOFSV = NRDOFSV + ndofV*NREQNV(icase)
      endif
      if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
        nvar = NREQNQ(icase)*NRCOMS
        allocate(NODES(Nod)%zdofQ(nvar, ndofQ))
        NODES(Nod)%zdofQ = ZERO
        NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(icase)
      endif
c  ...set active flag on
      NODES(Nod)%act = 1
c
      if (iprint.ge.1) then 
        write(*,*) 'activate : HAS ACTIVATED', Nod
      endif
c     
      end subroutine 
