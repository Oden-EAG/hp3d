c---------------------------------------------------------------------
c
c   routine name       - nodgen
c
c---------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - routine generates a new node 
c
c   arguments :
c
c     in:
c              Type    - type of a new node
c              Icase   - node case
c              Nbcond  - boundary condition flag
c              Nfath   - father of the node
c              Norder  - order of approximation for the new node
c              Nfilter - refinement filter
c              Iact    = 1  active node, allocate dof
c                      = 0  inactive node, DO NOT allocate dof
c              X       = coord for vertex       
c     out:
c              Nod     - number of the new node
c
c-----------------------------------------------------------------------
c
      subroutine nodgen(Type,Icase,
     .                  Nbcond,Nfath,Norder,Nfilter,Iact,X, Nod)
c
      use data_structure3D
c
#include "syscom.blk"
      dimension X(NDIMEN)
      character(len=4) :: Type
      dimension ncase(NR_PHYSA)
c
c-----------------------------------------------------------------------
c
      iprint=0
      if (iprint.eq.1) then
        write(*,7000) Type,Icase,Nbcond,Nfath,Norder,Iact
 7000   format(' nodgen: Type,Icase,Nbcond,Nfath,Norder,Iact = ',
     .                  a4,2x,i3,2x,i5,2x,i6,2x,i3,2x,i2)
      endif
c
      call decod(Icase,2,NR_PHYSA, ncase)
c
c  ...pointer to the first free entry in NODES array
      Nod=NPNODS
      if ( (Nod.eq.0) .or. (Nod.gt.MAXNODS) ) then
        write(*,7002)Nod,NRNODS,MAXNODS
7002    format(' nodgen: Nod,NRNODS,MAXNODS = ',3(i8,2x))
        write(*,*)'NO ROOM FOR A NEW NODE!!!'
        stop
      endif
c
c  ...update      
      NRNODS=NRNODS+1
      NPNODS=NODES(Nod)%bcond
c
c  ...store node information 
      NODES(Nod)%type  = Type
      NODES(Nod)%case  = Icase
      NODES(Nod)%order = Norder
      NODES(Nod)%bcond = Nbcond
c
      call set_index(Icase,Nbcond, NODES(Nod)%index)
      NODES(Nod)%ref_kind    = 0
      NODES(Nod)%ref_filter  = Nfilter
      NODES(Nod)%father      = Nfath
      NODES(Nod)%geom_interf = 0
      NODES(Nod)%visit       = 0
c
c  ...determine the number of H1, AND H(curl), AND H(div), AND L2 dofs
c     for "a" node with order of approximation NODES(Nod)%order
c     REMARK : same order of approximation is used for all the energy
c              spaces.
      call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
c     
c  ...printing
      if (iprint.eq.1) then
        write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
 7001   format('nodgen: Nod = ',i10,' ndofH,ndofE,ndofV,ndofQ = ',4i4)
      endif
c
c  ...allocate and initialize geometry dofs
      if (ndofH.gt.0) then
        allocate(NODES(Nod)%coord(NDIMEN,ndofH))
        NODES(Nod)%coord=0.d0
      endif
      if (Type.eq.'vert')  NODES(Nod)%coord(1:NDIMEN,1)=X(1:NDIMEN)
c
c  ...allocate and initialize solution dofs      
      if (Iact.eq.1) then
        if ((NREQNH(Icase).gt.0).and.(ndofH.gt.0)) then
          nvar = NREQNH(Icase)*NRCOMS
          allocate( NODES(Nod)%zdofH(nvar, ndofH))
          NODES(Nod)%zdofH = ZERO
          NRDOFSH = NRDOFSH + ndofH*NREQNH(Icase)
        endif
        if ((NREQNE(Icase).gt.0).and.(ndofE.gt.0)) then
          nvar = NREQNE(Icase)*NRCOMS
          allocate( NODES(Nod)%zdofE(nvar, ndofE))
          NODES(Nod)%zdofE = ZERO
          NRDOFSE = NRDOFSE + ndofE*NREQNE(Icase)
        endif
        if ((NREQNV(Icase).gt.0).and.(ndofV.gt.0)) then
          nvar = NREQNV(Icase)*NRCOMS
          allocate( NODES(Nod)%zdofV(nvar, ndofV))
          NODES(Nod)%zdofV = ZERO
          NRDOFSV = NRDOFSV + ndofV*NREQNV(Icase)
        endif
        if ((NREQNQ(Icase).gt.0).and.(ndofQ.gt.0)) then
          nvar = NREQNQ(Icase)*NRCOMS
          allocate( NODES(Nod)%zdofQ(nvar, ndofQ))
          NODES(Nod)%zdofQ = ZERO
          NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(Icase)
        endif
      endif
c
c  ...activation flag      
      NODES(Nod)%act=Iact
c      
c  ...printing
      if (iprint.eq.1) then
        write(*,*) 'nodgen: Nod ',Nod,'HAS BEEN GENERATED'
        write(*,*) '        Type ',Type, ', Norder = ',Norder
        call pause
      endif
c
c      
      end
