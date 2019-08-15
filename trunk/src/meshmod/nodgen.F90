!---------------------------------------------------------------------
!
!   routine name       - nodgen
!
!---------------------------------------------------------------------
!
!   latest revision    - May 10
!
!   purpose            - routine generates a new node
!
!   arguments :
!
!     in:
!              Type    - type of a new node
!              Icase   - node case
!              Nbcond  - boundary condition flag
!              Nfath   - father of the node
!              Norder  - order of approximation for the new node
!              Nfilter - refinement filter
!              Subd    - subdomain of node
!              Iact    = 1  active node, allocate dof
!                      = 0  inactive node, DO NOT allocate dof
!              X       = coord for vertex
!     out:
!              Nod     - number of the new node
!
!-----------------------------------------------------------------------
!
subroutine nodgen(Type,Icase,Nbcond,Nfath,Norder,Nfilter,Subd,Iact, Nod)
!
   use data_structure3D
!
   implicit none
!
   character(len=4), intent(in)  :: Type
   integer, intent(in)           :: Icase
   integer, intent(in)           :: Nbcond
   integer, intent(in)           :: Nfath
   integer, intent(in)           :: Norder
   integer, intent(in)           :: Nfilter
   integer, intent(in)           :: Subd
   integer, intent(in)           :: Iact
   integer, intent(out)          :: Nod
!
   integer :: ncase(NR_PHYSA)
!
   integer :: iprint,ndofH,ndofE,ndofV,ndofQ,nvar
!
!-----------------------------------------------------------------------
!
   iprint=0
   if (iprint.eq.1) then
      write(*,7000) Type,Icase,Nbcond,Nfath,Norder,Iact
 7000 format(' nodgen: Type,Icase,Nbcond,Nfath,Norder,Iact = ', &
                        a4,2x,i3,2x,i6,2x,i6,2x,i3,2x,i2)
   endif
!
   call decod(Icase,2,NR_PHYSA, ncase)
!
!..pointer to the first free entry in NODES array
   Nod=NPNODS
   if ( (Nod.eq.0) .or. (Nod.gt.MAXNODS) ) then
      write(*,7002)Nod,NRNODS,MAXNODS
7002  format(' nodgen: Nod,NRNODS,MAXNODS = ',3(i8,2x))
      write(*,*)'NO ROOM FOR A NEW NODE!!!'
      stop
   endif
!
!..update
   NRNODS=NRNODS+1
   NPNODS=NODES(Nod)%bcond
!
!..store node information
   NODES(Nod)%type  = Type
   NODES(Nod)%case  = Icase
   NODES(Nod)%order = Norder
   NODES(Nod)%bcond = Nbcond
!
   call set_index(Icase,Nbcond, NODES(Nod)%index)
   NODES(Nod)%ref_kind    = 0
   NODES(Nod)%ref_filter  = Nfilter
   NODES(Nod)%father      = Nfath
   NODES(Nod)%geom_interf = 0
   NODES(Nod)%visit       = 0
   NODES(Nod)%subd        = Subd
!
!..determine the number of H1, AND H(curl), AND H(div), AND L2 dofs
!  for "a" node with order of approximation NODES(Nod)%order
!  REMARK : same order of approximation is used for all the energy spaces.
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
!
!..printing
   if (iprint.eq.1) then
      write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
 7001 format(' nodgen: Nod = ',i10,' ndofH,ndofE,ndofV,ndofQ = ',4i4)
   endif
!
!..allocate and initialize geometry dofs
   if ((Iact.eq.1) .and. (ndofH.gt.0)) then
      allocate(NODES(Nod)%coord(NDIMEN,ndofH))
      NODES(Nod)%coord=0.d0
   endif
!
!..allocate and initialize solution dofs
   if (Iact.eq.1) then
      if ((NREQNH(Icase).gt.0).and.(ndofH.gt.0)) then
         nvar = NREQNH(Icase)*NRCOMS
         allocate( NODES(Nod)%zdofH(nvar,ndofH))
         NODES(Nod)%zdofH = ZERO
         NRDOFSH = NRDOFSH + ndofH*NREQNH(Icase)
      endif
      if ((NREQNE(Icase).gt.0).and.(ndofE.gt.0)) then
         nvar = NREQNE(Icase)*NRCOMS
         allocate( NODES(Nod)%zdofE(nvar,ndofE))
         NODES(Nod)%zdofE = ZERO
         NRDOFSE = NRDOFSE + ndofE*NREQNE(Icase)
      endif
      if ((NREQNV(Icase).gt.0).and.(ndofV.gt.0)) then
         nvar = NREQNV(Icase)*NRCOMS
         allocate( NODES(Nod)%zdofV(nvar,ndofV))
         NODES(Nod)%zdofV = ZERO
         NRDOFSV = NRDOFSV + ndofV*NREQNV(Icase)
      endif
      if ((NREQNQ(Icase).gt.0).and.(ndofQ.gt.0)) then
         nvar = NREQNQ(Icase)*NRCOMS
         allocate( NODES(Nod)%zdofQ(nvar,ndofQ))
         NODES(Nod)%zdofQ = ZERO
         NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(Icase)
      endif
   endif
!
!..activation flag
   NODES(Nod)%act=Iact
!
!..printing
   if (iprint.eq.1) then
      write(*,*) 'nodgen: Nod ',Nod,'HAS BEEN GENERATED'
      write(*,*) '        Type ',Type, ', Norder = ',Norder
      call pause
   endif
!
!
end subroutine nodgen
