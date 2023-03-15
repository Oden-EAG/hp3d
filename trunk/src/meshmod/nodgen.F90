!---------------------------------------------------------------------
!
!   routine name       - nodgen
!
!---------------------------------------------------------------------
!
!
!> @brief      routine generates a new node
!!
!> @param[in]  Ntype   - type of a new node
!> @param[in]  Icase   - node case
!> @param[in]  Nbcond  - boundary condition flag
!> @param[in]  Nfath   - father of the node
!> @param[in]  Norder  - order of approximation for the new node
!> @param[in]  Subd    - subdomain of node
!> @param[in]  Iact    = T  active node, allocate dof
!!                     = F  inactive node, DO NOT allocate dof
!> @param[out] Nod     - number of the new node
!!
!> @date       Feb 2023
!-----------------------------------------------------------------------
!
subroutine nodgen(Ntype,Icase,Nbcond,Nfath,Norder,Subd,Iact, Nod)
!
   use data_structure3D
   use environment, only: QUIET_MODE
   use mpi_param  , only: RANK,ROOT
!
   implicit none
!
   integer, intent(in)  :: Ntype
   integer, intent(in)  :: Icase
   integer, intent(in)  :: Nbcond
   integer, intent(in)  :: Nfath
   integer, intent(in)  :: Norder
   integer, intent(in)  :: Subd
   logical, intent(in)  :: Iact
   integer, intent(out) :: Nod
!
   integer :: ndofH,ndofE,ndofV,ndofQ,nvar
!
#if DEBUG_MODE
   integer :: ncase(NR_PHYSA)
   integer :: ibcnd(NRINDEX)
   integer :: iprint
   iprint=0
#endif
!
!-----------------------------------------------------------------------
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7000) S_Type(Ntype),Icase,Nbcond,Nfath,Norder,Iact
 7000 format(' nodgen: Type,Icase,Nbcond,Nfath,Norder,Iact = ', &
                        a4,2x,i3,2x,i6,2x,i6,2x,i3,2x,l2)
      call decod(Icase,2,NR_PHYSA, ncase)
      write(*,7001) ncase(1:NR_PHYSA)
 7001 format(' decoded Icase = ',10i1)
      call decod(Nbcond,2,NRINDEX, ibcnd)
      write(*,7002) ibcnd(1:NRINDEX)
 7002 format(' decoded Nbcond = ',30i1)
   endif
#endif
!
   if ((NPNODS.eq.0) .or. (NPNODS.gt.MAXNODS)) then
      if (.not. QUIET_MODE .and. RANK.eq.ROOT) then
         write(*,*) 'nodgen: increasing size of NODES array.'
      endif
      call increase_MAXNODS
   endif
!
!..pointer to the first free entry in NODES array
   Nod=NPNODS
!
!..update
   NRNODS=NRNODS+1
   NPNODS=NODES(Nod)%bcond
!
!..store node information
   NODES(Nod)%ntype = Ntype
   NODES(Nod)%case  = Icase
   NODES(Nod)%order = Norder
   NODES(Nod)%bcond = Nbcond
!
   NODES(Nod)%ref_kind    = 0
   NODES(Nod)%father      = Nfath
   NODES(Nod)%visit       = 0
   NODES(Nod)%subd        = Subd
!
!..determine the number of H1, AND H(curl), AND H(div), AND L2 dofs
!  for "a" node with order of approximation NODES(Nod)%order
!  REMARK : same order of approximation is used for all the energy spaces.
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
!
!..printing
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7005) Nod,ndofH,ndofE,ndofV,ndofQ
 7005 format(' nodgen: Nod = ',i10,' ndofH,ndofE,ndofV,ndofQ = ',4i4)
   endif
#endif
!
!..allocate dof data type
   if (Iact) then
      allocate(NODES(Nod)%dof)
      nullify(NODES(Nod)%dof%coord)
      nullify(NODES(Nod)%dof%zdofH)
      nullify(NODES(Nod)%dof%zdofE)
      nullify(NODES(Nod)%dof%zdofV)
      nullify(NODES(Nod)%dof%zdofQ)
   endif
!
!..allocate and initialize geometry dofs
   if (Iact .and. (ndofH.gt.0)) then
      allocate(NODES(Nod)%dof%coord(NDIMEN,ndofH))
      NODES(Nod)%dof%coord=0.d0
   endif
!
!..allocate and initialize solution dofs
   if (Iact) then
      if ((NREQNH(Icase).gt.0).and.(ndofH.gt.0)) then
         nvar = NREQNH(Icase)*NRCOMS
         allocate( NODES(Nod)%dof%zdofH(nvar,ndofH))
         NODES(Nod)%dof%zdofH = ZERO
         NRDOFSH = NRDOFSH + ndofH*NREQNH(Icase)
      endif
      if ((NREQNE(Icase).gt.0).and.(ndofE.gt.0)) then
         nvar = NREQNE(Icase)*NRCOMS
         allocate( NODES(Nod)%dof%zdofE(nvar,ndofE))
         NODES(Nod)%dof%zdofE = ZERO
         NRDOFSE = NRDOFSE + ndofE*NREQNE(Icase)
      endif
      if ((NREQNV(Icase).gt.0).and.(ndofV.gt.0)) then
         nvar = NREQNV(Icase)*NRCOMS
         allocate( NODES(Nod)%dof%zdofV(nvar,ndofV))
         NODES(Nod)%dof%zdofV = ZERO
         NRDOFSV = NRDOFSV + ndofV*NREQNV(Icase)
      endif
      if ((NREQNQ(Icase).gt.0).and.(ndofQ.gt.0)) then
         nvar = NREQNQ(Icase)*NRCOMS
         allocate( NODES(Nod)%dof%zdofQ(nvar,ndofQ))
         NODES(Nod)%dof%zdofQ = ZERO
         NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(Icase)
      endif
   endif
!
!..activation flag
   NODES(Nod)%act=Iact
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,*) 'nodgen: Nod ',Nod,'HAS BEEN GENERATED'
      write(*,*) '        Type ',S_Type(Ntype), ', Norder = ',Norder
      call pause
   endif
#endif
!
!
end subroutine nodgen
