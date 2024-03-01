!-------------------------------------------------------------------------
!> @brief         activate a node and initialize its dofs.
!> @note          if mesh is distributed, routine will only allocate dofs
!!                if the node is within subdomain (NODES(Nod)%subd==rank).
!!
!> @param[in]     Nod          - node to be activated
!> @param[inout]  NrdofH,E,V,Q - Allocated number of solution dofs
!!
!> @date          Sep 2023
!-------------------------------------------------------------------------
subroutine activate(Nod, NrdofH,NrdofE,NrdofV,NrdofQ)
   use data_structure3D
   use par_mesh  , only: DISTRIBUTED
   use mpi_param , only: RANK
!
   implicit none
!
   integer, intent(in)    :: Nod
   integer, intent(inout) :: NrdofH,NrdofE,NrdofV,NrdofQ
!
   integer :: icase,nvar,ndofH,ndofE,ndofV,ndofQ
   integer :: act_dof,subd
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!-------------------------------------------------------------------------
!
!..if node is already active, do nothing
   if (Is_active(Nod)) return
!
!..determine number of dofs, update geometry dofs
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
 7001 format('activate: Nod=',i10,', ndofH,ndofE,ndofV,ndofQ= ',4i4)
   endif
#endif
!
!..find out whether node is inside my subdomain
   act_dof = 1
   if (DISTRIBUTED) then
      call get_subd(Nod, subd)
      if (RANK .ne. subd) act_dof = 0
   endif
!
!..allocate dof data type
   if (act_dof.eq.1) then
      allocate(NODES(Nod)%dof)
      nullify(NODES(Nod)%dof%coord)
      nullify(NODES(Nod)%dof%zdofH)
      nullify(NODES(Nod)%dof%zdofE)
      nullify(NODES(Nod)%dof%zdofV)
      nullify(NODES(Nod)%dof%zdofQ)
   endif
!
!..allocate and initialize dofs
   if ((ndofH.gt.0) .and. (act_dof.eq.1)) then
      allocate(NODES(Nod)%dof%coord(3, ndofH))
      NODES(Nod)%dof%coord = 0.d0
   endif
   icase = NODES(Nod)%case
   if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
      nvar = NREQNH(icase)*NRRHS
      if(act_dof.eq.1) then
         allocate(NODES(Nod)%dof%zdofH(nvar, ndofH, NRCOMS))
         NODES(Nod)%dof%zdofH = ZERO
      endif
      NrdofH = NrdofH + ndofH*NREQNH(icase)
   endif
   if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
      nvar = NREQNE(icase)*NRRHS
      if(act_dof.eq.1) then
         allocate(NODES(Nod)%dof%zdofE(nvar, ndofE, NRCOMS))
         NODES(Nod)%dof%zdofE = ZERO
      endif
      NrdofE = NrdofE + ndofE*NREQNE(icase)
   endif
   if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
      nvar = NREQNV(icase)*NRRHS
      if(act_dof.eq.1) then
         allocate(NODES(Nod)%dof%zdofV(nvar, ndofV, NRCOMS))
         NODES(Nod)%dof%zdofV = ZERO
      endif
      NrdofV = NrdofV + ndofV*NREQNV(icase)
   endif
   if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
      nvar = NREQNQ(icase)*NRRHS
      if(act_dof.eq.1) then
         allocate(NODES(Nod)%dof%zdofQ(nvar, ndofQ, NRCOMS))
         NODES(Nod)%dof%zdofQ = ZERO
      endif
      NrdofQ = NrdofQ + ndofQ*NREQNQ(icase)
   endif
!
!..raise activation flag
   NODES(Nod)%act=.true.
!
#if HP3D_DEBUG
   if (iprint.ge.1) then
      write(*,*) 'activate : ACTIVATED Nod = ', Nod
   endif
#endif
!
end subroutine activate
!
!
!-------------------------------------------------------------------------
!> @brief deactivate a node.
!            note: dofs are only deallocated if node was in subdomain,
!                  i.e., if dof pointers were associated with data.
!
!> @date Aug 2019
!
!> @param[in]    Nod          - node to be deactivated
!> @param[inout] NrdofH,E,V,Q - Deallocated number of solution dofs
!-------------------------------------------------------------------------
subroutine deactivate(Nod, NrdofH,NrdofE,NrdofV,NrdofQ)
   use data_structure3D
!
   implicit none
!
   integer, intent(in)    :: Nod
   integer, intent(inout) :: NrdofH,NrdofE,NrdofV,NrdofQ
!
   integer :: icase,ndofH,ndofE,ndofV,ndofQ
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!-------------------------------------------------------------------------
!
!..if node is already inactive, do nothing
   if (Is_inactive(Nod)) return
!
!..find number of dofs associated to node, update number of gdofs
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
 7001 format('deactivate: Nod=',i10,'ndofH,ndofE,ndofV,ndofQ= ',4i4)
   endif
#endif
!
!..deallocate geometry dof
   if (associated(NODES(Nod)%dof)) then
      if (ndofH.gt.0 .and. associated(NODES(Nod)%dof%coord)) then
         deallocate(NODES(Nod)%dof%coord)
      endif
   endif
!
!..deallocate solution dofs
   icase=NODES(Nod)%case
!
   if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
      if (associated(NODES(Nod)%dof)) then
         if (associated(NODES(Nod)%dof%zdofH)) then
            deallocate(NODES(Nod)%dof%zdofH)
         endif
      endif
      NrdofH = NrdofH - ndofH*NREQNH(icase)
   endif
   if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
      if (associated(NODES(Nod)%dof)) then
         if (associated(NODES(Nod)%dof%zdofE)) then
            deallocate(NODES(Nod)%dof%zdofE)
         endif
      endif
      NrdofE = NrdofE - ndofE*NREQNE(icase)
   endif
   if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
      if (associated(NODES(Nod)%dof)) then
         if (associated(NODES(Nod)%dof%zdofV)) then
            deallocate(NODES(Nod)%dof%zdofV)
         endif
      endif
      NrdofV = NrdofV - ndofV*NREQNV(icase)
   endif
   if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
      if (associated(NODES(Nod)%dof)) then
         if (associated(NODES(Nod)%dof%zdofQ)) then
            deallocate(NODES(Nod)%dof%zdofQ)
         endif
      endif
      NrdofQ = NrdofQ - ndofQ*NREQNQ(icase)
   endif
   if (associated(NODES(Nod)%dof)) then
      deallocate(NODES(Nod)%dof)
   endif
!
!..lower activation flag
   NODES(Nod)%act=.false.
!
#if HP3D_DEBUG
   if (iprint.ge.1) then
      write(*,*) 'deactivate : DEACTIVATED Nod = ', Nod
   endif
#endif
!
end subroutine deactivate
