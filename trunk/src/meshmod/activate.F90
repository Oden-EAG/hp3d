!-------------------------------------------------------------------------
!> Purpose : activate a node and initialize its dofs.
!            note: if mesh is distributed, routine will only allocate dofs
!            if the node is within subdomain (NODES(Nod)%subd == rank).
!!
!> @date July 2019
!!
!> @param[in] Nod - node to be activated
!-------------------------------------------------------------------------
subroutine activate(Nod)
   use data_structure3D
   use par_mesh  , only: DISTRIBUTED
   use MPI_param , only: RANK
!
   implicit none
!
   integer, intent(in) :: Nod
   integer :: iprint, icase, nvar, ndofH, ndofE, ndofV, ndofQ
   integer :: act_dof, subd
!
!-------------------------------------------------------------------------
!
   iprint=0
!
!..if node is already active, do nothing
   if (NODES(Nod)%act.eq.1) then
      return
   endif
!
!..determine number of dofs, update geometry dofs
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
!
!..printing
   if (iprint.eq.1) then
      write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
 7001 format('activate: Nod=',i10,'ndofH,ndofE,ndofV,ndofQ= ',4i4)
   endif
!
!..find out whether node is inside my subdomain
   act_dof = 1
   if (DISTRIBUTED) then
      call get_subd(Nod, subd)
      if (RANK .ne. subd) act_dof = 0
   endif
!
!..allocate and initialize dofs
   if (ndofH.gt.0) then
      allocate(NODES(Nod)%coord(3, ndofH))
      NODES(Nod)%coord = 0.d0
   endif
   icase = NODES(Nod)%case
   if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
      nvar = NREQNH(icase)*NRCOMS
      if(act_dof.eq.1) then
         allocate(NODES(Nod)%zdofH(nvar, ndofH))
         NODES(Nod)%zdofH = ZERO
      endif
      NRDOFSH = NRDOFSH + ndofH*NREQNH(icase)
   endif
   if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
      nvar = NREQNE(icase)*NRCOMS
      if(act_dof.eq.1) then
         allocate(NODES(Nod)%zdofE(nvar, ndofE))
         NODES(Nod)%zdofE = ZERO
      endif
      NRDOFSE = NRDOFSE + ndofE*NREQNE(icase)
   endif
   if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
      nvar = NREQNV(icase)*NRCOMS
      if(act_dof.eq.1) then
         allocate(NODES(Nod)%zdofV(nvar, ndofV))
         NODES(Nod)%zdofV = ZERO
      endif
      NRDOFSV = NRDOFSV + ndofV*NREQNV(icase)
   endif
   if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
      nvar = NREQNQ(icase)*NRCOMS
      if(act_dof.eq.1) then
         allocate(NODES(Nod)%zdofQ(nvar, ndofQ))
         NODES(Nod)%zdofQ = ZERO
      endif
      NRDOFSQ = NRDOFSQ + ndofQ*NREQNQ(icase)
   endif
!
!..raise activation flag
   NODES(Nod)%act=1
!
!..printing
   if (iprint.ge.1) then
      write(*,*) 'activate : ACTIVATED Nod = ', Nod
   endif
!
end subroutine activate
!
!
!-------------------------------------------------------------------------
!> Purpose : deactivate a node.
!            note: dofs are only deallocated if node was in subdomain,
!                  i.e., if dof pointers were associated with data.
!!
!> @param[in] Nod - node to be deactivated
!-------------------------------------------------------------------------
subroutine deactivate(Nod)
   use data_structure3D
!
   implicit none
!
   integer, intent(in) :: Nod
!
   integer :: iprint, icase, ndofH, ndofE, ndofV, ndofQ,nn,nn1
!
!-------------------------------------------------------------------------
!
   iprint=0
!
!..if node is already inactive, do nothing
   if (NODES(Nod)%act.eq.0) then
      return
   endif
!
!..find number of dofs associated to node, update number of gdofs
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
!
!..printing
   if (iprint.eq.1) then
      write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
 7001 format('activate: Nod=',i10,'ndofH,ndofE,ndofV,ndofQ= ',4i4)
   endif
!
!..deallocate geometry dof
   if (ndofH.gt.0 .and. associated(NODES(Nod)%coord)) then
      deallocate(NODES(Nod)%coord)
   endif
!
!..deallocate solution dofs
   icase=NODES(Nod)%case
!
   if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
      if (associated(NODES(Nod)%zdofH)) then
         deallocate(NODES(Nod)%zdofH)
      endif
      NRDOFSH = NRDOFSH - ndofH*NREQNH(icase)
   endif
   if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
      if (associated(NODES(Nod)%zdofE)) then
         deallocate(NODES(Nod)%zdofE)
      endif
      NRDOFSE = NRDOFSE - ndofE*NREQNE(icase)
   endif
   if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
      if (associated(NODES(Nod)%zdofV)) then
         deallocate(NODES(Nod)%zdofV)
      endif
      NRDOFSV = NRDOFSV - ndofV*NREQNV(icase)
   endif
   if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
      if (associated(NODES(Nod)%zdofQ)) then
         deallocate(NODES(Nod)%zdofQ)
      endif
      NRDOFSQ = NRDOFSQ - ndofQ*NREQNQ(icase)
   endif
!
!..lower activation flag
   NODES(Nod)%act=0
!
!..printing
   if (iprint.ge.1) then
      write(*,*) 'deactivate : DEACTIVATED Nod = ', Nod
   endif
!
end subroutine deactivate
