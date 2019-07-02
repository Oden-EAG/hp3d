!-------------------------------------------------------------------------
!> Purpose : activate a node and initialize its dofs
!!
!> @param[in] Nod - node to be activated
!-------------------------------------------------------------------------
!
subroutine activate(Nod)
  use data_structure3D
  implicit none
  integer, intent(in) :: Nod
  integer :: iprint, icase, nvar, ndofH, ndofE, ndofV, ndofQ
!
!-------------------------------------------------------------------------
!
  iprint=0

  ! if node is already active, do nothing
  if (NODES(Nod)%act.eq.1) then
     return
  endif

  ! determine number of dofs, update geometry dofs
  call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)

  ! printing
  if (iprint.eq.1) then
     write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
7001 format('activate: Nod=',i10,'ndofH,ndofE,ndofV,ndofQ= ',4i4)
  endif

  ! allocate and initialize dofs
  if (ndofH.gt.0) then
     allocate(NODES(Nod)%coord(3, ndofH))
     NODES(Nod)%coord = 0.d0
  endif
  icase = NODES(Nod)%case
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

  ! raise activation flag
  NODES(Nod)%act=1

  ! printing
  if (iprint.ge.1) then
     write(*,*) 'activate : ACTIVATED Nod = ', Nod
  endif


end subroutine activate


!-------------------------------------------------------------------------
!> Purpose : deactivate a node
!!
!> @param[in] Nod - node to be deactivated
!-------------------------------------------------------------------------
!  REMARK (Sep 12) : this routine is not used, hence we never subtract
!                    geometry dofs pertaining to an element that has
!                    been refined
!-------------------------------------------------------------------------

subroutine deactivate(Nod)
  use data_structure3D
  implicit none
  integer, intent(in) :: Nod

  integer :: iprint, icase, ndofH, ndofE, ndofV, ndofQ,nn,nn1

!-------------------------------------------------------------------------

  iprint=0

  ! if node is already inactive, do nothing
  if (NODES(Nod)%act.eq.0) then
     return
  endif

  ! find number of dofs associated to node, update number of gdofs
  call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)

  ! printing
  if (iprint.eq.1) then
     write(*,7001) Nod,ndofH,ndofE,ndofV,ndofQ
7001 format('activate: Nod=',i10,'ndofH,ndofE,ndofV,ndofQ= ',4i4)
  endif

!!!!  ! deallocate geometry dof
  if (ndofH.gt.0) then
     deallocate(NODES(Nod)%coord)
  endif
  ! deallocate dofs
  icase=NODES(Nod)%case

  if ((NREQNH(icase).gt.0).and.(ndofH.gt.0)) then
     deallocate(NODES(Nod)%zdofH)
     NRDOFSH = NRDOFSH - ndofH*NREQNH(icase)
  endif
  if ((NREQNE(icase).gt.0).and.(ndofE.gt.0)) then
     deallocate(NODES(Nod)%zdofE)
     NRDOFSE = NRDOFSE - ndofE*NREQNE(icase)
  endif
  if ((NREQNV(icase).gt.0).and.(ndofV.gt.0)) then
     deallocate(NODES(Nod)%zdofV)
     NRDOFSV = NRDOFSV - ndofV*NREQNV(icase)
  endif
  if ((NREQNQ(icase).gt.0).and.(ndofQ.gt.0)) then
     deallocate(NODES(Nod)%zdofQ)
     NRDOFSQ = NRDOFSQ - ndofQ*NREQNQ(icase)
  endif

  ! lower activation flag
  NODES(Nod)%act=0

  ! printing
  if (iprint.ge.1) then
     write(*,*) 'deactivate : DEACTIVATED Nod = ', Nod
  endif


end subroutine deactivate
