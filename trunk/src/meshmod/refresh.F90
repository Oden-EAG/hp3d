!--------------------------------------------------------------------
!> @brief activate inactive nodes which are not constrained
!          but in the current mesh
!> @date Oct 2019
!--------------------------------------------------------------------
subroutine refresh
!
   use data_structure3D
   use par_mesh,        only: DISTRIBUTED,set_subd_elem
   use mpi_param,       only: RANK,ROOT
   use bitvisit
!
   implicit none
!
   integer :: ntype
   integer :: nodesl(27),norientl(27)
   integer :: i,iel,nod,nfath,mdle,subd,ierr
   integer :: nrdofH,nrdofE,nrdofV,nrdofQ
   integer :: nvef
!
#if DEBUG_MODE
   integer :: iprint = 0
#endif
!
!--------------------------------------------------------------------
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,*) 'refresh: Begin'
   endif
#endif
!
   call update_ELEM_ORDER
   call bitvisit_init(NRNODS)
!
   if (.not. DISTRIBUTED) then
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
!..set thread local dof counters
   nrdofH=0; nrdofE=0; nrdofV=0; nrdofQ=0
!
!--------------------------------------------------------------------
! Step 1 : raise visitation flag for vertex, edge and face nodes of |
!          all active elements                                      |
!--------------------------------------------------------------------
!$omp parallel default(shared)
!$omp do private(iel,mdle,subd,type,i,nodesl,norientl) schedule(guided)
   do iel=1,NRELES_SUBD
!
      mdle = ELEM_SUBD(iel)
!
      call elem_nodes(mdle, nodesl,norientl)
!
      ntype = NODES(mdle)%ntype
      nvef = nvert(ntype)+nedge(ntype)+nface(ntype)
!
      do i=1,nvef
!     ...bitflag stores visit in bit collection (better for MPI)
         call visit(nodesl(i))
!
!     ...if node is visited by an element within my subdomain,
!        add node to my subdomain (need its dofs). this flag will
!        indicate that dofs must be allocated in activation.
         if (DISTRIBUTED) then
            call set_subd(nodesl(i),RANK)
         endif
      enddo
!
   enddo
!$omp end do
!$omp end parallel
!
!..Reduce over node flag
   if (DISTRIBUTED) then
      call reduce_visit
   endif
!
!--------------------------------------------------------------------
! Step 2: activate all inactive edge and face nodes whose father    |
!         node was not visited (i.e., activate inactive edge and    |
!         face nodes that are unconstrained)                        |
!--------------------------------------------------------------------
!
!..loop over all nodes
!$omp parallel
!$omp do private(nfath) schedule(guided) &
!$omp reduction(+:nrdofH,nrdofE,nrdofV,nrdofQ)
   do nod=1,NRNODS
!
!  ...skip if a middle node
      select case(NODES(nod)%ntype)
         case(MDLB,MDLN,MDLP,MDLD) ; cycle
      end select
!
!  ...skip if the node has not been marked,
!     and deactivate if it is still active
      if (.not. visited(nod)) then
         if (Is_active(nod)) call deactivate(nod, nrdofH,nrdofE,nrdofV,nrdofQ)
         cycle
      endif
!
!  ...skip if active
      if (Is_active(nod)) cycle
!
      nfath=NODES(nod)%father
!
!  ...if father node has not been visited, then the node is unconstrained
!     thus activate the node
      if (.not. visited(nfath)) then
         call activate(nod, nrdofH,nrdofE,nrdofV,nrdofQ)
      endif
!
!..end of loop over nodes
   enddo
!$omp end do
!$omp end parallel
!
   call bitvisit_finalize
!
!..update global dof counters
   NRDOFSH = NRDOFSH + nrdofH
   NRDOFSE = NRDOFSE + nrdofE
   NRDOFSV = NRDOFSV + nrdofV
   NRDOFSQ = NRDOFSQ + nrdofQ
!
end subroutine refresh
