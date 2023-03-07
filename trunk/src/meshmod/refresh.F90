!--------------------------------------------------------------------
!> @brief activate inactive nodes which are not constrained
!          but in the current mesh
!> @date Oct 2019
!--------------------------------------------------------------------
subroutine refresh
!
   use data_structure3D
   use par_mesh  , only: DISTRIBUTED,set_subd_elem
   use mpi_param , only: RANK
!
   implicit none
!
   integer :: ntype
   integer :: nodesl(27),norientl(27)
   integer :: i,iel,nod,nfath,mdle,subd
   integer :: nrdofH,nrdofE,nrdofV,nrdofQ
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
!
!..set thread local dof counters
   nrdofH=0; nrdofE=0; nrdofV=0; nrdofQ=0
!
!..reset visitation flags for all nodes
!$OMP PARALLEL
!$OMP DO
   do i=1,NRNODS
      NODES(i)%visit = 0
   enddo
!$OMP END DO
!
!--------------------------------------------------------------------
! Step 1 : raise visitation flag for vertex, edge and face nodes of |
!          all active elements                                      |
!--------------------------------------------------------------------
!$OMP DO PRIVATE(mdle,subd,nodesl,norientl,ntype,i)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      call elem_nodes(mdle, nodesl,norientl)
      ntype=NODES(mdle)%ntype
      do i=1,nvert(ntype)+nedge(ntype)+nface(ntype)
         NODES(nodesl(i))%visit=1
!     ...if node is visited by an element within my subdomain,
!        add node to my subdomain (need its dofs). this flag will
!        indicate that dofs must be allocated in activation.
         if (DISTRIBUTED .and. (subd.eq.RANK)) then
            call set_subd(nodesl(i),subd)
         endif
      enddo
!
   enddo
!$OMP END DO
!
!--------------------------------------------------------------------
! Step 2: activate all inactive edge and face nodes whose father    |
!         node was not visited (i.e., activate inactive edge and    |
!         face nodes that are unconstrained)                        |
!--------------------------------------------------------------------
!
!..loop over all nodes
!$OMP DO PRIVATE(nfath) SCHEDULE(DYNAMIC) &
!$OMP REDUCTION(+:nrdofH,nrdofE,nrdofV,nrdofQ)
   do nod=1,NRNODS
!
!  ...skip if a middle node
      select case(NODES(nod)%ntype)
         case(MDLB,MDLN,MDLP,MDLD) ; cycle
      end select
!
!  ...skip if the node has not been marked,
!     and deactivate if it is still active
      if (NODES(nod)%visit.eq.0) then
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
      if (NODES(nfath)%visit.eq.0) then
         call activate(nod, nrdofH,nrdofE,nrdofV,nrdofQ)
      endif
!
!..end of loop over nodes
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
!..update global dof counters
   NRDOFSH = NRDOFSH + nrdofH
   NRDOFSE = NRDOFSE + nrdofE
   NRDOFSV = NRDOFSV + nrdofV
   NRDOFSQ = NRDOFSQ + nrdofQ
!
end subroutine refresh
