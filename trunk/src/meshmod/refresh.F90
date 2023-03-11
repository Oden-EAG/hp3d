!--------------------------------------------------------------------
!> @brief activate inactive nodes which are not constrained
!         but in the current mesh
!> @date Mar 2023
!--------------------------------------------------------------------
subroutine refresh
!
   use data_structure3D
   use par_mesh
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
   if ((.not. DISTRIBUTED) .or. HOST_MESH) goto 99
!
!-----------------------------------------------------------------------
!  Remark: We must make sure that unconstrained active nodes are only
!          marked by a subdomain if they are part of a modified element
!          of the subdomain.
!          Otherwise, this creates issues with determining node
!          ownership in solver interfaces without having to rerun
!          logic_nodes to obtain modified element nodes.
!          Without this fix, DOFs were also stored unnecessarily for
!          some active nodes outside of the subdomain.
!
!  An example where an edge was marked incorrectly:
!  (the same issue can happen for a face node)
!
!      The unconstrained, active       In href, edge0's kref already
!      edge0 is marked by subd 0       coincides with the desired href
!      b/c it is part of one of        but "break" always calls
!      its modified element lists;     "activate_sons" which sets subd
!      edge1 is constrained and        values for edge0's sons, incl.
!      only marked by subd 1 as        edge1, based on edge0's subd.
!      one of its local elem nodes.    So edge1 is marked by subd 0,
!                                      even though it is neither on
!                                      any modified element nor local
!                                      element node list of subd 0.
!
!                edge0                           edge0
!                  |            href ->            |
!                  |                               |
!                  v                               v
!      *-----------*-----------*       *-----------*-----------*
!      |           |           |       |           |           |
!      |  subd 0   |           |       |  subd 0   |  subd 1   |
!      |           |           |       |           |           |
!      *-----------*  subd 1   |       *-----------*-----------*
!      |           |           |       |           |           |
!      |  subd 1   | <-edge1   |       |  subd 1   |  subd 1   |
!      |           |           |       |           | <-edge1   |
!      *-----------*-----------*       *-----------*-----------*
!
!-----------------------------------------------------------------------
!
!$OMP PARALLEL PRIVATE(iel,mdle,nod,subd)
!
!..1. Reset subdomain values
!$OMP DO
   do nod=1,NRNODS
      call get_subd(nod, subd)
      if ((subd.eq.RANK) .and. (.not.Is_middle(nod))) then
         call set_subd(nod,-1)
      endif
   enddo
!$OMP END DO
!
!..2. Set subdomains for all (unconstrained) nodes within subdomain
!$OMP DO
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      call set_subd_elem(mdle)
   enddo
!$OMP END DO
!
!..3. Delete degrees of freedom for NODES outside of subdomain
!$OMP DO
   do nod=1,NRNODS
      call get_subd(nod, subd)
      if ((subd.ne.RANK) .and. Is_active(nod)) then
#if DEBUG_MODE
!     ...print nodes that had previously been incorrectly marked
         if ((iprint.eq.2) .and. associated(NODES(nod)%dof)) then
            write(*,2345) RANK,nod
            2345 format('[',I2,'] refresh: dealloc nod = ', I5)
         endif
#endif
!     ...delete solution degrees of freedom
         call dealloc_nod_dof(nod)
      endif
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
#if DEBUG_MODE
   call par_verify
#endif
!
  99 continue
!
end subroutine refresh
