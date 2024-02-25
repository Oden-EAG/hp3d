!--------------------------------------------------------------------
!> @brief activate inactive nodes which are not constrained
!         but in the current mesh
!> @date Mar 2023
!--------------------------------------------------------------------
subroutine refresh
!
   use data_structure3D
   use bitvisit
   use par_mesh
   use mpi_param, only: RANK,ROOT
!
   implicit none
!
   integer :: ntype
   integer :: nodesl(27),norientl(27)
   integer :: i,iel,nod,nfath,nvef,mdle,subd
   integer :: nrdofH,nrdofE,nrdofV,nrdofQ
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!--------------------------------------------------------------------
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,*) 'refresh: Begin'
   endif
#endif
!
   call update_ELEM_ORDER
   call bitvisit_init(NRNODS)
!
!..set thread local dof counters
   nrdofH=0; nrdofE=0; nrdofV=0; nrdofQ=0
!
!--------------------------------------------------------------------
! Step 1 : raise visitation flag for vertex, edge and face nodes of |
!          all active elements                                      |
!--------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(mdle,subd,nodesl,norientl,ntype,nvef) &
!$OMP    SCHEDULE(GUIDED)
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
   enddo
!$OMP END DO
!$OMP END PARALLEL
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
!$OMP PARALLEL
!$OMP DO PRIVATE(nfath) SCHEDULE(GUIDED) &
!$OMP    REDUCTION(+:nrdofH,nrdofE,nrdofV,nrdofQ)
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
!$OMP END DO
!$OMP END PARALLEL
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
!
!
!
!-----------------------------------------------------------------------
!> @brief Ensures that DOFs are allocated correctly within subdomains
!> @details This routines must be called after refinements are done and
!!          one-irregularity of the mesh was enforced by close_mesh.
!> @date Mar 2023
!-----------------------------------------------------------------------
!  Remarks
!-----------------------------------------------------------------------
!       1. We must make sure that unconstrained active nodes are only
!          marked by a subdomain if they are part of a modified element
!          of the subdomain.
!          Otherwise, this creates issues with determining node
!          ownership in solver interfaces without having to rerun
!          logic_nodes to obtain modified element nodes.
!          Without this fix, DOFs were also stored unnecessarily for
!          some active nodes outside of the subdomain.
!
!  An example where face sons were marked incorrectly:
!
!      The unconstrained, active       In href, face0's kref already
!      face0 is marked by subd 0       coincides with the desired href
!      b/c it is part of one of        but "break" always calls
!      its modified element lists;     "activate_sons" which sets subd
!      face1 is constrained and        values for face0's sons, incl.
!      only marked by subd 1 as        face1, based on face0's subd.
!      one of its local elem nodes.    So face1 is marked by subd 0,
!                                      even though it is neither on
!                                      any modified element nor local
!                                      element node list of subd 0.
!
!                face0                           face0
!                  |            href ->            |
!                  |                               |
!                  v                               v
!      *-----------*-----------*       *-----------*-----------*
!      |           |           |       |           |           |
!      |  subd 0   |           |       |  subd 0   |  subd 1   |
!      |           |           |       |           |           |
!      *-----------*  subd 1   |       *-----------*-----------*
!      |           |           |       |           |           |
!      |  subd 1   | <-face1   |       |  subd 1   |  subd 1   |
!      |           |           |       |           | <-face1   |
!      *-----------*-----------*       *-----------*-----------*
!
!-----------------------------------------------------------------------
!       2. We must also make sure that unconstrained active nodes marked
!          by a subdomain have in fact their DOFs allocated.
!          Background:
!          Suppose an edge that was previously refined, but has its sons
!          still constrained, is marked for refinement. In this case,
!          hp3D’s break routine does not activate the edge sons, unlike
!          when this happens for a face (which will activate its sons).
!          The reason for this behavior is when a face refinement is
!          requested and the existing face refinement coincides with the
!          requested one, then both elements (a face is only attached
!          to two elements) are now refined, hence we can activate the
!          face sons. But an edge may be requested to be refined by more
!          than one element and its sons could still be constrained (as
!          long as not all the elements attached to the edge have been
!          refined). Thus, break does not call activate_sons when the
!          requested edge refinement coincides with the existing one.
!          Issue:
!          This can in some rare cases lead to an issue where in the
!          distributed mesh case, the DOFs are not correctly allocated
!          for sons of an edge when the constrained edge sons are
!          located on the subdomain boundary.
!-----------------------------------------------------------------------
subroutine distr_refresh()
!
   use data_structure3D
   use par_mesh
   use mpi_param, only: RANK
!
   implicit none
!
   integer :: iel,mdle,nod,subd
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
   if ((.not. DISTRIBUTED) .or. HOST_MESH) return
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
!     and ensure that DOFs allocated for NODES inside subdomain
!$OMP DO
   do nod=1,NRNODS
      call get_subd(nod, subd)
      if ((subd.ne.RANK) .and. Is_active(nod)) then
!     ...delete solution degrees of freedom
         call dealloc_nod_dof(nod)
!
#if HP3D_DEBUG
!     ...print nodes that had previously been incorrectly marked
         if ((iprint.eq.1) .and. associated(NODES(nod)%dof)) then
            !$OMP CRITICAL
            write(*,2345) RANK,nod
            2345 format('[',I2,'] distr_refresh: dealloc nod = ', I5)
            !$OMP END CRITICAL
         endif
#endif
!
      elseif ((subd.eq.RANK) .and. Is_active(nod)) then
         if (.not.associated(NODES(nod)%dof)) then
!        ...allocate solution degrees of freedom (see Remark 2)
            call alloc_nod_dof(nod)
!
#if HP3D_DEBUG
!        ...print active subd nodes that had previously not been allocated
            if (iprint.eq.1) then
               !$OMP CRITICAL
               write(*,2346) RANK,nod
               2346 format('[',I2,'] distr_refresh:   alloc nod = ', I5)
               !$OMP END CRITICAL
            endif
#endif
         endif
      endif
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
end subroutine distr_refresh
