!
!-----------------------------------------------------------------------
!
!    routine   - par_ghost
!
!-----------------------------------------------------------------------
!
!    last      - Mar 2023
!
!    purpose   - Gets all elements sharing (modified) nodes with subdomain;
!                stores in ELEM_GHOST and ELEM_INTERF
!
!                ELEM_GHOST : ELEM_SUBD plus elements touching my subd
!                ELEM_INTERF: mdle nodes touching subd interface
!
!                both are stored in natural order of elements
!
!-----------------------------------------------------------------------
module par_ghost

implicit none
!
contains


!----------------------------------------------------------------------
!> @brief computes ghosted subdomain and interface elements in parallel
!> @date Feb 2023
!----------------------------------------------------------------------
   subroutine get_ghost_subd
!
   use data_structure3D
   use bitvisit
   use par_mesh, only: DISTRIBUTED, HOST_MESH
   use mpi_wrapper
!
   implicit none
!
   integer :: nodesl(27), norientl(27), nodm(MAXNODM)
   integer :: nverts(NRELES_SUBD), verts(16,NRELES_SUBD)
   integer :: my_elems(NRELES_SUBD)
   integer :: nrelem_procs(NUM_PROCS), displs(NUM_PROCS)
!
   integer :: iel, mdle, subd, nrnodm, nrvert, nrelem
   integer :: i, nod, nrelem_interf, ierr, n, ivert, off, sum
   integer :: ig, il, mdleI
!
   logical :: shared
!
   integer, allocatable :: vert_nprocs(:), interf_elems(:), elem_verts(:)
   integer, allocatable :: elem_nrvert(:), subd2interf(:), offsets(:)
!
!-----------------------------------------------------------------------
!
   if (allocated(ELEM_GHOST)) deallocate(ELEM_GHOST)
   if (allocated(ELEM_INTERF)) deallocate(ELEM_INTERF)
   NRELES_GHOST = 0
   NRELES_INTERF = 0
   
!
   if (.not. DISTRIBUTED .or. HOST_MESH) then
      allocate(ELEM_GHOST(NRELES)); ELEM_GHOST = ELEM_ORDER
      NRELES_GHOST = NRELES
      return
!  ...interface not defined on serial mesh
   endif
!
!-----------------------------------------------------------------------
! Step 1: Get list of all vertex nodes on my subd
!-----------------------------------------------------------------------
!
   call bitvisit_init(NRNODS)
!
!$omp parallel do default(shared)                              &
!$omp    private(iel,mdle,i,nodm,nrnodm,nodesl,norientl,nod)   &
!$omp    schedule(guided)
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!
!  ...get unconstrained nodes
      call get_connect_info(mdle, nodesl,norientl)
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
      nverts(iel) = 0
      verts(:,iel) = 0
!
!  ...Loop through active nodes
      do i = 1,nrnodm
         nod = nodm(i)
!
!     ...visit vertex nodes and store elem to vertex map
         if (NODES(nod)%ntype.eq.VERT) then
            call visit(nod)
            nverts(iel) = nverts(iel) + 1
            verts(nverts(iel),iel) = nod
         endif
      enddo
   enddo
!$omp end parallel do
!
   call reduce_visit
!
!..Make map from nodes to vertices
   call reset_visit
!
   nrvert = 0;
   do nod = 1,NRNODS
      if (visited(nod)) then
         nrvert = nrvert + 1
         NODES(nod)%visit = nrvert
      endif
   enddo
!
!-----------------------------------------------------------------------
! Step 2: Get list of only shared vertex functions
!-----------------------------------------------------------------------
!
!..Initialize array to count number of procs touching a vertex
   allocate(vert_nprocs(nrvert)); vert_nprocs = 0
!
!$omp parallel do default(shared)        &
!$omp    private(iel,mdle,i,nod,ivert)   &
!$omp    schedule(static)
!..Loop again, marking all verts
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!
!  ...Loop through element vertices
      do i = 1,nverts(iel)
         nod = verts(i,iel)
         ivert = NODES(nod)%visit
         vert_nprocs(ivert) = 1
      enddo
   enddo
!$omp end parallel do
!
   call MPI_ALLREDUCE(MPI_IN_PLACE,vert_nprocs,nrvert,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, ierr)
!
!..Renumber verts, only counting interface verts
   nrvert = 0;
   do nod = 1,NRNODS
      if (NODES(nod)%visit .gt. 0) then
         ivert = NODES(nod)%visit
         shared = vert_nprocs(ivert).gt.1
         nrvert = nrvert + merge(1,0,shared)
         NODES(nod)%visit = merge(nrvert,0,shared)
      endif
   enddo
   deallocate(vert_nprocs)
!
!-----------------------------------------------------------------------
! Step 3: Get list of subdomain elements touching interface
!-----------------------------------------------------------------------
!
   nrelem = 0;
!
!$omp parallel do default(shared)   &
!$omp    private(iel,mdle,i,nod)    &
!$omp    schedule(static)
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!
!  ...Loop through element vertices
      do i = 1,nverts(iel)
         nod = verts(i,iel)
!
         if (NODES(nod)%visit.ge.1) then
!$omp critical
            nrelem = nrelem + 1
            my_elems(nrelem) = mdle
!$omp end critical
            exit
         endif
      enddo
   enddo
!$omp end parallel do
!
!-----------------------------------------------------------------------
! Step 4: Collect global list of elements touching interface
!-----------------------------------------------------------------------
!
   call MPI_Allgather(nrelem,1,MPI_INTEGER,nrelem_procs,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)
!
   displs = 0
   nrelem_interf = 0
   do i = 1,NUM_PROCS
      displs(i) = nrelem_interf
      nrelem_interf = nrelem_interf + nrelem_procs(i)
   enddo
!
   allocate(interf_elems(nrelem_interf))
!
   call MPI_AllgatherV(my_elems,nrelem,MPI_INTEGER, interf_elems,nrelem_procs,displs,MPI_INTEGER, MPI_COMM_WORLD,ierr)
!
!..Mark interface elems
   call bitvisit_reset
   do iel=1,nrelem_interf
      call visit(interf_elems(iel))
   enddo
!
!..Put interface elems in natural elem order
   nrelem = 0;
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
!
      if (visited(mdle)) then
         nrelem = nrelem + 1
         interf_elems(nrelem) = mdle
      endif
   enddo
!
   if (nrelem .ne. nrelem_interf) then
      write(*,*) 'Error getting ghosted subdomain.'
      stop 5
   endif
!
!-----------------------------------------------------------------------
! Step 5: Get interface elem to vert connectivity (globally)
!-----------------------------------------------------------------------
!..Note: This whole step is to avoid computing nodes for elems not
!        in subdomain. We could simply compute vertex nodes for interface
!        elems for step 6 below but for large meshes that is still
!        expensive (~10 seconds for 1 million elements on 16 ranks).
!        The approach here instead collects info so is scalable and is
!        faster even on relatively few ranks (<1 second in same example).
!
   allocate(elem_nrvert(nrelem_interf)); elem_nrvert = 0
   allocate(subd2interf(NRELES_SUBD)); subd2interf = 0
!
!..Construct subd to interface element map
!..Find subdomain elems in interface elem array
   i = 0;
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!
!  ...Skip if not interface elem
      if (.not. visited(mdle)) cycle
!
!  ...Interface mdle
      i = i + 1
      mdleI = interf_elems(i)
!
!  ...Iterate until found (subd and interf follow same ordering)
      do while (mdle .ne. mdleI)
         i = i + 1
         mdleI = interf_elems(i)
      enddo
      subd2interf(iel) = i
   enddo
!
!$omp parallel do default(shared)         &
!$omp    private(iel,mdle,i,nod,ig)       &
!$omp    schedule(static)
!..Loop again, marking all verts
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!
!  ...Skip if not interface elem
      if (.not. visited(mdle)) cycle
!
      ig = subd2interf(iel)
!
!  ...Loop through element vertices
      do i = 1,nverts(iel)
         nod = verts(i,iel)
         if (NODES(nod)%visit .gt. 0) then
            elem_nrvert(ig) = elem_nrvert(ig) + 1
         endif
      enddo
   enddo
!$omp end parallel do
!
   call MPI_ALLREDUCE(MPI_IN_PLACE,elem_nrvert,nrelem_interf,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, ierr)
!
   allocate(offsets(nrelem_interf)); offsets = 0;
!
   sum = 0
!
!NOTE: Need Intel compiler 19.0+ to use simd for scan
!!$omp parallel do simd reduction(inscan, +:sum)
!..Exclusive scan to get offsets
   do iel = 1,nrelem_interf
      offsets(iel) = sum
!!$omp scan exclusive(sum)
      sum = sum + elem_nrvert(iel)
   enddo
!
!..Interface elem to vert map
   allocate(elem_verts(sum)); elem_verts = 0;
!
!..Loop again, loading vertices into array
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!
!  ...Skip if not interface elem
      if (.not. visited(mdle)) cycle
!
      ig = subd2interf(iel)
      off = offsets(ig)
      il = 0
!
!  ...Load element vertex map
      do i = 1,nverts(iel)
         nod = verts(i,iel)
         if (NODES(nod)%visit .gt. 0) then
            il = il + 1
            elem_verts(off + il) = nod
         endif
      enddo
   enddo
!
   call MPI_ALLREDUCE(MPI_IN_PLACE,elem_verts,sum,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, ierr)
!
!-----------------------------------------------------------------------
! Step 6: Add only interface elements touching my subdomain
!-----------------------------------------------------------------------
!
   NRELES_INTERF = 0
   NRELES_GHOST = 0
!
   allocate(ELEM_INTERF(nrelem_interf))
   n = min(NRELES,NRELES_SUBD + nrelem_interf)
   allocate(ELEM_GHOST(n))
!
!..Visit subdomain vertices
   call bitvisit_reset
   do iel = 1,NRELES_SUBD
      do i=1,nverts(iel)
         nod = verts(i,iel)
         call visit(nod)
      enddo
   enddo
!
!..Add subd elements to ghost
   NRELES_GHOST = NRELES_SUBD
   ELEM_GHOST(1:NRELES_SUBD) = ELEM_SUBD(1:NRELES_SUBD)
!
!$omp parallel do default(shared)                                 &
!$omp    private(iel,mdle,subd,i,nodm,nrnodm,nodesl,norientl,nod, &
!$omp            off,n)                                           &
!$omp    schedule(guided)
!..Loop interface elems, get only those touching my subd
   do iel=1,nrelem_interf
      mdle = interf_elems(iel)
      subd = NODES(mdle)%subd
!
!  ...In subd and on interface, add to interface
      if (subd.eq.RANK) then
         call add_elem_interf(mdle)
         cycle
      endif
!
!  ...get elem vertex nodes
      off = offsets(iel)
      n = elem_nrvert(iel)
!
!  ...Loop active nodes
      do i = off+1,off+n
         nod = elem_verts(i)
!
!     ...Not in subd but on interface, add to interface and ghost
         if (visited(nod)) then
            call add_elem_interf(mdle)
            call add_elem_ghost(mdle)
            exit
         endif
      enddo
   enddo
!$omp end parallel do
!
   deallocate(interf_elems,elem_verts,elem_nrvert,subd2interf,offsets)
!
!..Reorder ghost elements to follow natural order of elements
   call bitvisit_reset
!
   do iel = 1,NRELES_GHOST
      mdle = ELEM_GHOST(iel)
      call visit(mdle)
   enddo
!
   nrelem = 0
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      if (visited(mdle)) then
         nrelem = nrelem + 1
         ELEM_GHOST(nrelem) = mdle
      endif
   enddo
!
   call bitvisit_reset
!
   do iel = 1,NRELES_INTERF
      mdle = ELEM_INTERF(iel)
      call visit(mdle)
   enddo
!
   nrelem = 0
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      if (visited(mdle)) then
         nrelem = nrelem + 1
         ELEM_INTERF(nrelem) = mdle
      endif
   enddo
!
   call bitvisit_finalize
!
#if HP3D_DEBUG
   call check_ghost_elems
   call check_interf_elems
#endif
!
   contains
!
!  ...Subroutines defined to have coherent critical between different calls
      subroutine add_elem_interf(mdle)
         integer, intent(in) :: mdle
!$omp critical
         NRELES_INTERF = NRELES_INTERF + 1
         ELEM_INTERF(NRELES_INTERF) = mdle
!$omp end critical
      end subroutine add_elem_interf

      subroutine add_elem_ghost(mdle)
         integer, intent(in) :: mdle
!$omp critical
         NRELES_GHOST = NRELES_GHOST + 1
         ELEM_GHOST(NRELES_GHOST) = mdle
!$omp end critical
      end subroutine add_elem_ghost
!
   end subroutine get_ghost_subd




!----------------------------------------------------------------------
!> @brief checks interface elements computed by get_ghost_subd
!> @date Feb 2023
!----------------------------------------------------------------------
   subroutine check_interf_elems
!
   use data_structure3D
   use bitvisit
   use mpi_param
   use par_mesh, only: DISTRIBUTED
!
   implicit none
!
   integer :: nodesl(27), norientl(27), nodm(MAXNODM)
!
   integer :: iel, mdle, subd, nrnodm
   integer :: i, nod, nr_interf_elems
   integer, allocatable :: interf_elems(:)
!
!-----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) return
!
   allocate(interf_elems(NRELES));
   nr_interf_elems = 0
!
   call reset_visit
!
!$omp parallel do default(shared)                              &
!$omp    private(iel,mdle,i,nodm,nrnodm,nodesl,norientl,nod)   &
!$omp    schedule(guided)
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!
!  ...get unconstrained nodes
      call get_connect_info(mdle, nodesl,norientl)
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
!  ...Loop through active nodes
      do i = 1,nrnodm
         nod = nodm(i)
         NODES(nod)%visit = 1
      enddo
   enddo
!$omp end parallel do
!
!$omp parallel do default(shared)                                    &
!$omp    private(iel,mdle,subd,i,nodm,nrnodm,nodesl,norientl,nod)    &
!$omp    schedule(guided)
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      subd = NODES(mdle)%subd
!
      if (subd.eq.RANK) cycle
!
!  ...get unconstrained nodes
      call get_connect_info(mdle, nodesl,norientl)
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
!  ...Loop through active nodes
      do i = 1,nrnodm
         nod = nodm(i)
!$omp critical
         NODES(nod)%visit = merge(0,2,NODES(nod)%visit.eq.0)
!$omp end critical
      enddo
   enddo
!$omp end parallel do
!
!$omp parallel do default(shared)                                    &
!$omp    private(iel,mdle,i,nodm,nrnodm,nodesl,norientl,nod)    &
!$omp    schedule(guided)
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
!
!  ...get unconstrained nodes
      call get_connect_info(mdle, nodesl,norientl)
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
!  ...Loop through active nodes
      do i = 1,nrnodm
         nod = nodm(i)
         if (NODES(nod)%visit .gt. 1) then
!$omp critical
            nr_interf_elems = nr_interf_elems + 1
            interf_elems(nr_interf_elems) = mdle
!$omp end critical
            exit
         endif
      enddo
   enddo
!$omp end parallel do
!
   call bitvisit_init(NRNODS)
!
!..reorder ghost elems
   do iel = 1, nr_interf_elems
      mdle = interf_elems(iel)
      call visit(mdle)
   enddo
!
   i = 0;
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      if (visited(mdle)) then
         i = i + 1
         interf_elems(i) = mdle
      endif
   enddo
!
   if (nr_interf_elems .ne. NRELES_INTERF) then
      write(*,*) 'check_interf_elems: Number of interface elements incorrect on Rank', RANK, ': ', NRELES_INTERF, ' vs ', nr_interf_elems, 'expected.'
      stop 1
   endif
!
!..compare to parallel computation
   do iel = 1,NRELES_INTERF
      if (interf_elems(iel).ne.ELEM_INTERF(iel)) then
         write(*,*) 'check_interf_elems: Unexpected ghost element on Rank', RANK, ': ',ELEM_INTERF(iel), ' vs: ', interf_elems(iel)
         stop 1
      endif
   enddo
!
   deallocate(interf_elems)
!
   call bitvisit_finalize
!
   end subroutine check_interf_elems




!----------------------------------------------------------------------
!> @brief checks ghost elements computed by get_ghost_subd
!> @date Feb 2023
!----------------------------------------------------------------------
   subroutine check_ghost_elems
!
   use data_structure3D
   use bitvisit
   use mpi_param
   use par_mesh, only: DISTRIBUTED
!
   implicit none
!
   integer :: nodesl(27), norientl(27), nodm(MAXNODM)
!
   integer :: iel, mdle, subd, nrnodm
   integer :: i, nod, nr_ghost_elems
   integer, allocatable :: ghost_elems(:)
!
!-----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) return
!
   allocate(ghost_elems(NRELES)); ghost_elems(1:NRELES_SUBD) = ELEM_SUBD(1:NRELES_SUBD)
   nr_ghost_elems = NRELES_SUBD
!
   call bitvisit_init(NRNODS)
!
!$omp parallel do default(shared)                              &
!$omp    private(iel,mdle,i,nodm,nrnodm,nodesl,norientl,nod)   &
!$omp    schedule(guided)
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!
!  ...get unconstrained nodes
      call get_connect_info(mdle, nodesl,norientl)
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
!  ...Loop through active nodes
      do i = 1,nrnodm
         nod = nodm(i)
         call visit(nod)
      enddo
   enddo
!$omp end parallel do
!
!$omp parallel do default(shared)                                    &
!$omp    private(iel,mdle,subd,i,nodm,nrnodm,nodesl,norientl,nod)    &
!$omp    schedule(guided)
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      subd = NODES(mdle)%subd
!
      if (subd.eq.RANK) cycle
!
!  ...get unconstrained nodes
      call get_connect_info(mdle, nodesl,norientl)
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
!  ...Loop through active nodes
      do i = 1,nrnodm
         nod = nodm(i)
!
         if (visited(nod)) then
!$omp critical
            nr_ghost_elems = nr_ghost_elems + 1
            ghost_elems(nr_ghost_elems) = mdle
!$omp end critical
            exit
         endif
      enddo
   enddo
!$omp end parallel do
!
   call bitvisit_reset
!
!..reorder ghost elems
   do iel = 1, nr_ghost_elems
      mdle = ghost_elems(iel)
      call visit(mdle)
   enddo
!
   i = 0;
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      if (visited(mdle)) then
         i = i + 1
         ghost_elems(i) = mdle
      endif
   enddo
!
   if (nr_ghost_elems .ne. NRELES_GHOST) then
      write(*,*) 'check_ghost_elems: Number of elements in ghost subdomain incorrect on Rank', RANK, ': ', NRELES_GHOST, ' vs ', nr_ghost_elems, 'expected.'
      stop 1
   endif
!
!..compare to parallel computation
   do iel = 1,NRELES_GHOST
      if (ghost_elems(iel).ne.ELEM_GHOST(iel)) then
         write(*,*) 'check_ghost_elems: Unexpected ghost element on Rank', RANK, ': ',ELEM_GHOST(iel), ' vs: ', ghost_elems(iel)
         stop 1
      endif
   enddo
!
   deallocate(ghost_elems)
!
   call bitvisit_finalize
!
   end subroutine check_ghost_elems


end module par_ghost


