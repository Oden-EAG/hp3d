!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!     module:              par_mesh
!
!     last modified:       Oct 2019
!
!     purpose:             provides functionality for distributing
!                          degrees of freedom in the FE mesh
!
!----------------------------------------------------------------------
module par_mesh
!
   use data_structure3D
   use parameters,     only: NRCOMS
   use mpi_param,      only: RANK,ROOT,NUM_PROCS
   use MPI,            only: MPI_COMM_WORLD,MPI_STATUS_IGNORE, &
                             MPI_SUCCESS,MPI_COMPLEX16,MPI_REAL8
   use zoltan_wrapper, only: ZOLTAN_LB,zoltan_w_partition
!
   implicit none
!
!..F if current mesh is not distributed (ALL DOFs present on ALL procs)
!..T if current mesh is distributed (DOFs are distributed across procs)
   logical, save :: DISTRIBUTED = .false.
!
!..T indicates that dofs should be exchanged (send/rcv) in routine
   logical, save :: EXCHANGE_DOF = .true.
!
!..T indicates that the ROOT proc holds the entire mesh
   logical, save :: HOST_MESH = .true.
!
   contains
!
!----------------------------------------------------------------------
!     routine:    distr_mesh
!     purpose:    (re-)distribute the current mesh to MPI processes
!----------------------------------------------------------------------
subroutine distr_mesh()
!
!..MPI variables
   integer :: ierr, base, remainder
   integer :: tag, count, src, dest
   VTYPE, allocatable :: buf(:)
!
!..auxiliary variables
   integer :: subd_next(NRELES), nodm(MAXNODM), nodesl(27)
   integer :: iel, inod, iproc, mdle, nod, subd, subd_size
   integer :: nrnodm, nrdof_nod, proc, el
   integer :: iprint = 0
!
   real(8) :: MPI_Wtime,start_time,end_time
!
   if (iprint .eq. 1) then
      write(6,100) 'start distr_mesh, DISTRIBUTED = ', DISTRIBUTED
   endif
   100 format(A,L2)
!
   if (NUM_PROCS.eq.1) return
!
!..1. Determine new partition
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   if ((ZOLTAN_LB .eq. 0) .or. (.not. DISTRIBUTED)) then
!
      base = NRELES/NUM_PROCS
      remainder = mod(NRELES,NUM_PROCS)
!
!  ...element offset
      iel=0
!
!  ...first remainder procs get base+1, rest get base
      do iproc = 0,NUM_PROCS-1
         if (iproc .lt. remainder) then
            subd_size = base + 1
         else
            subd_size = base
         endif
         subd_next(iel+1:iel+subd_size) = iproc
         iel = iel + subd_size
      enddo
   elseif (ZOLTAN_LB .eq. 7) then
      if (RANK .eq. ROOT) write(*,*) 'calling (re-)partition_fiber...'
      !call partition_fiber(subd_next)
      call repartition_fiber(subd_next)
   else
      call zoltan_w_partition(subd_next)
   endif
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
   if(RANK .eq. ROOT) write(*,110) end_time - start_time
   110 format(' partition time: ',f12.5,' seconds')
!
!..2. Reset visit flags for all nodes to 0
   call reset_visit
!
!..3. Exchange degrees of freedom with other MPI processes according to partition
!     (note: this step can be skipped in the initial mesh distribution)
   if (.not. DISTRIBUTED) goto 50
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
!     if mdle node current subdomain is not equal its new subdomain, and
!     if current or new subdomain are my subdomain, then
      if (subd .eq. subd_next(iel)) cycle
      if ((subd .ne. RANK) .and. (subd_next(iel) .ne. RANK)) cycle
!  ...3a. get list of nodes associated with mdle node
      call get_elem_nodes(mdle, nodesl,nodm,nrnodm)
!  ...3b. iterate over list of nodes
      do inod=1,nrnodm
         nod = nodm(inod)
         if (.not. EXCHANGE_DOF) then
            if (subd_next(iel) .eq. RANK) call alloc_nod_dof(nod)
            cycle
         endif
!     ...3c. Calculate nodal degrees of freedom, and allocate buffer
         call get_dof_buf_size(nod, nrdof_nod)
         if (nrdof_nod .eq. 0) cycle
         allocate(buf(nrdof_nod))
         buf = ZERO
!     ...3d. if current subdomain is my subdomain, send data
         if (subd .eq. RANK) then
!        ...3e. pack buffer with DOF data
            call pack_dof_buf(nod,nrdof_nod, buf)
!        ...3f. send buffer to new subdomain
            if (iprint .eq. 1) then
               write(6,130) '[', RANK, ']: ', &
                  'Sending data to [',subd_next(iel),'], nod = ',nod
            endif
            count = nrdof_nod; dest = subd_next(iel); tag = mod(nod,200000)
            call MPI_SEND(buf,count,MPI_VTYPE,dest,tag,MPI_COMM_WORLD,ierr)
            if (ierr .ne. MPI_SUCCESS) then
               write(6,*) 'MPI_SEND failed. stop.'
               stop
            endif
!     ...3d. if new subdomain is my subdomain, receive data
         else if (subd_next(iel) .eq. RANK) then
!        ...3e. receive buffer from old subdomain
            if (iprint .eq. 1) then
               write(6,130) '[', RANK, ']: ',   &
                  'Receiving data from [',subd,'], nod = ',nod
            endif
            count = nrdof_nod; src = subd; tag = mod(nod,200000)
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            if (ierr .ne. MPI_SUCCESS) then
               write(6,*) 'MPI_RECV failed. stop.'
               stop
            endif
!        ...3f. unpack DOF data from buffer
            call alloc_nod_dof(nod)
            call unpack_dof_buf(nod,nrdof_nod,buf)
         endif
         deallocate(buf)
!
      enddo
   enddo
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
   if(RANK .eq. ROOT) write(*,120) end_time - start_time
   120 format(' migration time: ',f12.5,' seconds')
   130 format(A,I2,A,A,I4,A,I6)
!
   50 continue
!
!$OMP PARALLEL PRIVATE(iel,mdle,nod,subd)
!
!..4. Reset subdomain values for all nodes
!$OMP DO
   do nod=1,NRNODS
      call set_subd(nod,-1)
   enddo
!$OMP END DO
!
!..5. Set new subdomains for middle nodes (elements) everywhere
!$OMP DO
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call set_subd(mdle,subd_next(iel))
   enddo
!$OMP END DO
!
!..6. Set subdomain values for all (unconstrained) nodes within subdomain
!$OMP DO SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      if (subd_next(iel) .eq. RANK) then
         mdle = ELEM_ORDER(iel)
         call set_subd_elem(mdle)
      endif
   enddo
!$OMP END DO
!
!..7. Delete degrees of freedom for NODES outside of subdomain
!$OMP DO
   do nod=1,NRNODS
      call get_subd(nod, subd)
      if (subd .ne. RANK .and. Is_active(nod)) then
!     ...delete solution degrees of freedom
         call dealloc_nod_dof(nod)
      endif
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
   if (NUM_PROCS > 1) HOST_MESH = .false.
   DISTRIBUTED = .true.
   call update_ELEM_ORDER
!
   if ((.not. EXCHANGE_DOF) .and. (.not. HOST_MESH)) then
      call update_gdof
      call update_Ddof
   endif
!
   if (iprint .eq. 1) then
      write(6,100) 'end   distr_mesh, DISTRIBUTED = .true.'
   endif
!
end subroutine distr_mesh
!
!----------------------------------------------------------------------
!     routine:    get_elem_nodes
!     purpose:    get (unconstrained) nodes associated with an element
!----------------------------------------------------------------------
subroutine get_elem_nodes(Mdle, Nodesl,Nodm,Nrnodm)
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Nodesl(27)
   integer, intent(out) :: Nodm(MAXNODM)
   integer, intent(out) :: Nrnodm
   integer :: norientl(27)
   Nodesl(1:27) = 0
   Nodm(1:MAXNODM) = 0; Nrnodm = 0
   call get_connect_info(Mdle, Nodesl,norientl)
   call logic_nodes(Mdle,Nodesl, Nodm,Nrnodm)
end subroutine get_elem_nodes
!
!----------------------------------------------------------------------
!     routine:    set_subd_elem
!     purpose:    set subd values for (unconstrained) nodes associated
!                 with an element
!----------------------------------------------------------------------
subroutine set_subd_elem(Mdle)
   integer, intent(in) :: Mdle
   integer :: nodm(MAXNODM),nodesl(27)
   integer :: i,nrnodm,subd,ntype
   call get_elem_nodes(Mdle, nodesl,nodm,nrnodm)
   call get_subd(Mdle, subd)
   ntype = NODES(Mdle)%ntype
   do i=1,nvert(ntype)+nedge(ntype)+nface(ntype)
      call set_subd(nodesl(i),subd)
   enddo
   do i=1,nrnodm
      call set_subd(nodm(i),subd)
   enddo
end subroutine set_subd_elem
!
!----------------------------------------------------------------------
!     routine:    alloc_nod_dof
!     purpose:    allocate all dofs associated with a node
!----------------------------------------------------------------------
subroutine alloc_nod_dof(Nod)
   integer, intent(in) :: Nod
   integer :: ndofH,ndofE,ndofV,ndofQ
   integer :: nvarH,nvarE,nvarV,nvarQ
   integer :: icase
!..calculate ndof,nvar for this node
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
   icase = NODES(Nod)%case
   nvarH = NREQNH(icase)*NRCOMS
   nvarE = NREQNE(icase)*NRCOMS
   nvarV = NREQNV(icase)*NRCOMS
   nvarQ = NREQNQ(icase)*NRCOMS
!..allocate dof data type
   if (.not. associated(NODES(Nod)%dof)) then
      allocate(NODES(Nod)%dof)
      nullify(NODES(Nod)%dof%coord)
      nullify(NODES(Nod)%dof%zdofH)
      nullify(NODES(Nod)%dof%zdofE)
      nullify(NODES(Nod)%dof%zdofV)
      nullify(NODES(Nod)%dof%zdofQ)
   endif
!..allocate geometry DOFs
   if (.not. associated(NODES(Nod)%dof%coord) .and. (ndofH .gt. 0)) then
      allocate(NODES(Nod)%dof%coord(NDIMEN, ndofH))
      NODES(Nod)%dof%coord = 0.d0
   endif
!..allocate H1 DOFS
   if (.not. associated(NODES(Nod)%dof%zdofH) .and. (ndofH .gt. 0)) then
      allocate(NODES(Nod)%dof%zdofH(nvarH, ndofH))
      NODES(Nod)%dof%zdofH = ZERO
   endif
!..allocate H(curl) DOFs
   if (.not. associated(NODES(Nod)%dof%zdofE) .and. (ndofE .gt. 0)) then
      allocate(NODES(Nod)%dof%zdofE(nvarE, ndofE))
      NODES(Nod)%dof%zdofE = ZERO
   endif
!..allocate H(div) DOFs
   if (.not. associated(NODES(Nod)%dof%zdofV) .and. (ndofV .gt. 0)) then
      allocate(NODES(Nod)%dof%zdofV(nvarV, ndofV))
      NODES(Nod)%dof%zdofV = ZERO
   endif
!..allocate L2 DOFs
   if (.not. associated(NODES(Nod)%dof%zdofQ) .and. (ndofQ .gt.  0)) then
      allocate(NODES(Nod)%dof%zdofQ(nvarQ, ndofQ))
      NODES(Nod)%dof%zdofQ = ZERO
   endif
end subroutine alloc_nod_dof
!
!----------------------------------------------------------------------
!     routine:    dealloc_nod_dof
!     purpose:    deallocate all dofs associated with a node
!----------------------------------------------------------------------
subroutine dealloc_nod_dof(Nod)
   integer, intent(in) :: Nod
   if (associated(NODES(Nod)%dof)) then
     if (associated(NODES(Nod)%dof%coord)) deallocate(NODES(Nod)%dof%coord)
     if (associated(NODES(Nod)%dof%zdofH)) deallocate(NODES(Nod)%dof%zdofH)
     if (associated(NODES(Nod)%dof%zdofE)) deallocate(NODES(Nod)%dof%zdofE)
     if (associated(NODES(Nod)%dof%zdofV)) deallocate(NODES(Nod)%dof%zdofV)
     if (associated(NODES(Nod)%dof%zdofQ)) deallocate(NODES(Nod)%dof%zdofQ)
     deallocate(NODES(Nod)%dof)
   endif
end subroutine dealloc_nod_dof
!
!----------------------------------------------------------------------
!     routine:    get_dof_buf_size
!     purpose:    calculate buffer size for message passing of nodal dofs
!----------------------------------------------------------------------
subroutine get_dof_buf_size(Nod, Nrdof_nod)
   integer, intent(in) :: Nod
   integer, intent(out) :: Nrdof_nod
   integer :: ndofH,ndofE,ndofV,ndofQ
   integer :: nvarH,nvarE,nvarV,nvarQ
   integer :: icase
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
   icase = NODES(Nod)%case
   nvarH = NREQNH(icase)*NRCOMS
   nvarE = NREQNE(icase)*NRCOMS
   nvarV = NREQNV(icase)*NRCOMS
   nvarQ = NREQNQ(icase)*NRCOMS
   Nrdof_nod = ndofH*NDIMEN + ndofH*nvarH + ndofE*nvarE + ndofV*nvarV + ndofQ*nvarQ
end subroutine get_dof_buf_size
!
!----------------------------------------------------------------------
!     routine:    pack_dof_buf
!     purpose:    pack nodal dofs into contiguous buffer for message passing
!----------------------------------------------------------------------
subroutine pack_dof_buf(Nod,Nrdof_nod, Buf)
   integer, intent(in)  :: Nod
   integer, intent(in)  :: Nrdof_nod
   VTYPE  , intent(out) :: Buf(Nrdof_nod)
   integer :: ndofH,ndofE,ndofV,ndofQ
   integer :: nvarH,nvarE,nvarV,nvarQ
   integer :: icase,ivar,j
!
   Buf(1:Nrdof_nod) = ZERO
!
!..calculate ndof,nvar for this node
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
   icase = NODES(Nod)%case
   nvarH = NREQNH(icase)*NRCOMS
   nvarE = NREQNE(icase)*NRCOMS
   nvarV = NREQNV(icase)*NRCOMS
   nvarQ = NREQNQ(icase)*NRCOMS
!
   j = 0
!
!..add geometry dof to buffer
   if(ndofH .gt. 0) then
      do ivar=1,NDIMEN
         Buf(j+1:j+ndofH) = NODES(Nod)%dof%coord(ivar,1:ndofH)
         j = j + ndofH
      enddo
   endif
!..add H1 dof to buffer
   if(nvarH .gt. 0 .and. ndofH .gt. 0) then
      do ivar=1,nvarH
         Buf(j+1:j+ndofH) = NODES(Nod)%dof%zdofH(ivar,1:ndofH)
         j = j + ndofH
      enddo
   endif
!..add H(curl) dof to buffer
   if(nvarE .gt. 0 .and. ndofE .gt. 0) then
      do ivar=1,nvarE
         Buf(j+1:j+ndofE) = NODES(Nod)%dof%zdofE(ivar,1:ndofE)
         j = j + ndofE
      enddo
   endif
!..add H(div) dof to buffer
   if(nvarV .gt. 0 .and. ndofV .gt. 0) then
      do ivar=1,nvarV
         Buf(j+1:j+ndofV) = NODES(Nod)%dof%zdofV(ivar,1:ndofV)
         j = j + ndofV
      enddo
   endif
!..add L2 dof to buffer
   if(nvarQ .gt. 0 .and. ndofQ .gt. 0) then
      do ivar=1,nvarQ
         Buf(j+1:j+ndofQ) = NODES(Nod)%dof%zdofQ(ivar,1:ndofQ)
         j = j + ndofQ
      enddo
   endif
!..verify j equals Nrdof_nod
   if( Nrdof_nod .ne. j ) then
      write(*,*) 'pack_dof_buf: Nrdof_nod .ne. j, stop.'
      stop
   endif
end subroutine pack_dof_buf
!
!----------------------------------------------------------------------
!     routine:    unpack_dof_buf
!     purpose:    unpack nodal dofs from contiguous buffer
!----------------------------------------------------------------------
subroutine unpack_dof_buf(Nod,Nrdof_nod,Buf)
   integer, intent(in) :: Nod
   integer, intent(in) :: Nrdof_nod
   VTYPE  , intent(in) :: Buf(Nrdof_nod)
   integer :: ndofH,ndofE,ndofV,ndofQ
   integer :: nvarH,nvarE,nvarV,nvarQ
   integer :: icase,ivar,j
!
!..calculate ndof,nvar for this node
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
!
   icase = NODES(Nod)%case
   nvarH = NREQNH(icase)*NRCOMS
   nvarE = NREQNE(icase)*NRCOMS
   nvarV = NREQNV(icase)*NRCOMS
   nvarQ = NREQNQ(icase)*NRCOMS
!
   j = 0
!..extract geometry dof from buffer
   if(ndofH .gt. 0) then
      do ivar=1,NDIMEN
         NODES(Nod)%dof%coord(ivar,1:ndofH) = real(Buf(j+1:j+ndofH))
         j = j + ndofH
      enddo
   endif
!..extract H1 dof from buffer
   if(nvarH .gt. 0 .and. ndofH .gt. 0) then
      do ivar=1,nvarH
         NODES(Nod)%dof%zdofH(ivar,1:ndofH) = Buf(j+1:j+ndofH)
         j = j + ndofH
      enddo
   endif
!..extract H(curl) dof from buffer
   if(nvarE .gt. 0 .and. ndofE .gt. 0) then
      do ivar=1,nvarE
         NODES(Nod)%dof%zdofE(ivar,1:ndofE) = Buf(j+1:j+ndofE)
         j = j + ndofE
      enddo
   endif
!..extract H(div) dof from buffer
   if(nvarV .gt. 0 .and. ndofV .gt. 0) then
      do ivar=1,nvarV
         NODES(Nod)%dof%zdofV(ivar,1:ndofV) = Buf(j+1:j+ndofV)
         j = j + ndofV
      enddo
   endif
!..extract L2 dof from buffer
   if(nvarQ .gt. 0 .and. ndofQ .gt. 0) then
      do ivar=1,nvarQ
         NODES(Nod)%dof%zdofQ(ivar,1:ndofQ) = Buf(j+1:j+ndofQ)
         j = j + ndofQ
      enddo
   endif
end subroutine unpack_dof_buf
!
end module par_mesh





!
! -----------------------------------------------------------------------
!
!    routine   - get_ghost_subd
!
! -----------------------------------------------------------------------
!
!    last      - Dec 22
!
!    purpose   - Gets all mdle's sharing nodes with subdomain;
!                stores in ELEM_GHOST and ELEM_INTERF
!
!                ELEM_GHOST : ELEM_SUBD plus elements touching my subd
!                ELEM_INTERF: mdle nodes touching subd interface
!
!                both are stored in natural order of elements
!
!  NOTE: This can be done more simply by just looping all elements in
!        mesh, getting unconstrained elements is expensive though (for
!        a mesh of 1 million elements, it takes 12 threads ~100 seconds
!        to get nodes for all elemes). This method is much faster for
!        large distributed meshes.
!
!-----------------------------------------------------------------------
!
   subroutine get_ghost_subd
!
   use data_structure3D
   use bitvisit
   use par_mesh,           only: DISTRIBUTED, HOST_MESH
   use mpi_param
   use MPI,                only: MPI_IN_PLACE, MPI_INTEGER, &
                                 MPI_SUM, MPI_COMM_WORLD
!
   implicit none
!
   integer :: nodesl(27), norientl(27), nodm(MAXNODM)
   integer :: nverts(NRELES_SUBD), verts(16,NRELES_SUBD)
   integer :: my_elems(NRELES_SUBD)
   integer :: nrelem_procs(NUM_PROCS), displs(NUM_PROCS)
!
   integer :: iel, mdle, subd, nrnodm, nrvert, nrelem, loc
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
   call MPI_Allreduce(MPI_IN_PLACE,vert_nprocs,nrvert,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, ierr)
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
   call MPI_Allreduce(MPI_IN_PLACE,elem_nrvert,nrelem_interf,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, ierr)
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
   call MPI_Allreduce(MPI_IN_PLACE,elem_verts,sum,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, ierr)
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
#if DEBUG_MODE
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







!
!-----------------------------------------------------------------------
   subroutine check_interf_elems
!
   use data_structure3D
   use bitvisit
   use par_mesh,           only: DISTRIBUTED, HOST_MESH
   use mpi_param
!
   implicit none
!
   integer :: nodesl(27), norientl(27), nodm(MAXNODM)
!
   integer :: iel, mdle, subd, nrnodm, nrvert, nrelem, loc
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






!
!-----------------------------------------------------------------------
   subroutine check_ghost_elems
!
   use data_structure3D
   use bitvisit
   use par_mesh,           only: DISTRIBUTED, HOST_MESH
   use mpi_param
!
   implicit none
!
   integer :: nodesl(27), norientl(27), nodm(MAXNODM)
!
   integer :: iel, mdle, subd, nrnodm, nrvert, nrelem, loc
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

