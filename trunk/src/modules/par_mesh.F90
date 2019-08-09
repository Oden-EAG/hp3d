!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!     module:              par_mesh
!
!     last modified:       July 2019
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
   integer :: ierr
   integer :: tag, count, src, dest
   VTYPE, allocatable :: buf(:)
!
!..auxiliary variables
   integer :: subd_next(NRELES), nodm(MAXNODM)
   integer :: i, iel, inod, iproc, mdle, nod, subd, subd_size
   integer :: nrnodm, nrdof_nod
   integer :: iprint = 0
!
   real(8) :: MPI_Wtime,start_time,end_time
!
   if (iprint .eq. 1) then
      write(6,100) 'start distr_mesh, DISTRIBUTED = ', DISTRIBUTED
   endif
   100 format(A,L2)
!
!..1. Determine new partition
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   if ((ZOLTAN_LB .eq. 0) .or. (.not. DISTRIBUTED)) then
      iproc=0
      do iel=1,NRELES
!     ...decide subdomain for iel
         subd_next(iel) = iproc
         subd_size = (NRELES+NUM_PROCS-1)/NUM_PROCS
         iproc = iel/subd_size
      enddo
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
      call get_elem_nodes(mdle, nodm,nrnodm)
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
            count = nrdof_nod; dest = subd_next(iel); tag = nod
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
            count = nrdof_nod; src = subd; tag = nod
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
!..6. Set subdomain values for all nodes within subdomain
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
!     purpose:    get nodes associated with an element
!----------------------------------------------------------------------
subroutine get_elem_nodes(Mdle, Nodm,Nrnodm)
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Nodm(MAXNODM)
   integer, intent(out) :: Nrnodm
   integer :: nodesl(27),norientl(27)
   Nodm(1:MAXNODM) = 0; Nrnodm = 0
   call get_connect_info(Mdle, nodesl,norientl)
   call logic_nodes(Mdle,nodesl, Nodm,Nrnodm)
end subroutine get_elem_nodes
!
!----------------------------------------------------------------------
!     routine:    set_subd_elem
!     purpose:    set subd values for nodes associated with an element
!----------------------------------------------------------------------
subroutine set_subd_elem(Mdle)
   integer, intent(in) :: Mdle
   integer :: nodm(MAXNODM)
   integer :: i, nrnodm, subd
   call get_elem_nodes(Mdle, nodm,nrnodm)
   call get_subd(Mdle, subd)
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
!..allocate geometry DOFs
   if (.not. associated(NODES(Nod)%coord) .and. (ndofH .gt. 0)) then
      allocate( NODES(Nod)%coord(NDIMEN, ndofH))
   endif
   NODES(Nod)%coord = 0.d0
!..allocate H1 DOFS
   if (.not. associated(NODES(Nod)%zdofH) .and. (ndofH .gt. 0)) then
      allocate( NODES(Nod)%zdofH(nvarH, ndofH))
   endif
   NODES(Nod)%zdofH = ZERO
!..allocate H(curl) DOFs
   if (.not. associated(NODES(Nod)%zdofE) .and. (ndofE .gt. 0)) then
      allocate( NODES(Nod)%zdofE(nvarE, ndofE))
   endif
   NODES(Nod)%zdofE = ZERO
!..allocate H(div) DOFs
   if (.not. associated(NODES(Nod)%zdofV) .and. (ndofV .gt. 0)) then
      allocate( NODES(Nod)%zdofV(nvarV, ndofV))
   endif
   NODES(Nod)%zdofV = ZERO
!..allocate L2 DOFs
   if (.not. associated(NODES(Nod)%zdofQ) .and. (ndofQ .gt.  0)) then
      allocate( NODES(Nod)%zdofQ(nvarQ, ndofQ))
   endif
   NODES(Nod)%zdofQ = ZERO
end subroutine alloc_nod_dof
!
!----------------------------------------------------------------------
!     routine:    dealloc_nod_dof
!     purpose:    deallocate all dofs associated with a node
!----------------------------------------------------------------------
subroutine dealloc_nod_dof(Nod)
   integer, intent(in) :: Nod
   if (associated(NODES(Nod)%coord)) deallocate(NODES(Nod)%coord) 
   if (associated(NODES(Nod)%zdofH)) deallocate(NODES(Nod)%zdofH)
   if (associated(NODES(Nod)%zdofE)) deallocate(NODES(Nod)%zdofE)
   if (associated(NODES(Nod)%zdofV)) deallocate(NODES(Nod)%zdofV)
   if (associated(NODES(Nod)%zdofQ)) deallocate(NODES(Nod)%zdofQ)
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
         Buf(j+1:j+ndofH) = NODES(Nod)%coord(ivar,1:ndofH)
         j = j + ndofH
      enddo
   endif
!..add H1 dof to buffer
   if(nvarH .gt. 0 .and. ndofH .gt. 0) then
      do ivar=1,nvarH
         Buf(j+1:j+ndofH) = NODES(Nod)%zdofH(ivar,1:ndofH)
         j = j + ndofH
      enddo
   endif
!..add H(curl) dof to buffer
   if(nvarE .gt. 0 .and. ndofE .gt. 0) then
      do ivar=1,nvarE
         Buf(j+1:j+ndofE) = NODES(Nod)%zdofE(ivar,1:ndofE)
         j = j + ndofE
      enddo
   endif
!..add H(div) dof to buffer
   if(nvarV .gt. 0 .and. ndofV .gt. 0) then
      do ivar=1,nvarV
         Buf(j+1:j+ndofV) = NODES(Nod)%zdofV(ivar,1:ndofV)
         j = j + ndofV
      enddo
   endif
!..add L2 dof to buffer
   if(nvarQ .gt. 0 .and. ndofQ .gt. 0) then
      do ivar=1,nvarQ
         Buf(j+1:j+ndofQ) = NODES(Nod)%zdofQ(ivar,1:ndofQ)
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
         NODES(Nod)%coord(ivar,1:ndofH) = Buf(j+1:j+ndofH)
         j = j + ndofH
      enddo
   endif
!..extract H1 dof from buffer
   if(nvarH .gt. 0 .and. ndofH .gt. 0) then
      do ivar=1,nvarH
         NODES(Nod)%zdofH(ivar,1:ndofH) = Buf(j+1:j+ndofH)
         j = j + ndofH
      enddo
   endif
!..extract H(curl) dof from buffer
   if(nvarE .gt. 0 .and. ndofE .gt. 0) then
      do ivar=1,nvarE
         NODES(Nod)%zdofE(ivar,1:ndofE) = Buf(j+1:j+ndofE)
         j = j + ndofE
      enddo
   endif
!..extract H(div) dof from buffer
   if(nvarV .gt. 0 .and. ndofV .gt. 0) then
      do ivar=1,nvarV
         NODES(Nod)%zdofV(ivar,1:ndofV) = Buf(j+1:j+ndofV)
         j = j + ndofV
      enddo
   endif
!..extract L2 dof from buffer
   if(nvarQ .gt. 0 .and. ndofQ .gt. 0) then
      do ivar=1,nvarQ
         NODES(Nod)%zdofQ(ivar,1:ndofQ) = Buf(j+1:j+ndofQ)
         j = j + ndofQ
      enddo
   endif
end subroutine unpack_dof_buf
!
end module par_mesh
