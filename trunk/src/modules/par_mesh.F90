!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!> @name       par_mesh
!
!> @date       Sep 2023
!
!> @brief            provides functionality for distributing
!                          degrees of freedom in the FE mesh
!
!----------------------------------------------------------------------
module par_mesh
!
   use data_structure3D
   use environment,    only: QUIET_MODE
   use mpi_wrapper
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
!> @name    distr_mesh
!> @brief   (re-)distribute the current mesh to MPI processes
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
   integer :: nrnodm, nrdof_nod
   real(8) :: start_time,end_time
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
!
   if (iprint .eq. 1) then
      write(6,100) 'start distr_mesh, DISTRIBUTED = ', DISTRIBUTED
   endif
   100 format(A,L2)
#endif
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
   if(.not.QUIET_MODE .and. RANK.eq.ROOT) write(*,110) end_time - start_time
   110 format(' partition  : ',f12.5,'  seconds')
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
!     ...Calculate nodal degrees of freedom
         call get_dof_buf_size(nod, nrdof_nod)
         if ((.not.EXCHANGE_DOF) .or. (nrdof_nod.eq.0)) then
!        ...even if nrdof_nod=0, the node%dof pointer must be associated
            if (subd_next(iel) .eq. RANK) call alloc_nod_dof(nod)
            cycle
         endif
!     ...3c. Allocate exchange buffer
         allocate(buf(nrdof_nod))
         buf = ZERO
!     ...3d. if current subdomain is my subdomain, send data
         if (subd .eq. RANK) then
!        ...3e. pack buffer with DOF data
            call pack_dof_buf(nod,nrdof_nod, buf)
!        ...3f. send buffer to new subdomain
#if HP3D_DEBUG
            if (iprint .eq. 1) then
               write(6,130) '[', RANK, ']: ', &
                  'Sending data to [',subd_next(iel),'], nod = ',nod
            endif
#endif
            count = nrdof_nod; dest = subd_next(iel); tag = mod(nod,200000)
            call MPI_SEND(buf,count,MPI_VTYPE,dest,tag,MPI_COMM_WORLD,ierr)
            if (ierr .ne. MPI_SUCCESS) then
               write(6,*) 'MPI_SEND failed. stop.'
               stop
            endif
!     ...3d. if new subdomain is my subdomain, receive data
         else if (subd_next(iel) .eq. RANK) then
!        ...3e. receive buffer from old subdomain
#if HP3D_DEBUG
            if (iprint .eq. 1) then
               write(6,130) '[', RANK, ']: ',   &
                  'Receiving data from [',subd,'], nod = ',nod
            endif
#endif
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
   if(.not.QUIET_MODE .and. RANK.eq.ROOT) write(*,120) end_time - start_time
   120 format(' migration  : ',f12.5,'  seconds')
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
      if ((subd.ne.RANK) .and. Is_active(nod)) then
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
#if HP3D_DEBUG
   if (iprint .eq. 1) then
      write(6,100) 'end   distr_mesh, DISTRIBUTED = .true.'
   endif
#endif
!
end subroutine distr_mesh
!
!----------------------------------------------------------------------
!> @name    get_elem_nodes
!> @brief   get (unconstrained) nodes associated with an element
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
!> @name    set_subd_elem
!> @brief   set subd values for (unconstrained) nodes associated
!                 with an element
!----------------------------------------------------------------------
subroutine set_subd_elem(Mdle)
   integer, intent(in) :: Mdle
   integer :: nodm(MAXNODM),nodesl(27)
   integer :: i,nrnodm,subd,ntype
   call get_elem_nodes(Mdle, nodesl,nodm,nrnodm)
   call get_subd(Mdle, subd)
   ntype = NODES(Mdle)%ntype
!..marking local node list is necessary so that activate_sons works
!  correctly when breaking the element
   do i=1,nvert(ntype)+nedge(ntype)+nface(ntype)
      call set_subd(nodesl(i),subd)
   enddo
   do i=1,nrnodm
      call set_subd(nodm(i),subd)
   enddo
end subroutine set_subd_elem
!
!----------------------------------------------------------------------
!> @name    alloc_nod_dof
!> @brief   allocate all dofs associated with a node
!----------------------------------------------------------------------
subroutine alloc_nod_dof(Nod)
   integer, intent(in) :: Nod
   integer :: ndofH,ndofE,ndofV,ndofQ
   integer :: nvarH,nvarE,nvarV,nvarQ
   integer :: icase
!..calculate ndof,nvar for this node
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
   icase = NODES(Nod)%case
   nvarH = NREQNH(icase)*NRRHS
   nvarE = NREQNE(icase)*NRRHS
   nvarV = NREQNV(icase)*NRRHS
   nvarQ = NREQNQ(icase)*NRRHS
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
      allocate(NODES(Nod)%dof%zdofH(nvarH, ndofH, NRCOMS))
      NODES(Nod)%dof%zdofH = ZERO
   endif
!..allocate H(curl) DOFs
   if (.not. associated(NODES(Nod)%dof%zdofE) .and. (ndofE .gt. 0)) then
      allocate(NODES(Nod)%dof%zdofE(nvarE, ndofE, NRCOMS))
      NODES(Nod)%dof%zdofE = ZERO
   endif
!..allocate H(div) DOFs
   if (.not. associated(NODES(Nod)%dof%zdofV) .and. (ndofV .gt. 0)) then
      allocate(NODES(Nod)%dof%zdofV(nvarV, ndofV, NRCOMS))
      NODES(Nod)%dof%zdofV = ZERO
   endif
!..allocate L2 DOFs
   if (.not. associated(NODES(Nod)%dof%zdofQ) .and. (ndofQ .gt.  0)) then
      allocate(NODES(Nod)%dof%zdofQ(nvarQ, ndofQ, NRCOMS))
      NODES(Nod)%dof%zdofQ = ZERO
   endif
end subroutine alloc_nod_dof
!
!----------------------------------------------------------------------
!> @name    dealloc_nod_dof
!> @brief   deallocate all dofs associated with a node
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
!> @name    get_dof_buf_size
!> @brief   calculate buffer size for message passing of nodal dofs
!----------------------------------------------------------------------
subroutine get_dof_buf_size(Nod, Nrdof_nod)
   integer, intent(in)  :: Nod
   integer, intent(out) :: Nrdof_nod
   integer :: ndofH,ndofE,ndofV,ndofQ
   integer :: nvarH,nvarE,nvarV,nvarQ
   integer :: icase
   call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)
   icase = NODES(Nod)%case
   nvarH = NREQNH(icase)*NRRHS
   nvarE = NREQNE(icase)*NRRHS
   nvarV = NREQNV(icase)*NRRHS
   nvarQ = NREQNQ(icase)*NRRHS
   Nrdof_nod = ndofH*NDIMEN + ndofH*nvarH + ndofE*nvarE + ndofV*nvarV + ndofQ*nvarQ
end subroutine get_dof_buf_size
!
!----------------------------------------------------------------------
!> @name    pack_dof_buf
!> @brief   pack nodal dofs into contiguous buffer for message passing
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
   nvarH = NREQNH(icase)*NRRHS
   nvarE = NREQNE(icase)*NRRHS
   nvarV = NREQNV(icase)*NRRHS
   nvarQ = NREQNQ(icase)*NRRHS
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
         Buf(j+1:j+ndofH) = NODES(Nod)%dof%zdofH(ivar,1:ndofH,N_COMS)
         j = j + ndofH
      enddo
   endif
!..add H(curl) dof to buffer
   if(nvarE .gt. 0 .and. ndofE .gt. 0) then
      do ivar=1,nvarE
         Buf(j+1:j+ndofE) = NODES(Nod)%dof%zdofE(ivar,1:ndofE,N_COMS)
         j = j + ndofE
      enddo
   endif
!..add H(div) dof to buffer
   if(nvarV .gt. 0 .and. ndofV .gt. 0) then
      do ivar=1,nvarV
         Buf(j+1:j+ndofV) = NODES(Nod)%dof%zdofV(ivar,1:ndofV,N_COMS)
         j = j + ndofV
      enddo
   endif
!..add L2 dof to buffer
   if(nvarQ .gt. 0 .and. ndofQ .gt. 0) then
      do ivar=1,nvarQ
         Buf(j+1:j+ndofQ) = NODES(Nod)%dof%zdofQ(ivar,1:ndofQ,N_COMS)
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
!> @name    unpack_dof_buf
!> @brief   unpack nodal dofs from contiguous buffer
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
   nvarH = NREQNH(icase)*NRRHS
   nvarE = NREQNE(icase)*NRRHS
   nvarV = NREQNV(icase)*NRRHS
   nvarQ = NREQNQ(icase)*NRRHS
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
         NODES(Nod)%dof%zdofH(ivar,1:ndofH,N_COMS) = Buf(j+1:j+ndofH)
         j = j + ndofH
      enddo
   endif
!..extract H(curl) dof from buffer
   if(nvarE .gt. 0 .and. ndofE .gt. 0) then
      do ivar=1,nvarE
         NODES(Nod)%dof%zdofE(ivar,1:ndofE,N_COMS) = Buf(j+1:j+ndofE)
         j = j + ndofE
      enddo
   endif
!..extract H(div) dof from buffer
   if(nvarV .gt. 0 .and. ndofV .gt. 0) then
      do ivar=1,nvarV
         NODES(Nod)%dof%zdofV(ivar,1:ndofV,N_COMS) = Buf(j+1:j+ndofV)
         j = j + ndofV
      enddo
   endif
!..extract L2 dof from buffer
   if(nvarQ .gt. 0 .and. ndofQ .gt. 0) then
      do ivar=1,nvarQ
         NODES(Nod)%dof%zdofQ(ivar,1:ndofQ,N_COMS) = Buf(j+1:j+ndofQ)
         j = j + ndofQ
      enddo
   endif
end subroutine unpack_dof_buf
!
end module par_mesh
