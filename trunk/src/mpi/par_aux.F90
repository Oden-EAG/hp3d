!
!..auxiliary subroutines for distributed mesh (par_mesh module)
!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!     subroutine:          collect_dofs
!
!     last modified:       July 2019
!
!     purpose:             collect solution degrees of freedom from all
!                          processors into the ROOT processor
!
!----------------------------------------------------------------------
subroutine collect_dofs()
!
   use data_structure3D
   use par_mesh
   use mpi_param, only: ROOT,RANK
   use MPI      , only: MPI_COMM_WORLD,MPI_STATUS_SIZE,  &
                        MPI_COMPLEX16,MPI_REAL8
!
   implicit none
!
!..MPI variables
   integer :: ierr
   integer :: tag, count, src, dest
   integer :: stat(MPI_STATUS_SIZE)
   VTYPE, allocatable :: buf(:)
!
!..auxiliary variables
   integer :: nodm(MAXNODM)
   integer :: mdle, nod, inod, iel, subd, vis, nrnodm, nrdof_nod
   integer :: iprint = 0
!
!----------------------------------------------------------------------
!
   if ((.not. DISTRIBUTED) .or. HOST_MESH) then
      write(*,*) 'collect_dofs: mesh is not distributed (or on host).'
      goto 190
   endif
!
!..collect degrees of freedom from every element
!..use visit flags to avoid re-sending information
   call reset_visit
!
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
!     if mdle node is already in ROOT's subdomain
      if (subd .eq. ROOT) cycle
!  ...3a. get list of nodes associated with mdle node
      call get_elem_nodes(mdle, nodm,nrnodm)
!  ...3b. iterate over list of nodes
      do inod=1,nrnodm
         nod = nodm(inod)
         call get_visit(nod, vis)
         if (vis .eq. 1) goto 150
         if (RANK .ne. ROOT .and. RANK .ne. subd) goto 140
!     ...3c. Calculate nodal degrees of freedom, and allocate buffer
         call get_dof_buf_size(nod, nrdof_nod)
         if (nrdof_nod .eq. 0) cycle
         allocate(buf(nrdof_nod))
         buf = ZERO
!     ...3d. if current subdomain is my subdomain, send data
         if (RANK .eq. subd) then
!        ...3e. pack buffer with DOF data
            call pack_dof_buf(nod,nrdof_nod, buf)
!        ...3f. send buffer to new subdomain
            if (iprint .eq. 1) then
               write(6,1000) '[', RANK, ']: ', &
                  'Sending data to [',ROOT,'], nod = ',nod
            endif
            count = nrdof_nod; dest = ROOT; tag = nod
            call MPI_SEND(buf,count,MPI_VTYPE,dest,tag,MPI_COMM_WORLD,ierr)
!     ...3d. if new subdomain is my subdomain, receive data
         else if (RANK .eq. ROOT) then
!        ...3e. receive buffer from old subdomain
            if (iprint .eq. 1) then
               write(6,1000) '[', RANK, ']: ',   &
                  'Receiving data from [',subd,'], nod = ',nod
            endif
            count = nrdof_nod; src = subd; tag = nod
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,stat,ierr)
!        ...3f. unpack DOF data from buffer
            call alloc_nod_dof(nod)
            call unpack_dof_buf(nod,nrdof_nod,buf)
         endif
         deallocate(buf)
     140 continue
         call set_visit(nod)
         call set_subd(nod,ROOT)
     150 continue
         if ((RANK .ne. ROOT) .and. (RANK .eq. subd)) then
            call dealloc_nod_dof(nod)
         endif
      enddo
      call set_subd(mdle, ROOT)
   enddo
   1000 format(A,I2,A,A,I4,A,I6)
!
   HOST_MESH = .true.
   write(*,*) 'collect_dofs: mesh is now on host, HOST_MESH = .true.'
!
   190 continue
!
end subroutine collect_dofs
!
!
!----------------------------------------------------------------------
!
!     subroutine:          print_partition
!
!     last modified:       July 2019
!
!     purpose:             print current partition of distributed mesh
!
!----------------------------------------------------------------------
subroutine print_partition()
!
   use data_structure3D
   use par_mesh , only: DISTRIBUTED
   use mpi_param, only: RANK,ROOT,NUM_PROCS
!
   implicit none
!
   integer :: par(NRELES)
   integer :: iel,j,k,l,nreles_subd,mdle,subd
!
!----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'print_partition: mesh is not distributed.'
      goto 290
   endif
!
   if (RANK .ne. ROOT) goto 290
!
   do j=0,NUM_PROCS-1
      nreles_subd = 0
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         call get_subd(mdle, subd)
         if (j .eq. subd) then
            nreles_subd = nreles_subd+1
            par(nreles_subd) = mdle
         endif
      enddo
      write(6,2000) 'partition [', j, '] : '
      k = 0
      do while (k .lt. nreles_subd)
         l = MIN(k+10,nreles_subd)
         write(6,2010) '     ',par(k+1:l)
         k = l
      enddo
   enddo
  2000 format(A,I4,A)
  2010 format(A,<l-k>I6)
!
  290 continue
!
end subroutine print_partition
!
!
!----------------------------------------------------------------------
!
!     subroutine:          print_subd
!
!     last modified:       July 2019
!
!     purpose:             print current subdomain for each process
!
!----------------------------------------------------------------------
subroutine print_subd()
!
   use data_structure3D
   use par_mesh , only: DISTRIBUTED
   use mpi_param, only: RANK,ROOT
!
   implicit none
!
   integer :: sub(NRNODS)
   integer :: k,l,nod,nrnod_subd,subd,vis
!
!----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'print_subd: mesh is not distributed.'
      goto 390
   endif
!
   call reset_visit
   nrnod_subd = 0
   do nod=1,NRNODS
      call get_subd(nod, subd)
      if (RANK .ne. subd) cycle
      call get_visit(nod, vis)
      if (vis .eq. 1) cycle
      nrnod_subd = nrnod_subd + 1
      sub(nrnod_subd) = nod
      call set_visit(nod)
   enddo
   write(6,3000) 'subdomain [', RANK, '] : '
   k = 0
   do while (k .lt. nrnod_subd)
      l = MIN(k+10,nrnod_subd)
      write(6,3010) '     ',sub(k+1:l)
      k = l
   enddo
 3000 format(A,I3,A)
 3010 format(A,<l-k>I5)
!
  390 continue
!
end subroutine print_subd
