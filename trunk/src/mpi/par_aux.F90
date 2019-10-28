!
!..auxiliary subroutines for distributed mesh (par_mesh module)
!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!     subroutine:          partition_fiber
!
!     last modified:       Aug 2019
!
!     purpose:             partition a waveguide structure into slabs
!                          along z-axis (with equal size per proc)
!
!----------------------------------------------------------------------
subroutine partition_fiber(subd_next)
!
   use data_structure3D
   use mpi_param, only: RANK,ROOT,NUM_PROCS
   use MPI,       only: MPI_COMM_WORLD,MPI_INTEGER,MPI_REAL8,  &
                        MPI_MIN,MPI_MAX,MPI_IN_PLACE
!
   implicit none
!
   integer, intent(out) :: subd_next(NRELES)
!
   integer :: iel,mdle,subd,i,k,nrv
   real(8) :: xnod(NDIMEN,8)
   real(8) :: x3,x3_lo,x3_hi,x3_subd
!
   integer :: count,ierr
!
!----------------------------------------------------------------------
!
   subd_next(1:NRELES) = 0
!
   x3_lo = 100.d0; x3_hi = -100.d0
!
!$OMP PARALLEL DO PRIVATE(mdle,i,nrv,xnod,x3)   &
!$OMP REDUCTION(MIN:x3_lo) REDUCTION(MAX:x3_hi)
   do iel = 1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      call nodcor_vert(mdle, xnod)
      nrv = nvert(NODES(mdle)%type)
      do i = 1,nrv
         x3 = xnod(3,i)
         if (x3 .lt. x3_lo) x3_lo = x3
         if (x3 .gt. x3_hi) x3_hi = x3
      enddo
   enddo
!$OMP END PARALLEL DO
!
   count = 1
   call MPI_ALLREDUCE(MPI_IN_PLACE,x3_lo,count,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,x3_hi,count,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
!
!   write(*,100) '1. x3_lo = ', x3_lo
!   write(*,100) '1. x3_hi = ', x3_hi
   x3_subd = (x3_hi-x3_lo)/NUM_PROCS
!   write(*,100) ' x3_subd = ', x3_subd
!   x3_lo = x3_subd * RANK
!   x3_hi = x3_subd * (RANK+1)
!   write(*,100) '2. x3_lo = ', x3_lo
!   write(*,100) '2. x3_hi = ', x3_hi
 100 format(A,F8.2)
!
!$OMP PARALLEL DO PRIVATE(mdle,subd,i,nrv,xnod,x3) SCHEDULE(DYNAMIC)
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      if (subd .ne. RANK) cycle
      call nodcor_vert(mdle, xnod)
      nrv = nvert(NODES(mdle)%type)
      x3 = 0.d0
      do i = 1,nrv
         x3 = x3 + xnod(3,i)
      enddo
      x3 = x3 / nrv
      subd_next(iel) = int(x3/x3_subd)
   enddo
!$OMP END PARALLEL DO
!
   count = NRELES
   call MPI_ALLREDUCE(MPI_IN_PLACE,subd_next,count,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
!
end subroutine partition_fiber
!
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
      call set_subd(mdle,ROOT)
   enddo
   1000 format(A,I2,A,A,I4,A,I6)
!
   call update_ELEM_ORDER
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
   integer :: iel,j,k,l,nrelem,mdle,subd
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
      nrelem = 0
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         call get_subd(mdle, subd)
         if (j .eq. subd) then
            nrelem = nrelem+1
            par(nrelem) = mdle
         endif
      enddo
      write(6,2000) 'partition [', j, '] : '
      k = 0
      do while (k .lt. nrelem)
         l = MIN(k+10,nrelem)
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
!
!
!----------------------------------------------------------------------
!
!     subroutine:          print_coord
!
!     last modified:       July 2019
!
!     purpose:             print current subdomain for each process
!
!----------------------------------------------------------------------
subroutine print_coord()
!
   use data_structure3D
   use par_mesh , only: DISTRIBUTED
   use mpi_param, only: RANK,ROOT
!
   implicit none
!
   integer :: par(NRELES)
   integer :: iel,i,mdle,nrv
   real*8  :: x(NDIMEN), xnod(NDIMEN,8)
!
!----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'print_partition: mesh is not distributed.'
      return
   endif
!
   write(6,4000) 'partition [', RANK, '] : '
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      call nodcor_vert(mdle, xnod)
      nrv = nvert(NODES(mdle)%type)
      x(1:3) = 0.d0
      do i = 1,nrv
         x(1:3) = x(1:3) + xnod(1:3,i)
      enddo
      x(1:3) = x(1:3) / nrv
      write(*,4010) 'Mdle = ', mdle,', Coords = ', x(1:3)
   enddo
  4000 format(A,I4,A)
  4010 format(A,I5,A,3F6.2)
!
end subroutine print_coord
