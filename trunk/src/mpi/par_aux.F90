!
!..auxiliary subroutines for distributed mesh (par_mesh module)
!
#include "typedefs.h"
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
   integer :: iel,mdle,subd,i,nrv
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
      nrv = nvert(NODES(mdle)%ntype)
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
      nrv = nvert(NODES(mdle)%ntype)
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
subroutine repartition_fiber(Subd_next)
!
   use data_structure3D
   use mpi_param, only: RANK,ROOT,NUM_PROCS
   use MPI,       only: MPI_COMM_WORLD,MPI_INTEGER,MPI_REAL8,  &
                        MPI_SUM,MPI_IN_PLACE
   use stc,       only: stc_get_nrdof
!
   implicit none
!
   integer, intent(out) :: Subd_next(NRELES)
!
   integer :: iel,mdle,subd,i,k,nrv,proc
   real(8) :: xnod(NDIMEN,8)
!
   integer :: count,ierr
!
!..dynamic load balancing
   real(8) :: x3_mdle(NRELES)
   real(8) :: x3_subd(NUM_PROCS)
   integer :: iel_load(NRELES)
   integer :: iel_array(NRELES)
   integer :: nrdofi(NR_PHYSA),nrdofb(NR_PHYSA)
   integer :: load
!
!----------------------------------------------------------------------
!
   Subd_next(1:NRELES) = 0
!
   x3_mdle (1:NRELES) = 0.d0
   iel_load(1:NRELES) = 0
!
!$OMP PARALLEL DO PRIVATE(mdle,subd,i,nrv,xnod,nrdofi,nrdofb)
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      if (subd .ne. RANK) cycle
      call nodcor_vert(mdle, xnod)
      nrv = nvert(NODES(mdle)%ntype)
      x3_mdle(iel) = maxval(xnod(3,1:nrv))
      call stc_get_nrdof(mdle, nrdofi,nrdofb)
      iel_load(iel) = sum(nrdofi) + sum(nrdofb)
   enddo
!$OMP END PARALLEL DO
!
   count = NRELES
   call MPI_ALLREDUCE(MPI_IN_PLACE,x3_mdle,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,iel_load,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!..sort elements by z-coordinate
   do iel=1,NRELES
      iel_array(iel) = iel
   enddo
   call qsort_duplet_z(iel_array,x3_mdle,NRELES,1,NRELES)
!
!..add elements to processor based on z-coordinate and load
   load = (sum(iel_load)+NUM_PROCS-1) / NUM_PROCS
   x3_subd(1:NUM_PROCS) = 0.d0
   proc = 0; k = 0
   do i = 1,NRELES
      iel = iel_array(i)
      k = k + iel_load(iel)
      x3_subd(proc+1) = x3_mdle(i)
      Subd_next(iel) = proc
      if ((k .ge. (proc+1)*load) .and. (i < NRELES)) then
         if (x3_mdle(i+1) > (x3_mdle(i)+1.0d-6)) proc = proc+1
      endif
   enddo
!
   if (RANK.eq.ROOT) then
      do i=1,NUM_PROCS
         write(*,130) 'rank, x3_subd = ', i-1,x3_subd(i)
      enddo
  130 format(A,I4,F10.4)
      !do i=1,NRELES
      !   iel = iel_array(i)
      !   write(*,140) 'i,iel,x3_mdle(i),iel_load(iel) = ',i,iel,x3_mdle(i),iel_load(iel)
      !enddo
  140 format(A,I6,',',I6,',',F10.4,',',I6)
   endif
!
end subroutine repartition_fiber
!
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Oct 2019
!
!    purpose:           - sorts an array of duplets (iel,residual) with
!                         residual (sort key) in ascending order
!                         (initial call needs: First = 1, Last = N)
!
!    arguments:
!           in/out
!                       - Iel_array: 1D integer array (element indices)
!                       - Residuals: 1D real    array (residual values)
!           in
!                       - N        : size of array
!                       - First    : first index of current partition
!                       - Last     : last  index of current partition
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
recursive subroutine qsort_duplet_z(Iel_array,Residuals,N,First,Last)
!
   implicit none
!..declare variables
   integer , intent(in)    :: N,First,Last
   integer , intent(inout) :: Iel_array(N)
   real(8) , intent(inout) :: Residuals(N)
!
   real(8) :: x,pivot
   integer :: i,j,l
!
   pivot = Residuals((First+Last) / 2)
   i = First
   j = Last
!..iterate through the array to be sorted
   do
!  ...find first element from the left that needs to be swapped
      do while ((Residuals(i) < pivot))
         i = i + 1
      end do
!  ...find first element from the right that needs to be swapped
      do while ((pivot < Residuals(j)))
         j = j - 1
      end do
!  ...end loop if no elements need to be swapped
      if (i >= j) exit
!  ...swap the elements
      l = Iel_array(i); Iel_array(i) = Iel_array(j); Iel_array(j) = l
      x = Residuals(i); Residuals(i) = Residuals(j); Residuals(j) = x
      i = i + 1
      j = j - 1
   end do
   if (First < i-1) call qsort_duplet_z(Iel_array,Residuals,N,First,i-1 )
   if (j+1 < Last)  call qsort_duplet_z(Iel_array,Residuals,N,j+1,  Last)
!
end subroutine qsort_duplet_z
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
   integer :: nodm(MAXNODM), nodesl(27)
   integer :: mdle, nod, inod, iel, subd, vis, nrnodm, nrdof_nod
!
   integer :: iprint
   iprint=0
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
      call get_elem_nodes(mdle, nodesl,nodm,nrnodm)
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
   character(16) :: fmt
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
         write(fmt,'("(", I0, "I6)")') l-k
         write(6,fmt) par(k+1:l)
         k = l
      enddo
   enddo
  2000 format(A,I4,A)
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
!     last modified:       Mar 2023
!
!     purpose:             print current subdomain for each process
!
!----------------------------------------------------------------------
subroutine print_subd()
!
   use data_structure3D
   use par_mesh , only: DISTRIBUTED
   use mpi_param, only: RANK,ROOT,NUM_PROCS
   use mpi      , only: MPI_COMM_WORLD
!
   implicit none
!
   integer :: sub(NRNODS)
   integer :: ierr,j,k,l,nod,nrnod_subd,subd
   character(16) :: fmt
!
!----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'print_subd: mesh is not distributed.'
      goto 390
   endif
!
   nrnod_subd = 0
   do nod=1,NRNODS
      call get_subd(nod, subd)
      if (RANK .ne. subd) cycle
      nrnod_subd = nrnod_subd + 1
      sub(nrnod_subd) = nod
   enddo
   do j=0,NUM_PROCS-1
      if (j .ne. RANK) then
         call MPI_BARRIER (MPI_COMM_WORLD, ierr);
         cycle
      endif
      write(6,3000) 'subdomain [', RANK, '] : '
 3000 format(A,I3,A)
      k = 0
      do while (k .lt. nrnod_subd)
         l = MIN(k+10,nrnod_subd)
         write(fmt,'("(",I0,"I6)")') l-k
         write(6,fmt) sub(k+1:l)
         k = l
      enddo
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
   enddo
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
!     last modified:       Mar 2023
!
!     purpose:             print current subdomain for each process
!
!----------------------------------------------------------------------
subroutine print_coord()
!
   use data_structure3D
   use par_mesh , only: DISTRIBUTED
   use mpi_param, only: RANK,ROOT,NUM_PROCS
   use mpi      , only: MPI_COMM_WORLD
!
   implicit none
!
   integer :: ierr,iel,i,j,mdle,nrv
   real(8) :: x(NDIMEN), xnod(NDIMEN,8)
!
!----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'print_partition: mesh is not distributed.'
      return
   endif
!
   do j=0,NUM_PROCS-1
      if (j .ne. RANK) then
         call MPI_BARRIER (MPI_COMM_WORLD, ierr);
         cycle
      endif
      write(6,4000) 'partition [', RANK, '] : '
      do iel=1,NRELES_SUBD
         mdle = ELEM_SUBD(iel)
         call nodcor_vert(mdle, xnod)
         nrv = nvert(NODES(mdle)%ntype)
         x(1:3) = 0.d0
         do i = 1,nrv
            x(1:3) = x(1:3) + xnod(1:3,i)
         enddo
         x(1:3) = x(1:3) / nrv
         write(*,4010) 'Mdle = ', mdle,', Coords = ', x(1:3)
      enddo
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
   enddo
!
  4000 format(A,I4,A)
  4010 format(A,I5,A,3F6.2)
!
end subroutine print_coord
