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
   use MPI_param, only: ROOT,RANK
   use MPI      , only: MPI_COMM_WORLD,MPI_SUCCESS,MPI_STATUS_SIZE,  &
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
   integer :: iprint = 1
!
!----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'collect_dofs: mesh is not distributed.'
      goto 90
   endif
!
!..collect degrees of freedom from every element
!..use visit flags to avoid re-sending information
   call reset_visit
!
   mdle = 0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      call get_subd(mdle, subd)
!     if mdle node is already in ROOT's subdomain
      if (subd .eq. ROOT) cycle
!  ...3a. get list of nodes associated with mdle node
      call get_elem_nodes(mdle, nodm,nrnodm)
!  ...3b. iterate over list of nodes
      do inod=1,nrnodm
         nod = nodm(inod)
         call get_visit(nod, vis)
         if (vis .eq. 1) goto 50
         if (RANK .ne. ROOT .and. RANK .ne. subd) goto 40
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
               write(6,100) '[', RANK, ']: ', &
                  'Sending data to [',ROOT,'], nod = ',nod
            endif
            count = nrdof_nod; dest = ROOT; tag = nod
            call MPI_SEND(buf,count,MPI_VTYPE,dest,tag,MPI_COMM_WORLD,ierr)
            if (ierr .ne. MPI_SUCCESS) then
               write(6,*) 'MPI_SEND failed. stop.'
               stop
            endif
!     ...3d. if new subdomain is my subdomain, receive data
         else if (RANK .eq. ROOT) then
!        ...3e. receive buffer from old subdomain
            if (iprint .eq. 1) then
               write(6,100) '[', RANK, ']: ',   &
                  'Receiving data from [',subd,'], nod = ',nod
            endif
            count = nrdof_nod; src = subd; tag = nod
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,stat,ierr)
            if (ierr .ne. MPI_SUCCESS) then
               write(6,*) 'MPI_RECV failed. stop.'
               stop
            endif
!        ...3f. unpack DOF data from buffer
            call alloc_nod_dof(nod)
            call unpack_dof_buf(nod,nrdof_nod,buf)
         endif
         deallocate(buf)
      40 continue
         call set_visit(nod)
      50 continue
         if (RANK .eq. subd) then
            call dealloc_nod_dof(nod)
         endif
      enddo
      call set_subd(mdle, ROOT)
   enddo
   100 format(A,I2,A,A,I4,A,I6)

!..Delete degrees of freedom for NODES if not ROOT processor
   if (RANK .ne. ROOT) then
      do nod=1,NRNODS
         call get_subd(nod, subd)
         if (subd .ne. RANK .and. Is_active(nod)) then

         endif
      enddo
   endif
!
  90 continue
!
end subroutine collect_dofs
