!
!..verification subroutines
!
#include "implicit_none.h"
!----------------------------------------------------------------------
!
!     subroutine:          par_verify
!
!     last modified:       Oct 2019
!
!     purpose:             verify the par_mesh module functionality
!
!----------------------------------------------------------------------
subroutine par_verify()
!
   use mpi_param, only: RANK,ROOT
   use par_mesh , only: DISTRIBUTED
   use MPI      , only: MPI_COMM_WORLD,MPI_INTEGER,MPI_MIN
!
   implicit none
!
   integer :: ipass,ierr,buf,count
!
   if (.not. DISTRIBUTED) then
      if (RANK.eq.ROOT) write(*,*) 'par_verify: mesh is not distributed.'
      goto 190
   endif
!
   !if (RANK.eq.ROOT) write(*,*) 'par_verify: Checking mesh consistency...'
   call mesh_consistency(ipass)
   count = 1
   call MPI_ALLREDUCE(ipass,buf,count,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
   if (buf .eq. 1) then
      if (RANK.eq.ROOT)  write(*,*) 'par_verify: MESH CONSISTENCY TEST PASSED.'
   else
      if (RANK.eq.ROOT)  write(*,*) 'par_verify: MESH CONSISTENCY TEST FAILED.'
      stop
   endif
!
!   write(*,*) 'par_verify: Checking collection and distribution of dofs...'
!   call verify_dof(ipass)
!   if (ipass .eq. 1) then
!      write(*,*) 'par_verify: COLLECT/DISTRIBUTE TEST PASSED.'
!   else
!      write(*,*) 'par_verify: COLLECT/DISTRIBUTE TEST FAILED.'
!      call pause
!   endif
!
  190 continue
!
end subroutine par_verify
!
!
!----------------------------------------------------------------------
!
!     subroutine:          mesh_consistency
!
!     last modified:       July 2019
!
!     purpose:             check mesh consistency for distributed mesh
!
!     arguments:
!         out              ipass = 0 (test passed)
!                                = 1 (test failed)
!
!----------------------------------------------------------------------
subroutine mesh_consistency(ipass)
!
   use data_structure3D
   use par_mesh , only: DISTRIBUTED
   use mpi_param, only: ROOT,RANK
   use MPI      , only: MPI_INTEGER,MPI_COMM_WORLD
!
   implicit none
!
!..test success parameter
   integer, intent(out) :: ipass
!
!..auxiliary variables
   integer :: i,iel,mdle,nod,nrnodm,count,src,ierr
   integer, allocatable :: val(:),buf(:)
!
!..element nodes
   integer, dimension(27)      :: nodesl, norientl
   integer, dimension(MAXNODM) :: nodm
!
!----------------------------------------------------------------------
!
   ipass = 1
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'mesh_consistency: mesh is not distributed.'
      goto 290
   endif
!
!  --------------------------------
!  1. check general mesh parameters
!  --------------------------------
   allocate(val(6),buf(6))
   buf(1:6) = 0
   val(1) = NRELES ; val(2) = NRNODS
   val(3) = NRDOFSH; val(4) = NRDOFSE
   val(5) = NRDOFSV; val(6) = NRDOFSQ
!
   if (rank .eq. ROOT) then
      buf(1:6) = val(1:6)
   endif
!
   count = 6; src = ROOT
   call MPI_BCAST (buf,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
!
   if (rank .ne. ROOT) then
      do i=1,count
         if (buf(i) .ne. val(i)) then
            write(6,1210) '[', rank, ']: ','mesh inconsistency (1): ', &
                          'i,val,buf = ', i, val(i), buf(i)
            ipass = 0
            deallocate(val,buf); return
         endif
      enddo
   endif
 1210 format(A,I3,A,A,A,I7,I10,I10)
!
   deallocate(val,buf)
!
!  --------------------------------------------
!  2. check mdle node subdomains (partitioning)
!  --------------------------------------------
   count = 0
   do nod=1,NRNODS
      if (Is_middle(nod)) then
         count = count+1
      endif
   enddo
   allocate(val(count),buf(count))
   buf(1:count) = 0
   i = 0
   do nod=1,NRNODS
      if (Is_middle(nod)) then
         i = i+1
         val(i) = NODES(nod)%subd
      endif
   enddo
   if (rank .eq. ROOT) then
      buf(1:count) = val(1:count)
   endif
!
   src = ROOT
   call MPI_BCAST (buf,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
!
   if (rank .ne. ROOT) then
      do i=1,count
         if (buf(i) .ne. val(i)) then
            write(6,1210) '[', rank, ']: ','mesh inconsistency (2): ', &
                          'i,val,buf = ', i, val(i), buf(i)
            ipass = 0
            deallocate(val,buf); return
         endif
      enddo
   endif
!
   deallocate(val,buf)
!
!  ----------------------------------------
!  3. check activation flags in NODES array
!  ----------------------------------------
   !goto 290
   allocate(val(NRNODS),buf(NRNODS))
   buf(1:NRNODS) = 0
!
   do nod=1,NRNODS
      if(Is_active(nod)) then
         val(nod) = 1
      else
         val(nod) = 0
      endif
   enddo
   if (rank .eq. ROOT) then
      buf(1:NRNODS) = val(1:NRNODS)
   endif
!
   count = NRNODS; src = ROOT
   call MPI_BCAST (buf,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
!
   if (rank .ne. ROOT) then
      do i=1,count
         if (buf(i) .ne. val(i)) then
            write(6,1210) '[', rank, ']: ','mesh inconsistency (3): ', &
                          'i,val,buf = ', i, val(i), buf(i)
            ipass = 0
            deallocate(val,buf); return
         endif
      enddo
   endif
!
   deallocate(val,buf)
!
!  ------------------------------------------------
!  4. check modified element nodes within subdomain
!  ------------------------------------------------
!
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      call get_connect_info(mdle, nodesl,norientl)
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
      do i=1,nrnodm
         nod = nodm(i)
         if (Is_inactive(nod)) then
            write(6,1240) '[', rank, ']: ','mesh inconsistency (4): ', &
                           'inactive nod = ', nod
            ipass = 0
            call result
         endif
         if (NODES(nod)%subd .ne. RANK) then
            write(6,1241) '[', rank, ']: ','mesh inconsistency (4): ', &
                          'nod,subd = ', nod,NODES(nod)%subd
            ipass = 0
            call result
         endif
         if (.not.associated(NODES(nod)%dof)) then
            write(6,1240) '[', rank, ']: ','mesh inconsistency (4): ', &
                           'dof not associated, nod = ', nod
            ipass = 0
            call result
         endif
      enddo
   enddo
 1240 format(A,I3,A,A,A,I7)
 1241 format(A,I3,A,A,A,I7,I4)
!
  290 continue
!
end subroutine mesh_consistency
!
!----------------------------------------------------------------------
!
!     subroutine:          verify_dof
!
!     last modified:       July 2019
!
!     purpose:             verify collection and distribution of dofs
!
!     arguments:
!         out              ipass = 0 (test passed)
!                                = 1 (test failed)
!
!----------------------------------------------------------------------
subroutine verify_dof(ipass)
!
   use par_mesh , only: DISTRIBUTED,distr_mesh
   use mpi_param, only: ROOT,RANK
   use MPI      , only: MPI_INTEGER,MPI_COMM_WORLD
!
   implicit none
!
   integer, intent(out) :: ipass
!
!----------------------------------------------------------------------
!
   ipass = 1
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'verify_dof: mesh is not distributed.'
      goto 390
   endif
!
!..1. Run mesh distribution to obtain current partition
   call distr_mesh
!
!..2. Every processor: store subdomain dofs locally

!
!..3. Collect all dofs on ROOT processor, delete from other processors
   call collect_dofs
!
!..4. Distribute dofs from ROOT processor by mesh distribution
   call distr_mesh
!
!..5. Compare locally stored dofs to received dofs (fail test if different)

!
  390 continue
!
end subroutine verify_dof
