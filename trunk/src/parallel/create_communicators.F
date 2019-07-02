      subroutine create_communicators

      use communicators
      use parallelism

#include "syscom.blk"      
#include "mpif.h"

      integer :: group_comm_world
      integer :: processors(MAXPROC)
      integer :: comm_myrank_local

      iprint=0
#if PARALLEL_MODE
      if(iprint==1)then
        write(*,*)MYRANK,'NRPROC',NRPROC
      endif
c  ...get the group underlying MPI_COMM_WORLD
      call mpi_comm_group(MPI_COMM_WORLD,group_comm_world,i1)
      if ((i1)/=0) then
        write(*,*)MYRANK,':',
     .  'main: error calling mpi_comm_group!'
        call abort
      endif
      if(iprint.eq.1)then
        write(*,*)MYRANK,'got group',group_comm_world
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c  ...create the new group
c  ...insert ONLY current process
      do i=1,NRPROC
        processors(1)=i-1
        call mpi_group_incl(group_comm_world,1,processors,
     .    MYRANK_GROUP(i),i2)
        if ((i2)/=0) then
          write(*,*)MYRANK,':',
     .    'main: error calling mpi_group_incl!'
          call abort
        endif
        if(iprint.eq.1)then        
          write(*,*)MYRANK,'processors',processors(1)
        endif
      enddo
      if(iprint.eq.1)then
        write(*,*)MYRANK,'MYRANK_GROUP',MYRANK_GROUP(1:NRPROC+1)
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c  ...create the new communicator
      do i=1,NRPROC
        call mpi_comm_create(MPI_COMM_WORLD,MYRANK_GROUP(i),
     .    comm_myrank_local,i3)
        COMM_MYRANK(i)=comm_myrank_local
        if ((i3)/=0) then
          write(*,*)MYRANK,':',
     .    'main: error calling mpi_com_create!'
          call abort
        endif
      enddo
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(iprint.eq.1)then
        write(*,*)MYRANK,'COMM_MYRANK',COMM_MYRANK(1:NRPROC+1)
      endif
#endif
      end
