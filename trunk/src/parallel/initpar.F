c----------------------------------------------------------------------
c
c   file name          - initpar
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jun 07
c
c   purpose            - initialization of parameters for parallel run
c
c   usage              - main
c
c   arguments          - none
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine initpar

      use parallelism

#include "mpif.h"
      character(4) :: temp_printrank

#if PARALLEL_MODE
c  ...initialize mpi
      call mpi_init(i1)
      call mpi_comm_size(MPI_COMM_WORLD,NRPROC,i2)
      call mpi_comm_rank(MPI_COMM_WORLD,MYRANK,i3)

      if ((i1+i2+i3)/=0) then
        write(*,*)MYRANK,':',
     .  'main: error initializing MPI!'
        call abort
      endif
      call create_communicators
#else
      NRPROC=1
      MYRANK=0
#endif
      call int2str(MYRANK, temp_printrank)
      if(MYRANK.lt.10)then
        PRINTRANK = "000"//temp_printrank
      elseif(MYRANK.lt.100)then
        PRINTRANK = "00"//temp_printrank
      elseif(MYRANK.lt.1000)then
        PRINTRANK = "0"//temp_printrank
      elseif(MYRANK.lt.10000)then
        PRINTRANK = temp_printrank
      else
        write(*,*)'more then 10000 processors'
        stop
      endif
      
      end
