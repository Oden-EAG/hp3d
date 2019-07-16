!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!     module:              MPI_param
!     last modified:       July 2019
!
!----------------------------------------------------------------------
module MPI_param
!
   use MPI, only: MPI_THREAD_FUNNELED,MPI_INIT_THREAD,         &
                  MPI_COMM_RANK,MPI_COMM_WORLD,MPI_COMM_SIZE,  &
                  MPI_BARRIER,MPI_FINALIZE
!
   implicit none
!
   integer, parameter :: ROOT = 0
   integer, save      :: RANK = -1
   integer, save      :: NUM_PROCS = -1
   logical, save      :: MPI_PARAM_IS_INIT = .false.
!
   contains
!
!----------------------------------------------------------------------
!     routine:    MPI_param_init
!     purpose:    initialize MPI environment, and set parameters
!----------------------------------------------------------------------
   subroutine MPI_param_init()
!
      integer :: ierr,req,ret
!
!  ...MPI initialization
      req = MPI_THREAD_FUNNELED
      call MPI_INIT_THREAD ( req, ret, ierr )
      if (ret < req) then
         write(*,*) 'MPI_param: MPI_INIT_THREAD, ', &
            'provided thread level < requested thread level. stop.'
         stop
      endif
!
!  ...find out my RANK (process id), and number of processors
      call MPI_COMM_RANK (MPI_COMM_WORLD, RANK, ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, NUM_PROCS, ierr)
!
!  ...set initialization flag
      MPI_PARAM_IS_INIT = .true.
!
   end subroutine MPI_param_init
!
!
!----------------------------------------------------------------------
!     routine:    MPI_param_finalize
!     purpose:    close MPI environment
!----------------------------------------------------------------------
   subroutine MPI_param_finalize()
      integer :: ierr
      call MPI_FINALIZE ( ierr )
      MPI_PARAM_IS_INIT = .false.
   end subroutine MPI_param_finalize
!
!
end module MPI_param
