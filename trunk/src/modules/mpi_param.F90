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
   integer, save      :: INIT = 0
!
   contains
!
!----------------------------------------------------------------------
!     routine:    MPI_param_init
!     purpose:    initialize MPI environment, and set parameters
!----------------------------------------------------------------------
   subroutine MPI_param_init()
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
      INIT = 1
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
      INIT = 0
   end subroutine MPI_param_finalize
!
!
!----------------------------------------------------------------------
!     function:   Is_init
!     purpose:    return initialization status
!----------------------------------------------------------------------
   function Is_init()
      logical Is_init
      if (INIT .eq. 1) then
         Is_init = .true.
      else
         Is_init = .false.
      endif
   end function Is_init
!
end module MPI_param
