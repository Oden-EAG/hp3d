!----------------------------------------------------------------------
!     module:              MPI_param
!     last modified:       July 2019
!----------------------------------------------------------------------
module mpi_param
!
   integer, parameter :: ROOT = 0
   integer, save      :: RANK = -1
   integer, save      :: NUM_PROCS = -1
!
end module mpi_param

