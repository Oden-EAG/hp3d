!
#include "typedefs.h"
! -----------------------------------------------------------------------
!
!    module name        - par_mumps
!
! -----------------------------------------------------------------------
!
!    latest revision    - Mar 2020
!
!    purpose            - module sets up required workspace for
!                         interfacing with distributed MUMPS solver
!
! -----------------------------------------------------------------------
module par_mumps
!
   use MPI      , only: MPI_COMM_WORLD,MPI_COMM_SELF
   use mpi_param, only: RANK,NUM_PROCS
!
   implicit none
!
#if C_MODE
   include 'zmumps_struc.h'
   type (ZMUMPS_STRUC) mumps_par
   type (ZMUMPS_STRUC) mumps_bub
#else
   include 'dmumps_struc.h'
   type (DMUMPS_STRUC) mumps_par
   type (DMUMPS_STRUC) mumps_bub
#endif
!
   contains
!
! -----------------------------------------------------------------------
!
subroutine mumps_start_par
   call mumps_start(mumps_par,MPI_COMM_WORLD)
end subroutine mumps_start_par
!
subroutine mumps_destroy_par
   call mumps_destroy(mumps_par)
end subroutine mumps_destroy_par
!
subroutine mumps_start(mumps,mpi_comm)
!
#if C_MODE
   type (ZMUMPS_STRUC), intent(inout) :: mumps
#else
   type (DMUMPS_STRUC), intent(inout) :: mumps
#endif
   integer, intent(in) :: mpi_comm
!
!..NOTE:
!..If you plan to run MUMPS sequentially from a parallel MPI application,
!  you need to install the parallel version of MUMPS and pass a
!  communicator containing a single processor to the MUMPS library
!
!..Define a communicator for the package.
   mumps%COMM = mpi_comm
!
!..PAR
!     0 : host is not involved in factorization/solve phases
!     1 : host involved in factorization/solve phases
   mumps%PAR = 1
!
!..SYM (for LU factorization, no Hermitian option in version 5.2.1)
!     0: unsymmetric
!     1: symmetric positive definite
!     2: symmetric
   mumps%SYM = 0
!
!..initialize an instance of the package
   mumps%JOB = -1
!
#if C_MODE
   call zmumps(mumps)
#else
   call dmumps(mumps)
#endif
!
!..diagnostics
   !write(*,100) 'par_mumps: MYID           = ', mumps_par%MYID
   !write(*,110) 'par_mumps: VERSION_NUMBER = ', mumps_par%VERSION_NUMBER
 100 format(A,I3)
 110 format(A,A)
!
!..output stream for error messages
!     0: suppress messages
!     6: stdout
   mumps%icntl(1) = 6
!
!..output stream for diagnostic/statistics/warning messages
!     0: suppress messages
!     6: stdout
   mumps%icntl(2) = 0
!
!..output stream for global information
!     0: suppress messages
!     6: stdout
   mumps%icntl(3) = 0
!
!..print level for error/warning/diagnostic messages
!     0: no output messages
!     1: only error messages
!     2: errors, warnings, and main stats
!     3: errors, warnings, and terse diagnostics
!     4: errors, warnings, and info on input/output params
   mumps%icntl(4) = 1
!
!..icntl(5): matrix input format
!     0: assembled input format
!     1: elemental input format
   mumps%icntl(5) = 0
!
!..icntl(7): choice of sequential pivot ordering tool (not used if icntl(28)=2)
!     0: Approximate Minimum Degree (AMD)
!     1: pivot order set by the user
!     2: Approximate Minimum Fill (AMF)
!     3: Scotch
!     4: PORD
!     5: Metis
!     6: AMD w/ quasi-dense row detection
!     7: automatic value
   mumps%icntl(7) = 5
!
!..icntl(14): percentage increase in estimated workspace
!     [default: 20] - 20% increase in workspace
   mumps%icntl(14) = 80
!
!..icntl(18): distribution strategy of the input matrix
!     0: centralized on host
!     1: user provides matrix structure at analysis, MUMPS returns mapping for entries
!     2: user provides matrix structure at analysis, entries at factorization
!     3: user provides distributed matrix
   mumps%icntl(18) = 3
!
!..icntl(21): determines distribution of the solution vectors
!     0: centralized on host
!     1: distributed
   mumps%icntl(21) = 0
!
!..icntl(28): sequential or parallel analysis
!     0: automatic choice
!     1: sequential computation
!     2: parallel computation
   if (NUM_PROCS < 16) then
      mumps_par%icntl(28) = 1
   else
      mumps_par%icntl(28) = 2
   endif
!
!..icntl(29): choice of parallel pivot ordering tool (not used if icntl(28)=1)
!     0: automatic choice
!     1: PT-Scotch
!     2: ParMetis
   mumps%icntl(29) = 1
   ! parmetis appears to cause 'integer divide by zero error'
   ! when using partitions for small problem
   ! (e.g., try 32 mpi ranks on 64 cubes)
   ! but parmetis is faster than pt-scotch, at least on uniformly refined cube
!
end subroutine mumps_start
!
! -----------------------------------------------------------------------
!
subroutine mumps_destroy(mumps)
!
#if C_MODE
   type (ZMUMPS_STRUC), intent(inout) :: mumps
#else
   type (DMUMPS_STRUC), intent(inout) :: mumps
#endif
!
   if (associated(mumps%A_loc))   deallocate(mumps%A_loc)
   if (associated(mumps%IRN_loc)) deallocate(mumps%IRN_loc)
   if (associated(mumps%JCN_loc)) deallocate(mumps%JCN_loc)
   if (associated(mumps%RHS))     deallocate(mumps%RHS)
!
!..Destroy the instance (deallocate internal data structures)
   mumps%JOB = -2
!
#if C_MODE
   call zmumps(mumps)
#else
   call dmumps(mumps)
#endif
!
end subroutine mumps_destroy
!
! -----------------------------------------------------------------------
!
subroutine mumps_start_subd
!
!..NOTE:
!..this routine initiates sparse solver for interior of subdomain
!
!..Define a communicator for the package.
   mumps_bub%COMM = MPI_COMM_SELF
!
!..PAR
!     0 : host is not involved in factorization/solve phases
!     1 : host involved in factorization/solve phases
   mumps_bub%PAR = 1
!
!..SYM (for LU factorization, no Hermitian option in version 5.2.1)
!     0: unsymmetric
!     1: symmetric positive definite
!     2: symmetric
   mumps_bub%SYM = 0
!
!..initialize an instance of the package
   mumps_bub%JOB = -1
!
#if C_MODE
   call zmumps(mumps_bub)
#else
   call dmumps(mumps_bub)
#endif
!
!..diagnostics
   !write(*,100) 'mumps_bub: MYID           = ', mumps_bub%MYID
   !write(*,110) 'mumps_bub: VERSION_NUMBER = ', mumps_bub%VERSION_NUMBER
 100 format(A,I3)
 110 format(A,A)
!
!..output stream for error messages
!     0: suppress messages
!     6: stdout
   mumps_bub%icntl(1) = 6
!
!..output stream for diagnostic/statistics/warning messages
!     0: suppress messages
!     6: stdout
   mumps_bub%icntl(2) = 0
!
!..output stream for global information
!     0: suppress messages
!     6: stdout
   mumps_bub%icntl(3) = 0
!
!..print level for error/warning/diagnostic messages
!     0: no output messages
!     1: only error messages
!     2: errors, warnings, and main stats
!     3: errors, warnings, and terse diagnostics
!     4: errors, warnings, and info on input/output params
   mumps_bub%icntl(4) = 1
!
!..icntl(5): matrix input format
!     0: assembled input format
!     1: elemental input format
   mumps_bub%icntl(5) = 0
!
!..icntl(7): choice of sequential pivot ordering tool (not used if icntl(28)=2)
!     0: Approximate Minimum Degree (AMD)
!     1: pivot order set by the user
!     2: Approximate Minimum Fill (AMF)
!     3: Scotch
!     4: PORD
!     5: Metis
!     6: AMD w/ quasi-dense row detection
!     7: automatic value
   mumps_bub%icntl(7) = 5
!
!..icntl(14): percentage increase in estimated workspace
!     [default: 20] - 20% increase in workspace
   mumps_bub%icntl(14) = 30
!
end subroutine mumps_start_subd
!
! -----------------------------------------------------------------------
!
subroutine mumps_destroy_subd
!
   if (associated(mumps_bub%A))   deallocate(mumps_bub%A)
   if (associated(mumps_bub%IRN)) deallocate(mumps_bub%IRN)
   if (associated(mumps_bub%JCN)) deallocate(mumps_bub%JCN)
   if (associated(mumps_bub%RHS)) deallocate(mumps_bub%RHS)
!
!..Destroy the instance (deallocate internal data structures)
   mumps_bub%JOB = -2
!
#if C_MODE
   call zmumps(mumps_bub)
#else
   call dmumps(mumps_bub)
#endif
!
end subroutine mumps_destroy_subd
!
end module par_mumps
