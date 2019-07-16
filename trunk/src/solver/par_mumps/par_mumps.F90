!
#include "implicit_none.h"
! -----------------------------------------------------------------------
!
!    module name        - par_mumps
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 2019
!
!    purpose            - module sets up required workspace for 
!                         interfacing with distributed MUMPS solver
!
! -----------------------------------------------------------------------
module par_mumps
!
   use MPI      , only: MPI_COMM_WORLD
   use MPI_param, only: RANK,NUM_PROCS
!
   implicit none
!
#if C_MODE
   include 'zmumps_struc.h'
   type (ZMUMPS_STRUC) mumps_par
#else
   include 'dmumps_struc.h'
   type (DMUMPS_STRUC) mumps_par
#endif
!
   contains
!
! -----------------------------------------------------------------------
!
subroutine mumps_start
!
!..NOTE:
!..If you plan to run MUMPS sequentially from a parallel MPI application,
!  you need to install the parallel version of MUMPS and pass a
!  communicator containing a single processor to the MUMPS library
!
!..Define a communicator for the package.
   mumps_par%COMM = MPI_COMM_WORLD
!
!..PAR
!     0 : host is not involved in factorization/solve phases
!     1 : host involved in factorization/solve phases
   mumps_par%PAR = 1
!
!..SYM (for LU factorization, no Hermitian option in version 5.2.1)
!     0: unsymmetric
!     1: symmetric positive definite
!     2: symmetric
   mumps_par%SYM = 0
!
!..initialize an instance of the package
   mumps_par%JOB = -1
!
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
!
!..verify MYID == RANK
   if (mumps_par%MYID .ne. RANK) then
      write(*,*) 'par_mumps: MYID .ne. RANK'
      call pause
   endif
!
!..diagnostics
   write(*,100) 'par_mumps: MYID           = ', mumps_par%MYID
   write(*,110) 'par_mumps: VERSION_NUMBER = ', mumps_par%VERSION_NUMBER
 100 format(A,I3)
 110 format(A,A)
!
!..output for error messages
   mumps_par%icntl(1) = 6
!
!..output for diagnostic/statistics/warning messages
   mumps_par%icntl(2) = 0
!
!..output for global information
   mumps_par%icntl(3) = 0
!
!..print level for error/warning/diagnostic messages
   mumps_par%icntl(4) = 6
!
!..icntl(5): matrix input format
!     0: assembled input format
!     1: elemental input format
   mumps_par%icntl(5) = 0
!
!..icntl(18): distribution strategy of the input matrix
!     0: centralized on host
!     1: user provides matrix structure at analysis, MUMPS returns mapping for entries
!     2: user provides matrix structure at analysis, entries at factorization
!     3: user provides distributed matrix
   mumps_par%icntl(18) = 3
!
!..icntl(21): determines distribution of the solution vectors
!     0: centralized on host
!     1: distributed
   mumps_par%icntl(21) = 0
!
!..icntl(28): sequential or parallel analysis
!     0: automatic choice
!     1: sequential computation
!     2: parallel computation
   mumps_par%icntl(28) = 2
!
!..icntl(29): choice of parallel pivot ordering tool (not used if icntl(28)=1)
!     0: automatic choice
!     1: PT-Scotch
!     2: ParMetis
   mumps_par%icntl(29) = 2
!
end subroutine mumps_start
!
! -----------------------------------------------------------------------
!
subroutine mumps_destroy
!
   if (associated(mumps_par%IRN_loc)) deallocate(mumps_par%IRN_loc)
   if (associated(mumps_par%JCN_loc)) deallocate(mumps_par%JCN_loc)
   if (associated(mumps_par%A_loc))   deallocate(mumps_par%A_loc)
   if (associated(mumps_par%RHS))     deallocate(mumps_par%RHS)
!
!..Destroy the instance (deallocate internal data structures)
   mumps_par%JOB = -2
!
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
! 
end subroutine mumps_destroy
!
!
end module par_mumps
