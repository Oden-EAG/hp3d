!
#include "implicit_none.h"
! -----------------------------------------------------------------------
!
!    module name        - mumps
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 2019
!
!    purpose            - module sets up required workspace for 
!                         interfacing with OpenMP MUMPS solver
!
! -----------------------------------------------------------------------
module mumps
!
   use MPI, only: MPI_COMM_SELF
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
   mumps_par%COMM = MPI_COMM_SELF
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
!..output stream for error messages
!     0: suppress messages
!     6: stdout
   mumps_par%icntl(1) = 6
!
!..output stream for diagnostic/statistics/warning messages
!     0: suppress messages
!     6: stdout
   mumps_par%icntl(2) = 0
!
!..output stream for global information
!     0: suppress messages
!     6: stdout
   mumps_par%icntl(3) = 0
!
!..print level for error/warning/diagnostic messages
!     0: no output messages
!     1: only error messages
!     2: errors, warnings, and main stats
!     3: errors, warnings, and terse diagnostics
!     4: errors, warnings, and info on input/output params
   mumps_par%icntl(4) = 1
!
!..icntl(5): matrix input format
!     0: assembled input format
!     1: elemental input format
   mumps_par%icntl(5) = 0
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
   mumps_par%icntl(7) = 5
!
!..icntl(14): percentage increase in estimated workspace
!     [default: 20] - 20% increase in workspace
   mumps_par%icntl(14) = 30
!
end subroutine mumps_start
!
! -----------------------------------------------------------------------
!
subroutine mumps_destroy
!
   if (associated(mumps_par%IRN)) deallocate(mumps_par%IRN)
   if (associated(mumps_par%JCN)) deallocate(mumps_par%JCN)
   if (associated(mumps_par%A  )) deallocate(mumps_par%A)
   if (associated(mumps_par%RHS)) deallocate(mumps_par%RHS)
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
end module mumps
