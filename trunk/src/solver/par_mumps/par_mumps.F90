! -----------------------------------------------------------------------
!
!    module name        - mumps
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 2019
!
!    purpose            - module sets up required workspace for 
!                         interfacing with distributed MUMPS solver
!
! -----------------------------------------------------------------------
#include "implicit_none.h"
module par_mumps
!
   use MPI      , only: MPI_COMM_SELF,MPI_COMM_WORLD
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
!
! -----------------------------------------------------------------------
!
! 
subroutine mumps_start
!
!..NOTE:
!..If you plan to run MUMPS sequentially from a parallel MPI application,
!  you need to install the parallel version of MUMPS and pass a
!  communicator containing a single processor to the MUMPS library
!
!..Define a communicator for the package.
   if (NUM_PROCS > 1) then
!  ...this assumes that this program is only executed by the master
      mumps_par%COMM = MPI_COMM_SELF
   else
      mumps_par%COMM = MPI_COMM_WORLD
   endif
!
!..1 => host involved in factorization/solve phases, 0 otherwise
   mumps_par%PAR = 1
!..Initialize an instance of the package
   mumps_par%JOB = -1
!..for L U factorization (sym = 0, with working host)
   mumps_par%SYM = 0
!
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
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
!..1 => element input format, 0 => assembled input format
   mumps_par%icntl(5) = 0
! 
!..corresponds to the percentage increase in the estimated working space
   ! mumps_par%icntl(14)  = 100
!
!  0 => Approximate Minimum Degree (AMD)
!  2 => Approximate Minimum Fill (AMF)
!  4 => PORD
!  5 => METIS
!  6 => AMD w/ quasi-dense row detection
!  7 => automatic value
   mumps_par%icntl(7) = 5
!
end subroutine mumps_start
! 
!
! -----------------------------------------------------------------------
!
! 
subroutine mumps_destroy
!
   if (mumps_par%MYID .eq. 0) then
      deallocate(mumps_par%IRN)
      deallocate(mumps_par%JCN)
      deallocate(mumps_par%A)
      deallocate(mumps_par%RHS)
   endif
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
   !call mpi_finalize(IERR)
! 
end subroutine mumps_destroy
! 
! 
end module par_mumps
