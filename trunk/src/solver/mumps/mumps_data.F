c------------------------------------------------------------------
c
c   module name        - mumps_data
c
c------------------------------------------------------------------
c
c   latest revision    - June 11
c
c   purpose            - module stores a [dz]mumps_struc object and
c                        interface routines to control MUMPS
c
c-----------------------------------------------------------------------
c
      module mumps_data

      use parallelism
      use data_structure3D, ONLY: MAXNODS
      use solvermod
c
c     use non-preprocessed include statement since preprocessor objects
c     to un-paired quotation marks
#if C_MODE
#include "zmumps_struc.h"
#else
#include "dmumps_struc.h"
#endif
c
      save
c
#if C_MODE
      type(zmumps_struc) :: mumps_par
#else
      type(dmumps_struc) :: mumps_par
#endif
c
c
      parameter (MAX_PROC_OWNERS=10)

      type node_mumps_data_type
        integer :: visited
      endtype node_mumps_data_type

      type(node_mumps_data_type),allocatable :: NODE_MUMPS_DATA(:)
c-> MP 21.4.2014
      type(node_mumps_data_type),allocatable :: NODE_MUMPS_DATA_H1(:)
      type(node_mumps_data_type),allocatable :: NODE_MUMPS_DATA_Hdiv(:)
      type(node_mumps_data_type),allocatable :: NODE_MUMPS_DATA_Hcurl(:)
      type(node_mumps_data_type),allocatable :: NODE_MUMPS_DATA_L2(:)
c<- MP 21.4.2014

      parameter (ZERO_MUMPS = 1.d-16)

      contains
c
c
c-----------------------------------------------------------------------
c INITIALIZE MUMPS FOR THE FIRST TIME (BEGINNING OF MAIN)
c-----------------------------------------------------------------------
      subroutine start_mumps
      use control, only: ISYM_FLAG, SYMMETRIC, UNSYMMETRIC
#include "syscom.blk"
#include "surfsc1.blk"
#include "mpif.h"
c

#if PARALLEL_MODE
#else
! here we assume that PetscInitialize was called during startup, which initializes MPI,
! therefore skip it here:
#ifndef __PETSC_USE__
c  ...initialize MPI
      !call mpi_init(ierr)
#endif
#endif
c
c  ...set an MPI communicator for all MUMPS processors
      mumps_par%comm = MPI_COMM_WORLD
c
c  ...indicate that host processor should participate in the
c     factorization and solve phases (this is a sequential code)
      mumps_par%par = 1
c
c     unsymmetric (0) , symmetric positive definite (1), symmetric (2)
      select case (ISYM_FLAG)
      case(SYMMETRIC);    mumps_par%sym = 2
      case(UNSYMMETRIC);  mumps_par%sym = 0
      end select
c
      end subroutine start_mumps
c
c
c-----------------------------------------------------------------------
c    ALLOCATE MEMORY AND INITIALIZE MUMPS
c-----------------------------------------------------------------------
      subroutine alloc_mumps

c
c  ...initialize an instance of MUMPS
      mumps_par%job = -1
#if C_MODE
      write(*,*)'call zmumps, job=-1'
      call zmumps(mumps_par)
#else
      write(*,*)'call dmumps, job=-1'
      call dmumps(mumps_par)
#endif
c
c  ...set control flags
#if C_MODE
      mumps_par%cntl(1) = 0.0
c     error output stream (non-positive to suppress)
      mumps_par%icntl(1)  = 1 !1
c     diagnostic, statistics and warnings
      mumps_par%icntl(2)  = 1 !1
c     global information
      mumps_par%icntl(3)  = 6 !6
c     printing level
      mumps_par%icntl(4)  = 3 !3
c     input matrix in element format (1) or assembled format (0)
      mumps_par%icntl(5)  = 0
c     column permutation for zero-free diagonal (automatic)
      mumps_par%icntl(6)  = 0
c     pivot order (automatic)
      mumps_par%icntl(7)  = 5
c     scaling (automatic)
      mumps_par%icntl(8)  = 0
c     no transpose
      mumps_par%icntl(9)  = 1
c     max steps for iterative refinement
      mumps_par%icntl(10) = 0
c     statistics info
      mumps_par%icntl(11) = 0
c     controls parallelism
      mumps_par%icntl(12) = 1
c     use ScaLAPACK for root node
      mumps_par%icntl(13) = 0
c     percentage increase in estimated workspace
      mumps_par%icntl(14) = 100
c
c     matrix distribution for assembled input
#if PARALLEL_MODE
      mumps_par%icntl(18) = 3 !0 when submit from host 3 when distributed
#else
      mumps_par%icntl(18) = 0
#endif
c     nonzero for Schur complement
      mumps_par%icntl(19) = 0
#else
c     error output stream (non-positive to suppress)
      mumps_par%icntl(1)  = 1 !1
c     diagnostic, statistics and warnings
      mumps_par%icntl(2)  = 1 !1
c     global information
      mumps_par%icntl(3)  = 6 !6
c     printing level
      mumps_par%icntl(4)  = 3 !3
c     input matrix in element (1) or assembled (0) format
      mumps_par%icntl(5)  = 0
c     column permutation for zero-free diagonal (automatic)
      mumps_par%icntl(6)  = 0
c     pivot order (automatic)
      mumps_par%icntl(7)  = 5
c     scaling (automatic)
      mumps_par%icntl(8)  = 0
c     no transpose
      mumps_par%icntl(9)  = 1
c     max steps for iterative refinement
c     mumps_par%icntl(10) = 0
c     statistics info
      mumps_par%icntl(11) = 0
c     controls parallelism
      mumps_par%icntl(12) = 1
c     use ScaLAPACK for root node
c     mumps_par%icntl(13) = 0
c     percentage increase in estimated workspace
      mumps_par%icntl(14) = 100
c
c     matrix distribution for assembled input
#if PARALLEL_MODE
      mumps_par%icntl(18) = 3 ! 0 host 3 distributed storage
#else
      mumps_par%icntl(18) = 0
#endif
c     nonzero for Schur complement
c     mumps_par%icntl(19) = 0
#endif
c
      !mumps_par%icntl(22) = 1 ! out of core
      end subroutine alloc_mumps
c
c
c-----------------------------------------------------------------------
c    RUN MUMPS SOLVE (LU FACTORIZATION)
c-----------------------------------------------------------------------
      subroutine run_mumps_solve
#include "syscom.blk"
#include "mpif.h"
c
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call start_clock(iclock)
c
c  ...analysis phase
      mumps_par%job = 1
#if C_MODE
      write(*,*)'call zmumps, job=1'
      call zmumps(mumps_par)
#else
      write(*,*)'call dmumps, job=1'
      call dmumps(mumps_par)
#endif

c  ...check for errors
      if (mumps_par%info(1).ne.0) then
        write (*,*) MYRANK,'mumps_par%job=',mumps_par%job
        write (*,*) MYRANK,'mumps_par%info=',mumps_par%info
        if(mumps_par%info(1).lt.0)stop 1
      endif

      ifirst=0
 222  continue
      call stop_clock(dtime,iclock)
      ANALYSIS_TIME = dtime
      if(ifirst.eq.0)write(*,*)MYRANK,'ANALYSIS TIME:',dtime
      ifirst=1
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call start_clock(iclock)
c
c
c  ...factorization phase
      mumps_par%job = 2
      write(*,*)'run_mumps_solve230:mumps_par%NZ',mumps_par%NZ
#if C_MODE
      call zmumps(mumps_par)
#else
      call dmumps(mumps_par)
#endif
c
c  ...if there is an error (stop the code)
      if (mumps_par%info(1).ne.0) return
c
c  ...time and memory
      call stop_clock(dtime,iclock)
      FACTOR_TIME = dtime
      imemory=mumps_par%info(16)
      write(*,*)MYRANK,'FACTORIZATRION',dtime,
     .  'sec; memory',imemory,' Mb'

#if PARALLEL_MODE
      call mpi_allreduce(dtime,dtime_out,
     .     1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      dtime=dtime_out
      call mpi_allreduce(imemory,imemory_out,
     .     1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      imemory=imemory_out
#endif
      if (MYRANK.eq.0)
     . write (*,*) 'MAX MUMPS:',dtime,' sec',imemory,' Mb'
c
      if (mumps_par%icntl(11)>0) then
c
c       print statistical information
        write(*,1010) mumps_par%rinfog(10)
        write(*,1011) mumps_par%rinfog(11)
 1010   format('cond1(A) = ',e15.8)
 1011   format('cond2(A) = ',e15.8)
      endif
c
      end subroutine run_mumps_solve
c
c-----------------------------------------------------------------------
c     BACKWARD ELIMINATION
c-----------------------------------------------------------------------

      subroutine run_mumps_rhs
#include "syscom.blk"

      ifirst=0
 223  continue
      if(ifirst.eq.1)call stop_clock(dtme,iclock_rhs)
      call start_clock(iclock_rhs)
#if C_MODE
      mumps_par%job = 3
      call zmumps(mumps_par)
#else
      mumps_par%job = 3
      call dmumps(mumps_par)
#endif

c
c  ...if there is an error (or stop the code)
      if (mumps_par%info(1).ne.0) return
c
      call stop_clock(dtime,iclock_rhs)
      write(*,*)MYRANK,'SOLUTION',dtime
      SOLVE_TIME = dtime
      end subroutine run_mumps_rhs
c
c-----------------------------------------------------------------------
c    TERMINATE MUMPS
c-----------------------------------------------------------------------
      subroutine dealloc_mumps
c

c  ...terminate an instance of MUMPS
      mumps_par%job = -2
#if C_MODE
      call zmumps(mumps_par)
#else
      call dmumps(mumps_par)
#endif
c
      end subroutine dealloc_mumps
c
c
      end module mumps_data
