!
#include "typedefs.h"
!
! -----------------------------------------------------------------------
!
!    routine name       - par_solve
!
! -----------------------------------------------------------------------
!
!    latest revision    - Mar 2020
!
!    purpose            - interface for distributed MUMPS solver
!                       - routine perform distributed parallel solve
!                         of the MUMPS instance passed to the input
!
! -----------------------------------------------------------------------
subroutine par_solve(mumps)
!
   use assembly_sc, only: IPRINT_TIME
   use mpi_wrapper
   use par_mumps
!
   implicit none
!
#if HP3D_COMPLEX
   type (ZMUMPS_STRUC), intent(inout) :: mumps
#else
   type (DMUMPS_STRUC), intent(inout) :: mumps
#endif
!
!..MPI communicator
#if HP3D_USE_MPI_F08
   type(MPI_Comm) :: mumps_comm
#else
   integer        :: mumps_comm
#endif
!
!..aux variables
   integer :: ierr,mRANK,mNUM_PROCS
!
!..timer
   real(8) :: time_stamp,start_time
!
!..info (verbose output if true)
   logical :: info = .false.
!
! -----------------------------------------------------------------------
!
!..mumps%COMM uses F90 binding
#if HP3D_USE_MPI_F08
   mumps_comm%MPI_VAL = mumps%COMM
#else
   mumps_comm         = mumps%COMM
#endif
!
   call MPI_COMM_RANK(mumps_comm, mRANK     ,ierr)
   call MPI_COMM_SIZE(mumps_comm, mNUM_PROCS,ierr)
!
   if ((mRANK.eq.ROOT) .and. info) then
      write(*,5010) '[',RANK,'] par_solve: mNUM_PROCS = ', mNUM_PROCS
      write(*,*) ' - solving distributed sparse problem: mNUM_PROCS = ', mNUM_PROCS
      write(*,*) ' - mumps%icntl(28) = ', mumps%icntl(28)
      write(*,*) ' - mumps%icntl(29) = ', mumps%icntl(29)
 5010 format(A,I4,A,I4)
   endif
!
!..MUMPS analysis
   mumps%JOB = 1
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_comm, ierr)
      time_stamp = MPI_Wtime()
      start_time = time_stamp
   endif
!
#if HP3D_COMPLEX
   call zmumps(mumps)
#else
   call dmumps(mumps)
#endif
   if (mumps%INFOG(1) .ne. 0) then
      call mumps_destroy(mumps)
      if (mRANK.eq.ROOT) write(*,*) 'analysis: mumps%INFOG(1) .ne. 0'
      stop
   endif
   if (IPRINT_TIME .eq. 1 .and. info) then
      call MPI_BARRIER(mumps_comm, ierr)
      time_stamp = MPI_Wtime()-time_stamp
      if (mRANK .eq. ROOT) write(*,3001) time_stamp
 3001 format(' - Analysis : ',f12.5,'  seconds')
      if (mRANK .eq. ROOT .and. info) then
         write(*,1100) '   - MAX estimated size in GB = ',mumps%INFOG(16)/1000.d0
         write(*,1100) '   - SUM estimated size in GB = ',mumps%INFOG(17)/1000.d0
         write(*,1200) '   - Seq/parallel analysis    = ',mumps%INFOG(32)
         write(*,1200) '   - Ordering method used     = ',mumps%INFOG(7)
    1100 format(A,F11.3)
    1200 format(A,I1)
      endif
   endif
!
  55 continue
!
!..MUMPS factorization
   mumps%JOB = 2
!
   if (IPRINT_TIME .eq. 1 .and. info) then
      call MPI_BARRIER(mumps_comm, ierr)
      time_stamp = MPI_Wtime()
   endif
#if HP3D_COMPLEX
   call zmumps(mumps)
#else
   call dmumps(mumps)
#endif
   if (mumps%INFOG(1) .ne. 0) then
      if (mRANK.eq.ROOT) write(*,*) 'factorization: mumps%INFOG(1) .ne. 0'
      if (mumps%INFOG(1) .eq. -9) then
         if (mRANK.eq.ROOT) write(*,*) 'Increasing workspace, trying factorization again...'
         mumps%icntl(14) = mumps%icntl(14) + 30 ! increase workspace by 30 percent
         goto 55
      endif
      call mumps_destroy(mumps)
      stop
   endif
   if (IPRINT_TIME .eq. 1 .and. info) then
      call MPI_BARRIER(mumps_comm, ierr)
      time_stamp = MPI_Wtime()-time_stamp
      if (mRANK .eq. ROOT) write(*,3002) time_stamp
 3002 format(' - Factorize: ',f12.5,'  seconds')
      if (mRANK .eq. ROOT .and. info) then
         write(*,1100) '   - MAX memory used in GB    = ',mumps%INFOG(21)/1000.d0
         write(*,1100) '   - SUM memory used in GB    = ',mumps%INFOG(22)/1000.d0
      endif
   endif
!
!..MUMPS solve
   mumps%JOB = 3
!
  if (IPRINT_TIME .eq. 1 .and. info) then
      call MPI_BARRIER(mumps_comm, ierr)
      time_stamp = MPI_Wtime()
   endif
#if HP3D_COMPLEX
   call zmumps(mumps)
#else
   call dmumps(mumps)
#endif
   if (mumps%INFOG(1) .ne. 0) then
      call mumps_destroy(mumps)
      if (mRANK.eq.ROOT) write(*,*) 'solve: mumps%INFOG(1) .ne. 0'
      stop
   endif
   if (IPRINT_TIME .eq. 1 .and. info) then
      call MPI_BARRIER(mumps_comm, ierr)
      time_stamp = MPI_Wtime()-time_stamp
      if (mRANK .eq. ROOT) write(*,3003) time_stamp
 3003 format(' - Solve    : ',f12.5,'  seconds')
   endif
!
   if (IPRINT_TIME .eq. 1 .and. info) then
      call MPI_BARRIER(mumps_comm, ierr)
      time_stamp = MPI_Wtime()-start_time
      if (mRANK .eq. ROOT) write(*,3004) '[',RANK,'] par_solve: ',time_stamp,' seconds'
 3004 format(A,I4,A,f12.5,A)
   endif
!
end subroutine par_solve


























