!------------------------------------------------------------------------------
!> @brief      Example routine for automating code (no input). Does a number
!!             of solves and refinements; repartitioning in between
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine exec_job
!
      use commonParam
      use data_structure3D
      use MPI,             only: MPI_COMM_WORLD
      use mpi_param
      use par_mesh,        only: EXCHANGE_DOF,distr_mesh
      use zoltan_wrapper,  only: zoltan_w_set_lb,zoltan_w_eval
!
      implicit none
!
      integer :: flag(2),iParAttr(2)
      integer :: physNick, nstop
!
      integer :: i, ierr, numPts, fld
      real(8) :: MPI_Wtime, start_time, end_time
!
      logical :: adaptive_refs = .false.
!
!----------------------------------------------------------------------
!
      EXCHANGE_DOF = .false.
!
      if(RANK .eq. ROOT) then
         write(*,*) '=================='
         write(*,*) 'exec_job: starting'
         write(*,*) '=================='
      endif
!
      if (NUM_PROCS.gt.1) then
!     ...distribute mesh initially
         call distr_mesh
!     ...set Zoltan graph-based partitioner
         call zoltan_w_set_lb(6)
      endif
!
      do i=1,IMAX
!
         if (RANK .eq. ROOT) write(*,*) 'Beginning iteration i = ', i
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr);
         if (adaptive_refs) then
            call refine_DPG(IADAPTIVE,1,0.25d0, nstop)
         else
            call refine_DPG(IUNIFORM,1,0.25d0, nstop)
         endif
!
!     ...set partitioner for load balancing
         if (NUM_PROCS .gt. 1) then
            call distr_mesh
         endif
!
!     ...solve problem with par_mumps (MPI MUMPS)
         if (NUM_PROCS .eq. 1) then
            call par_mumps_sc('H')
         else
            call par_mumps_sc('H')
         endif
!
!     ...write paraview output
         iParAttr = (/2,6/)
         call my_paraview_driver(iParAttr)
!
      enddo
!
!  ...compute residual/error on last mesh
      call refine_DPG(INOREFINEMENT,1,0.25d0, nstop)
!
      if(RANK .eq. ROOT) then
         write(*,*)
         write(*,*) '=================='
         write(*,*) 'exec_job: finished'
         write(*,*) '=================='
         write(*,*)
      endif
!
end subroutine exec_job
