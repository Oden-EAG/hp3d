!----------------------------------------------------------------------
!> @brief       Simple job script for solving problem and outputting
!!              solution to paraview
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine exec_job
!
      use common_prob_data_UW
      use data_structure3D
      use mpi_param
      use MPI,             only: MPI_COMM_WORLD
      use par_mesh,        only: EXCHANGE_DOF, distr_mesh
      use paraview,        only: VLEVEL
      use zoltan_wrapper,  only: zoltan_w_set_lb,zoltan_w_eval,zz,zoltan_w_handle_err
      use zoltan
!
      implicit none
!
      integer :: i, ierr
      integer :: iParAttr(3)
      real(8) :: MPI_Wtime, start_time, end_time
!
!----------------------------------------------------------------------
!
      EXCHANGE_DOF = .true.
!
      if(RANK .eq. ROOT) then
         write(*,*) '=================='
         write(*,*) 'exec_job: starting'
         write(*,*) '=================='
      endif
!
!  ...h-refine so all procs can have elements
      if (NUM_PROCS .ge. NRELES) then
         if (RANK .eq. ROOT) write(*,*) 'global h-refinement...'
         call global_href
         if (IBC_PROB.eq.3 .or. IBC_PROB.eq.4 .or. IBC_PROB.eq.6) call propagate_flag(2,3)
         call update_gdof
         call update_Ddof
      endif
!
!  ...distribute mesh using graph-based partitioner
      if (NUM_PROCS.gt.1) then
!     ...initial distribution
         call distr_mesh
!
!     ...set graph-based partitioner
         call zoltan_w_set_lb(6)
         ierr = Zoltan_Set_Param(zz,'LB_APPROACH','PARTITION')
         call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
!     ...execute graph-based partition
         call distr_mesh
!
!     ...set graph algorithm to less-expensive `repartition` (for future repartitioning)
         ierr = Zoltan_Set_Param(zz,'LB_APPROACH','REPARTITION')
         call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
      endif
!
!  ...call MUMPS
      call par_mumps_sc('H')
!
!  ...set which attributes to output
      iParAttr = (/1, 1, 4/)
!
!  ...write paraview output
      call my_paraview_driver(iParAttr)
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
