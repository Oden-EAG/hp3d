!----------------------------------------------------------------------
! exec_job
!----------------------------------------------------------------------
subroutine exec_job
!
   use common_prob_data_UW
   use data_structure3D
   use MPI           , only: MPI_COMM_WORLD
   use mpi_param     , only: RANK,ROOT,NUM_PROCS
   use par_mesh      , only: EXCHANGE_DOF,distr_mesh
   use paraview      , only: VLEVEL
   use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_eval,zz,zoltan_w_handle_err
   use zoltan
!
   implicit none
!
   integer :: i,ierr
   real(8) :: MPI_Wtime,start_time,end_time
   integer :: iParAttr(3)
   character(len=2) :: vis_level
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
!..h-refine so all procs can have elements
   if (NUM_PROCS .ge. 8) then
      if (RANK .eq. ROOT) write(*,*) 'global h-refinement...'
      call global_href
      if (IBC_PROB.eq.3 .or. IBC_PROB.eq.4 .or. IBC_PROB.eq.6) call propagate_flag(2,3)
      call update_gdof
      call update_Ddof
   endif
!
!..distribute mesh
   if (NUM_PROCS.gt.1) then
      call distr_mesh
      call zoltan_w_set_lb(6)
!
      ierr = Zoltan_Set_Param(zz,'LB_APPROACH','PARTITION')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
      call distr_mesh
!
      ierr = Zoltan_Set_Param(zz,'LB_APPROACH','REPARTITION')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
   endif
!
   iParAttr = (/1, 1, 4/)
!
!..call Multigrid solver
   call mg_driver(1,1,0.2d0,1,20,50,20,1.d-4,0.2d0,.false.,iParAttr)
!
!..write paraview output
!   if(RANK .eq. ROOT) write(*,200) '8. writing paraview output...'
!   VLEVEL='2'
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!   call my_paraview_driver(iParAttr)
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!   if(RANK .eq. ROOT) write(*,300) end_time - start_time
!




!   do i=1,IMAX
!!
!      if(RANK .eq. ROOT) write(*,100) 'Beginning iteration i = ', i
!!
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!      if (i .le. IMAX) then
!!     ...single uniform h-refinement
!         if(RANK .eq. ROOT) write(*,200) '1. global uniform h-refinement...'
!         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!         call uniform_href(IUNIFORM,1,0.25d0)
!         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!      else
!!     ...single anisotropic (in z) h-refinement
!         if(RANK .eq. ROOT) write(*,200) '1. global anisotropic h-refinement...'
!         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!         call uniform_href(IUNIFORM,2,0.25d0)
!         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!      endif
!      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!!
!      if (NUM_PROCS .eq. 1) goto 30
!!
!!      if (i .eq. IMAX-5) then
!!         call zoltan_w_set_lb(7)
!!      elseif (i .gt. IMAX-5) then
!!         goto 30
!!      endif
!!  ...distribute mesh
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!      if(RANK .eq. ROOT) write(*,200) '2. distributing mesh...'
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!      call distr_mesh
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!!
!   30 continue
!      !if (i .le. IMAX-2) cycle
!!
!!  ...print current partition (elems)
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!      if(RANK .eq. ROOT) write(*,200) '3. printing current partition (elems)...'
!      if (i .le. 2) then
!         call print_partition
!      else
!         if(RANK .eq. ROOT) write(*,*) '   ... skipping for a large number of elements.'
!      endif
!!
!      if (NUM_PROCS .eq. 1) goto 50
!      goto 50
!!
!!  ...evaluate current partition
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!      if(RANK .eq. ROOT) write(*,200) '4. evaluating current partition...'
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!      call zoltan_w_eval
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!!
!   50 continue
!      goto 60
!!
!!  ...run mesh verification routines
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!      if(RANK .eq. ROOT) write(*,200) '5. verify distributed mesh consistency...'
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!      call par_verify
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!!
!   60 continue
!!
!!  ...solve problem with par_mumps (MPI MUMPS)
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!      if(RANK .eq. ROOT) write(*,200) '6. calling PETSc (MPI) solver...'
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!      !call par_mumps_sc('G')
!      !call par_nested('G')
!      call petsc_solve('G')
!      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!!
!   enddo
!!
!!..compute error on last mesh
!   if(RANK .eq. ROOT) write(*,200) '7. compute error on last mesh...'
!   call uniform_href(INOREFINEMENT,1,0.25d0)
!
  100 format(/,'/////////////////////////////////////////////////////////////', &
             /,'             ',A,I2,/)
  200 format(/,'-------------------------------------------------------------', &
             /,'  --> ',A,/)
  300 format(' timer: ',f12.5,' seconds')
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
