!----------------------------------------------------------------------
! exec_job
!----------------------------------------------------------------------
subroutine exec_job
!
   use common_prob_data
   use data_structure3D
   use mpi_wrapper
   use par_mesh      , only: EXCHANGE_DOF,distr_mesh
   use zoltan_wrapper
!
   implicit none
!
   integer :: i,ierr
   real(8) :: start_time,end_time
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
!..distribute mesh initially
   call distr_mesh
!..set Zoltan partitioner
   call zoltan_w_set_lb(ZOLTAN_LB_DEFAULT)
!
   do i=1,IMAX
!
      if(RANK .eq. ROOT) write(*,100) 'Beginning iteration i = ', i
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if (i .le. IMAX) then
!     ...single uniform h-refinement
         if(RANK .eq. ROOT) write(*,200) '1. global uniform h-refinement...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call uniform_href(IUNIFORM,1,0.25d0)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      else
!     ...single anisotropic (in z) h-refinement
         if(RANK .eq. ROOT) write(*,200) '1. global anisotropic h-refinement...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call uniform_href(IUNIFORM,2,0.25d0)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      endif
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
      if (NUM_PROCS .eq. 1) goto 30
!
!      if (i .eq. IMAX-5) then
!         call zoltan_w_set_lb(ZOLTAN_LB_FIBER)
!      elseif (i .gt. IMAX-5) then
!         goto 30
!      endif
!  ...distribute mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '2. distributing mesh...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call distr_mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   30 continue
      !if (i .le. IMAX-2) cycle
!
!  ...print current partition (elems)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '3. printing current partition (elems)...'
      if (i .le. 2) then
         call print_partition
      else
         if(RANK .eq. ROOT) write(*,*) '   ... skipping for a large number of elements.'
      endif
!
      if (NUM_PROCS .eq. 1) goto 50
      goto 50
!
!  ...evaluate current partition
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '4. evaluating current partition...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call zoltan_w_eval
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   50 continue
      goto 60
!
!  ...run mesh verification routines
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '5. verify distributed mesh consistency...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call par_verify
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   60 continue
!
!  ...solve problem with par_mumps (MPI MUMPS)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
      if(RANK .eq. ROOT) write(*,200) '6. calling solver...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call par_mumps_sc('G')
      !call par_nested('G')
      !call petsc_solve('G')
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   enddo
!
!..compute error on last mesh
   if(RANK .eq. ROOT) write(*,200) '7. compute error on last mesh...'
   call uniform_href(INOREFINEMENT,1,0.25d0)
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
