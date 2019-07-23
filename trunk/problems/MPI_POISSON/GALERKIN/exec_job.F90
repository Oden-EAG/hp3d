!----------------------------------------------------------------------
! exec_job
!----------------------------------------------------------------------
subroutine exec_job
!
   use common_prob_data
   use data_structure3D
   use MPI           , only: MPI_COMM_WORLD
   use mpi_param     , only: RANK,ROOT
   use par_mesh      , only: EXCHANGE_DOF,distr_mesh
   use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_eval
!
   implicit none
!
   integer :: i,imax,ierr
   real(8) :: MPI_Wtime,start_time,end_time
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
!..distribute mesh initially
   call distr_mesh
!..set Zoltan partitioner
   call zoltan_w_set_lb(1)
!
!..set number of iterations
   imax = 4
!
   do i=1,imax
!
      if(RANK .eq. ROOT) write(*,100) 'Beginning iteration i = ', i
!
!  ...single uniform h-refinement
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '1. global h-refinement...'
      call uniform_href(IUNIFORM,1,0.25d0)
      if (i .le. 2) cycle
!
!  ...distribute mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '2. distributing mesh...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call distr_mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
!  ...print current partition (elems)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if (i .le. 2) then
         if(RANK .eq. ROOT) write(*,200) '3. printing current partition (elems)...'
         call print_partition
      endif
!
!  ...evaluate current partition
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '4. evaluating current partition...'
      call zoltan_w_eval
!
!  ...run mesh verification routines
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '5. verify distributed mesh consistency...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call par_verify
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
!  ...solve problem with par_mumps (MPI MUMPS)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '6. calling MUMPS (MPI) solver...'
      call par_mumps_sc('G')
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
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
