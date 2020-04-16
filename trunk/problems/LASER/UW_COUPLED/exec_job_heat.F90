!----------------------------------------------------------------------
! exec_job: Linear Heat problem simulation
!----------------------------------------------------------------------
subroutine exec_job_heat
!
   use commonParam
   use control
   use data_structure3D
   use MPI           , only: MPI_COMM_WORLD
   use mpi_param     , only: RANK,ROOT,NUM_PROCS
   use par_mesh      , only: EXCHANGE_DOF,distr_mesh
   use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_eval
!
   implicit none
!
   integer :: flag(6)
   integer :: physNick,nstop
   logical :: ires
!
   integer :: i,ierr,numPts,fld
   real(8) :: MPI_Wtime,start_time,end_time
!
!----------------------------------------------------------------------
!
   if(NEXACT.ne.1) then
      write(*,*) 'NEXACT must be 1 for linear Heat problem. stop.'
      stop
   endif
!
   EXCHANGE_DOF = .false.
!
   NO_PROBLEM = 1
   call set_physAm(NO_PROBLEM, physNick,flag)
   ires = .true.
!
   if(RANK .eq. ROOT) then
      write(*,*) '======================='
      write(*,*) 'exec_job_heat: starting'
      write(*,*) '======================='
   endif
!
!..distribute mesh initially
   call distr_mesh
!..set Zoltan partitioner
   call zoltan_w_set_lb(0)
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
         call refine_DPG(IUNIFORM,1,0.25d0,flag,physNick,ires, nstop)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      else
!     ...single anisotropic (in z) h-refinement
         if(RANK .eq. ROOT) write(*,200) '1. global anisotropic h-refinement...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call refine_DPG(IANISOTROPIC,1,0.25d0,flag,physNick,ires, nstop)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      endif
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
      if (NUM_PROCS .eq. 1) goto 30
!
      if (i .eq. IMAX-3) then
         call zoltan_w_set_lb(7)
      elseif (i .gt. IMAX-3) then
         goto 30
      endif
!
!  ...distribute mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '2. distributing mesh...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call distr_mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   30 continue
      !if (i .le. IMAX-1) cycle
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
      if(RANK .eq. ROOT) write(*,200) '6. calling MUMPS (MPI) solver...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      if (NUM_PROCS .eq. 1) then
         call pardiso_sc('H')
      else
         !call par_mumps_sc('H')
         call par_nested('H')
         !call par_hybrid('H')
      endif
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   enddo
!
!..compute error on last mesh
   if(RANK .eq. ROOT) write(*,200) '7. compute error on last mesh...'
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   call refine_DPG(INOREFINEMENT,1,0.25d0,flag,physNick,ires, nstop)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!
  100 format(/,'/////////////////////////////////////////////////////////////', &
             /,'             ',A,I2,/)
  200 format(/,'-------------------------------------------------------------', &
             /,'  --> ',A,/)
  300 format(' timer: ',f12.5,' seconds')
!
   if(RANK .eq. ROOT) then
      write(*,*)
      write(*,*) '======================='
      write(*,*) 'exec_job_heat: finished'
      write(*,*) '======================='
      write(*,*)
   endif
!
end subroutine exec_job_heat
