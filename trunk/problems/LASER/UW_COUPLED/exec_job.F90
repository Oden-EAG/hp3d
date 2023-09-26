!----------------------------------------------------------------------
! exec_job: Linear Maxwell fiber simulation
!----------------------------------------------------------------------
subroutine exec_job
!
   use commonParam
   use data_structure3D
   use MPI           , only: MPI_COMM_WORLD,MPI_Wtime
   use mpi_param     , only: RANK,ROOT,NUM_PROCS
   use par_mesh      , only: EXCHANGE_DOF,distr_mesh
   use paraview      , only: paraview_select_attr
   use zoltan_wrapper
!
   implicit none
!
   integer :: flag(7)
   logical :: iPvAttr(7)
   integer :: physNick,nstop
   logical :: ires
!
   integer :: i,ierr,numPts,fld
   real(8) :: start_time,end_time
!
!----------------------------------------------------------------------
!
   EXCHANGE_DOF = .false.
!
   NO_PROBLEM = 3
   call set_physAm(NO_PROBLEM, physNick,flag)
   ires = .true.
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
   do i=1,IMAX+JMAX
!
      if(RANK .eq. ROOT) write(*,100) 'Beginning iteration i = ', i
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if (i .le. 0) then
!     ...single uniform h-refinement
         if(RANK .eq. ROOT) write(*,200) '1. global uniform h-refinement...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call refine_DPG(IUNIFORM,1,0.25d0,flag,physNick,ires, nstop)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      elseif  (i .le. IMAX) then
!     ...single anisotropic (in z) h-refinement
         if(RANK .eq. ROOT) write(*,200) '1. global anisotropic h-refinement...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call refine_DPG(IANISOTROPIC,1,0.25d0,flag,physNick,ires, nstop)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      else
!     ...adaptive h-refinement
         if(RANK .eq. ROOT) write(*,200) '1. adaptive h-refinement...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call refine_DPG(ICLAD,1,0.25d0,flag,physNick,ires, nstop)
         !call refine_DPG(IADAPTIVE,1,0.25d0,flag,physNick,ires, nstop)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      endif
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
      if (NUM_PROCS .eq. 1) goto 30
!
!  ...set partitioner for load balancing, redistributes mesh in 'distr_mesh'
      if (i .eq. IMAX-2) then
         call zoltan_w_set_lb(ZOLTAN_LB_FIBER) ! fiber partitioner
      elseif (i .gt. IMAX) then
         !call zoltan_w_set_lb(ZOLTAN_LB_GRAPH) ! graph partitioner
         goto 30 ! NO LOAD BALANCING
      else
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
      if (i .lt. IMAX) cycle
      if (i .eq. IMAX .and. JMAX .gt. 0) cycle
      if (NUM_PROCS .eq. 1) goto 60
!
!  ...skip printing partition
      goto 40
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
   40 continue
!
!  ...skip evaluating partition
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
      if(RANK .eq. ROOT) write(*,200) '6. calling MUMPS/PARDISO solver...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      if (NUM_PROCS .eq. 1) then
         !call pardiso_sc('H')
         call par_mumps_sc('H')
      else
         !call par_mumps_sc('H')
         call par_nested('H')
         !call par_hybrid('H')
      endif
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   70 continue
!
!  ...skip computing power
      !goto 80
!
      !if (i .lt. IMAX) cycle
      if (i .lt. IMAX+JMAX) cycle
!
!  ...compute power in fiber for signal field
      if(RANK .eq. ROOT) write(*,200) '7. computing power...'
      numPts = 2**min(i,IMAX); fld = 1
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      !call get_power(fld,numPts,-1) ! write to stdout
      !call get_power(fld,numPts,i-IMAX) ! write to file
      call get_power(fld,numPts,i-(IMAX+JMAX)) ! write to file
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   80 continue
!
!  ...skip writing paraview output
      !goto 90
!
      !if (i .lt. IMAX) cycle
      if (i .lt. IMAX+JMAX) cycle
!
!  ...write paraview output
      if(RANK .eq. ROOT) write(*,200) '8. writing paraview output...'
      iPvAttr = (/.false.,.false.,.false.,.false.,.true.,.false.,.false./)
      call paraview_select_attr(iPvAttr)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call my_paraview_driver
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   90 continue
   enddo
!
!..compute residual/error on last mesh
   if(RANK .eq. ROOT) write(*,200) 'Compute error/residual on last mesh...'
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   call refine_DPG(INOREFINEMENT,1,0.25d0,flag,physNick,ires, nstop)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()!
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
