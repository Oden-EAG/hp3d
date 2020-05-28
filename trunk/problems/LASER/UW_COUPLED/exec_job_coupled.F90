!----------------------------------------------------------------------
! exec_job_coupled: Coupled Maxwell/Heat fiber simulation
!----------------------------------------------------------------------
subroutine exec_job_coupled
!
   use assembly
   use assembly_sc
   use environment
   use commonParam
   use laserParam
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
!..variables for nonlinear loop
   integer :: No,No1,No2
   real(8) :: L2NormDiff,stopEpsilon
   real(8) :: L2NormDiffIter(100)
   real(8) :: FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!
   integer :: i,j,ierr,numPts,fld,time_step
   real(8) :: MPI_Wtime,start_time,end_time
!
!----------------------------------------------------------------------
!
   if(NONLINEAR_FLAG.ne.1) then
      write(*,*) 'NONLINEAR_FLAG must be 1 for coupled problem. stop.'
      stop
   endif
   if(HEAT_FLAG.ne.1) then
      write(*,*) 'HEAT_FLAG must be 1 for coupled problem. stop.'
      stop
   endif
   if(NEXACT.ne.0) then
      write(*,*) 'NEXACT must be 0 for coupled Maxwell/Heat problem. stop.'
      stop
   endif
!
   EXCHANGE_DOF = .false.
!
   if(RANK .eq. ROOT) then
      write(*,*) '=========================='
      write(*,*) 'exec_job_coupled: starting'
      write(*,*) '=========================='
   endif
!
!..distribute mesh initially
   call distr_mesh
!..set Zoltan partitioner
   call zoltan_w_set_lb(0)
!
!..refine the mesh anisotropically
   do i=1,IMAX
!
      if(RANK .eq. ROOT) write(*,100) 'Mesh refinement step: i = ', i
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!     ...single anisotropic (in z) h-refinement
      if(RANK .eq. ROOT) write(*,200) '1. global anisotropic h-refinement...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call global_href_aniso(0,1) ! refine in z
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
      if (NUM_PROCS .eq. 1) goto 30
!
      if (i .eq. IMAX-2) then
         call zoltan_w_set_lb(7)
      elseif (i .gt. IMAX-2) then
         goto 30
      endif
!
!  ...distribute mesh
      call update_gdof
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if(RANK .eq. ROOT) write(*,200) '2. distributing mesh...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call distr_mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
   30 continue
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
   enddo
!
!..update geometry degrees of freedom for current mesh
   call update_gdof
!
!..set stopping criterion for nonlinear iteration
   stopEpsilon = 1.0d-4
   L2NormDiff = 1.d0
   FieldNormQ = 1.d0
   L2NormDiffIter = 0.d0
!
!..set components (1: signal, 2: pump)
   No1 = 1; No2 = 2

!..start time stepping
   if (RANK.eq.ROOT) write(*,*) 'Begin time stepping..'
   do time_step = 1,NSTEPS
      TIMESTEP = time_step !setting global variable
      if (RANK.eq.ROOT) write(*,*) '---------------------------------------------'
      if (RANK.eq.ROOT) write(*,4220) ' Proceeding with time step = ', time_step
!  ...activate to compute maxwell only in first time step
      if(time_step .gt. 1 .and. time_step .le. 100) goto 420
!
!  ...first solve Nonlinear Maxwell loop for SIGNAL and PUMP
      physNick = 1
      L2NormDiff = 1.d0
      FieldNormQ = 1.d0
!              ...do until stopping criterion is satisfied
      if (RANK.eq.ROOT) write(*,4200) ' Beginning nonlinear iterations..'
 4200 format(/,A)
 4201 format(A,/)
!
      i = 0
!
!  ...do until stopping criterion is satisfied
      do
!
!     ...check stopping criterion
         L2NormDiffIter(i+1) = L2NormDiff/FieldNormQ
         if((L2NormDiff/FieldNormQ).lt.(stopEpsilon)) then
            if (RANK .eq. ROOT) then
               write(*,4230) '   ', L2NormDiff/FieldNormQ, ' < ', stopEpsilon
               write(*,*) ' Stopping criterion satisfied.'
               write(*,*) '---------------------------------------------'
               write(*,4210) ' Ending nonlinear loop after ', i, ' iterations.'
               write(*,*) '---------------------------------------------'
   4210        format(A,I3,A)
            endif
            exit
         elseif (RANK .eq. ROOT) then
            write(*,4220) ' Stopping criterion not yet satisfied. i = ', i
   4220     format(A,I3)
            write(*,4230) '   ', L2NormDiff/FieldNormQ, ' > ', stopEpsilon
   4230     format(A, F11.5,A,F11.5)
            write(*,*) '---------------------------------------------'
            write(*,4201) ' Proceed with nonlinear iterations..'
         endif
         if (i .ge. MAX_ITER) then
            if (RANK .eq. ROOT) write(*,4210) ' Ending nonlinear loop after', &
                                    MAX_ITER, ' iterations, no convergence.'
            exit
         endif
!
!     ...solve for signal first
         if (RANK.eq.ROOT) write(*,*) '   Signal solve..'
         NO_PROBLEM = 3
         call set_physAm(NO_PROBLEM, physNick,flag)
         call update_Ddof
         if (NUM_PROCS .eq. 1) then
            call pardiso_sc('H')
         else
            call par_nested('H')
         endif
!        QUIET_MODE = .true.; IPRINT_TIME = 0
!        if (RANK.eq.ROOT) write(*,*)
!        if (RANK.eq.ROOT) write(*,*) '   Signal residual:'
!        call residual
         if (RANK.eq.ROOT) write(*,*)
!        if (RANK.eq.ROOT) write(*,*)
!        QUIET_MODE = .false.; IPRINT_TIME = 1
!
!     ...next solve for pump
         if (RANK.eq.ROOT) write(*,*) '   Pump solve..'
         NO_PROBLEM = 4
         call set_physAm(NO_PROBLEM, physNick,flag)
         call update_Ddof
         if (NUM_PROCS .eq. 1) then
            call pardiso_sc('H')
         else
            call par_nested('H')
         endif
!        QUIET_MODE = .true.; IPRINT_TIME = 0
         if (RANK.eq.ROOT) write(*,*)
!        if (RANK.eq.ROOT) write(*,*) '   Pump residual:'
!        call residual
!        QUIET_MODE = .false.; IPRINT_TIME = 1
!
!     ...copy components and calculate norm corresponding to signal
         NO_PROBLEM = 3
         call set_physAm(NO_PROBLEM, physNick,flag)
         if (i .gt. 0) then
            call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
            call get_Norm(flag,No2, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
         endif
         if (RANK.eq.ROOT) write(*,4240) '   L2NormDiff = ', L2NormDiff
         if (RANK.eq.ROOT) write(*,4240) '   FieldNormQ = ', FieldNormQ
  4240   format(A,F10.4)
         if (RANK.eq.ROOT) write(*,*)
         if (RANK.eq.ROOT) write(*,*) 'copy_coms...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call copy_coms(No1,No2)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
         if (RANK.eq.ROOT) write(*,300) end_time - start_time
         if (RANK.eq.ROOT) write(*,*)
         i = i + 1
!  ...end nonlinear loop
      enddo
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!  ...calculating power
      !if (time_step .eq. 1) then
         if (RANK.eq.ROOT) write(*,*) 'Computing power..'
         numPts=32
         !call get_power(2,numPts,time_step-1)
         call get_power(2,numPts,-1)
         write(*,*) ''
      !endif
 420  continue
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!  ...now solve heat equation
      if (RANK.eq.ROOT) write(*,*) ''
      if (RANK.eq.ROOT) write(*,4210) ' Computing time-step ', time_step, &
                                      ' for heat equation..'
      if (RANK.eq.ROOT) write(*,*) ''
      NO_PROBLEM = 2
      call set_physAm(NO_PROBLEM, physNick,flag)
      call update_Ddof
      if (NUM_PROCS .eq. 1) then
         call pardiso_sc('H')
      else
         call par_nested('H')
      endif
!  ...calculating temperature
      !call get_avgTemp(numPts,time_step-1)
      call get_avgTemp(numPts,-1)
!
!..end time-stepping loop
   enddo
!
!!..Last step only display (no refinement)
!   QUIET_MODE = .true.; IPRINT_TIME = 0
!   if (RANK.eq.ROOT) write(*,*)
!   if (RANK.eq.ROOT) write(*,*) '   Pump residual:'
!   NO_PROBLEM = 4
!   call set_physAm(NO_PROBLEM, physNick,flag)
!   call residual
!!
!   if (RANK.eq.ROOT) write(*,*)
!   if (RANK.eq.ROOT) write(*,*) '   Signal residual:'
!   NO_PROBLEM = 3
!   call set_physAm(NO_PROBLEM, physNick,flag)
!   call update_Ddof
!   if (NUM_PROCS .eq. 1) then
!      call pardiso_sc('H')
!   else
!      call par_nested('H')
!   endif
!   call residual
!   QUIET_MODE = .false.; IPRINT_TIME = 1
!!
!   if (RANK.eq.ROOT) then
!      write(*,*) 'L2NormDiff/FieldNormQ:'
!      do j=1,i+1
!         write(*,4241) L2NormDiffIter(j)
! 4241    format(es14.5)
!      enddo
!   endif
!!
!!..compute power in fiber for signal and pump field
!   if(RANK .eq. ROOT) write(*,200) '7. computing power...'
!   numPts = 2**IMAX; fld = 2
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!   call get_power(fld,numPts,-1)
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
!   if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
!
  100 format(/,'/////////////////////////////////////////////////////////////', &
             /,'             ',A,I2,/)
  200 format(/,'-------------------------------------------------------------', &
             /,'  --> ',A,/)
  300 format(' timer: ',f12.5,' seconds')
!
   if(RANK .eq. ROOT) then
      write(*,*)
      write(*,*) '=========================='
      write(*,*) 'exec_job_coupled: finished'
      write(*,*) '=========================='
      write(*,*)
   endif
!
end subroutine exec_job_coupled
