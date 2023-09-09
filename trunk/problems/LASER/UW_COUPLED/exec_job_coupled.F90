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
   use MPI           , only: MPI_COMM_WORLD,MPI_Wtime
   use mpi_param     , only: RANK,ROOT,NUM_PROCS
   use par_mesh      , only: EXCHANGE_DOF,distr_mesh
   use paraview      , only: paraview_select_attr
   use zoltan_wrapper
!
   implicit none
!
   integer :: flag(6)
   logical :: iPvAttr(6)
   integer :: physNick,nstop,nskip
   logical :: ires
!
!..variables for nonlinear loop
   integer :: No,No1,No2
   real(8) :: L2NormDiff,stopEpsilon
   real(8) :: L2NormDiffIter(100)
   real(8) :: res
   real(8) :: FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!
   integer :: i,j,ierr,numPts,fld,time_step
   real(8) :: start_time,end_time
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
   NO_PROBLEM = 3
   call set_physAm(NO_PROBLEM, physNick,flag)
!
!..set .true. to compute the residual in all steps
   ires = .true.
!
!..set number of time steps for which to skip nonlinear Maxwell computation
!  (e.g., to wait with Maxwell until heat curve nears steady-state)
!  note: Maxwell solution is always computed in the initial time step
   nskip = 10
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
   call zoltan_w_set_lb(ZOLTAN_LB_DEFAULT)
!
!..refine the mesh anisotropically
   do i=1,IMAX+JMAX
!
      if(RANK .eq. ROOT) write(*,100) 'Mesh refinement step: i = ', i
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!     ...single anisotropic (in z) h-refinement
      if (i .le. IMAX) then
         if (RANK.eq.ROOT) write(*,200) '1. global anisotropic h-refinement...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call global_href_aniso(0,1) ! refine in z
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      else
!     ...adaptive h-refinement
      if (RANK.eq.ROOT) write(*,200) '1. adaptive h-refinement...'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call refine_DPG(ICLAD,1,0.25d0,flag,physNick,ires, nstop)
         !call refine_DPG(IADAPTIVE,1,0.25d0,flag,physNick,ires, nstop)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      endif
      if(RANK .eq. ROOT) write(*,300) end_time - start_time
!
      if (NUM_PROCS .eq. 1) cycle
!
!  ...set partitioner for load balancing, redistributes mesh in 'distr_mesh'
      if (i .eq. IMAX-2) then
         call zoltan_w_set_lb(ZOLTAN_LB_FIBER) ! fiber partitioner
      elseif (i .gt. IMAX-2) then
         goto 30 ! no load balancing
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
      if (i .le. IMAX .and. JMAX .gt. 0) cycle
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
!
!..set max number of nonlinear iterations per time step
   MAX_ITER = 10
!
!..start time stepping
   if (RANK.eq.ROOT) write(*,*) 'Begin time stepping..'
   do time_step = 0,NSTEPS ! initial step has ambient temperature (no solve)
      TIMESTEP = time_step ! setting global variable
      if (RANK.eq.ROOT) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
      if (RANK.eq.ROOT) write(*,4220) ' Proceeding with time step = ', time_step
      if (RANK.eq.ROOT) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
!
!  ...solve heat equation
      if (time_step .eq. 0) then
         if (RANK.eq.ROOT) write(*,4200) ' Skipping heat solve at initial time...'
      else
         if (RANK.eq.ROOT) write(*,4200) ' Solving the heat equation...'
      endif
!
      NO_PROBLEM = 2
      call set_physAm(NO_PROBLEM, physNick,flag)
      if (RANK.eq.ROOT) write(*,4200) ' Updating heat Dirichlet DOFs...'
      call update_Ddof
      if (time_step .gt. 0) then
         if (NUM_PROCS .eq. 1) then
            call pardiso_sc('H')
         else
            !call par_mumps_sc('H')
            call par_nested('H')
         endif
      endif
!  ...calculating temperature
      numPts = 2**IMAX
      call get_avgTemp(numPts,time_step)
      !call get_avgTemp(numPts,-1)
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!
!  ...Set default paraview output; all components = (/1,2,2,1,6,6/)
!     output temperature field and signal laser field
      iPvAttr = (/.true.,.false.,.false.,.false.,.true.,.false./)
!
!  ...Optionally, skip computing Maxwell problem for nskip time steps
!     (except in the initial time step)
      if (time_step.ge.1 .and. time_step.le.nskip) then
         if (RANK.eq.ROOT) write(*,4200) ' Skipping Maxwell solve in this time step...'
         iPvAttr = (/.true.,.false.,.false.,.false.,.false.,.false./) ! output only the heat solution
         goto 420
      endif
!
!  ...solve nonlinear Maxwell loop for SIGNAL and PUMP
      L2NormDiff = 1.d0
      FieldNormQ = 1.d0
!
!  ...update Dirichlet DOFs for signal field
      if (RANK.eq.ROOT) write(*,4200) ' Updating signal Dirichlet DOFs...'
      NO_PROBLEM = 3
      call set_physAm(NO_PROBLEM, physNick,flag)
      call update_Ddof
!
!  ...update Dirichlet DOFs for pump field
      if (FAKE_PUMP .ne. 1) then
         if (RANK.eq.ROOT) write(*,4200) ' Updating pump Dirichlet DOFs...'
         NO_PROBLEM = 4
         call set_physAm(NO_PROBLEM, physNick,flag)
         call update_Ddof
      endif
!
!  ...do until stopping criterion is satisfied
      if (RANK.eq.ROOT) write(*,4200) '---------------------------------------------'
      if (RANK.eq.ROOT) write(*,*)    ' Beginning nonlinear iterations...'
      if (RANK.eq.ROOT) write(*,4201) '---------------------------------------------'
 4200 format(/,A)
 4201 format(A,/)
!
      i = 0
!
!  ...do until stopping criterion is satisfied
      do
!
!     ...check stopping criterion
         if (i .eq. 0) goto 405
         L2NormDiffIter(i) = L2NormDiff/FieldNormQ
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
   405   continue
!
!     ...solve for signal first
         if (RANK.eq.ROOT) write(*,*) '   Signal solve...'
         NO_PROBLEM = 3
         call set_physAm(NO_PROBLEM, physNick,flag)
         if (NUM_PROCS .eq. 1) then
            call pardiso_sc('H')
         else
            call par_nested('H')
         endif
!
!     ...compute signal residual
         if (ires) then
            QUIET_MODE = .true.; IPRINT_TIME = 0
            if (RANK.eq.ROOT) write(*,4200) '   Signal residual:'
            call residual(res)
            QUIET_MODE = .false.; IPRINT_TIME = 1
         endif
!
!     ...if assuming a pump plane wave, skip the pump field computation
         if (FAKE_PUMP .eq. 1) then
            if (RANK.eq.ROOT) write(*,*) '   Assuming pump plane wave...'
            goto 410
         endif
!
!     ...next solve for pump
         if (RANK.eq.ROOT) write(*,*) '   Pump solve...'
         NO_PROBLEM = 4
         call set_physAm(NO_PROBLEM, physNick,flag)
         if (NUM_PROCS .eq. 1) then
            call pardiso_sc('H')
         else
            call par_nested('H')
         endif
!
!     ...compute pump residual
         if (ires) then
            QUIET_MODE = .true.; IPRINT_TIME = 0
            if (RANK.eq.ROOT) write(*,4200) '   Pump residual:'
            call residual(res)
            QUIET_MODE = .false.; IPRINT_TIME = 1
         endif
!
   410   continue
!
!     ...calculate norm corresponding to signal
         NO_PROBLEM = 3
         call set_physAm(NO_PROBLEM, physNick,flag)
         if (i .gt. 0 .or. time_step .gt. 0) then
            call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
            call get_Norm(flag,No2, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
         endif
         if (RANK.eq.ROOT) write(*,4240) '   L2NormDiff = ', L2NormDiff
         if (RANK.eq.ROOT) write(*,4240) '   FieldNormQ = ', FieldNormQ
  4240   format(A,F10.4)
!
!     ...copy current solution components of all fields into previous solution
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
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!
!  ...calculating power
      if (RANK.eq.ROOT) write(*,*) 'Computing power...'
      numPts = 2**IMAX; fld = 2
      if (FAKE_PUMP .eq. 1) fld = 1
      call get_power(fld,numPts,time_step)
      !call get_power(2,numPts,-1)
      if (RANK.eq.ROOT) write(*,*) ''
!
 420  continue
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!
!  ...only write paraview output every X time steps
      if (time_step.le.nskip .and. MOD(time_step,10).ne.0) goto 425
!
!  ...write paraview output
      if(RANK.eq.ROOT) write(*,200) ' Writing paraview output...'
      call paraview_select_attr(iPvAttr)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call my_paraview_driver
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if(RANK.eq.ROOT) write(*,300) end_time - start_time
 425  continue
!
!..end time-stepping loop
   enddo
!
!..compute final residual values (if not previously computed)
   if (.not. ires) then
      if (FAKE_PUMP .ne. 1) then
         QUIET_MODE = .true.; IPRINT_TIME = 0
         if (RANK.eq.ROOT) write(*,4200) '   Pump residual:'
         NO_PROBLEM = 4
         call set_physAm(NO_PROBLEM, physNick,flag)
         call residual(res)
      endif
!
      if (RANK.eq.ROOT) write(*,4200) '   Signal residual:'
      NO_PROBLEM = 3
      call set_physAm(NO_PROBLEM, physNick,flag)
      call residual(res)
      QUIET_MODE = .false.; IPRINT_TIME = 1
   endif
!
!..display stats
   if (RANK.eq.ROOT) then
      write(*,*) 'L2NormDiff/FieldNormQ:'
      do j=1,i
         write(*,4241) L2NormDiffIter(j)
 4241    format(es14.5)
      enddo
   endif
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
