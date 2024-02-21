!----------------------------------------------------------------------
! exec_job_nl: Nonlinear gain fiber simulation
!----------------------------------------------------------------------
subroutine exec_job_nl
!
   use assembly
   use assembly_sc
   use environment
   use commonParam
   use laserParam
   use control
   use data_structure3D
   use mpi_wrapper
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
!..variables for nonlinear loop
   integer :: currAttr,prevAttr
   real(8) :: L2NormDiff,stopEpsilon
   real(8) :: L2NormDiffIter(100),SignalRes(100),PumpRes(100)
   real(8) :: fieldNormQ
!
   integer :: i,j,ierr,numPts,fld
   real(8) :: start_time,end_time,aux_time
!
!----------------------------------------------------------------------
!
   if(NONLINEAR_FLAG.ne.1) then
      write(*,*) 'NONLINEAR_FLAG must be 1 for nonlinear problem. stop.'
      stop
   endif
   if(NEXACT.ne.0) then
      write(*,*) 'NEXACT must be 0 for nonlinear Maxwell problem. stop.'
      stop
   endif
!
!..Initiate pump field ODE solution
   if(PLANE_PUMP.eq.2) then
      numPts = 2**IMAX
      call pump_ode_alloc(numPts)
   endif
!
   EXCHANGE_DOF = .false.
!
   NO_PROBLEM = 3
   call set_physAm(NO_PROBLEM, physNick,flag)
!
!..set .true. to compute the residual in all steps
   ires = .false.
!
   if (RANK.eq.ROOT) then
      write(*,*) '====================='
      write(*,*) 'exec_job_nl: starting'
      write(*,*) '====================='
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
      if (RANK.eq.ROOT) write(*,100) 'Mesh refinement step: i = ', i
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if (i .le. IMAX) then
!     ...single anisotropic (in z) h-refinement
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
      if (RANK.eq.ROOT) write(*,300) end_time - start_time
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
      if (RANK.eq.ROOT) write(*,200) '2. distributing mesh...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call distr_mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if (RANK.eq.ROOT) write(*,300) end_time - start_time
!
   30 continue
      if (i .le. IMAX .and. JMAX .gt. 0) cycle
!
!  ...skip printing partition
      goto 40
!
!  ...print current partition (elems)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if (RANK.eq.ROOT) write(*,200) '3. printing current partition (elems)...'
      if (i .le. 2) then
         call print_partition
      else
         if (RANK.eq.ROOT) write(*,*) '   ... skipping for a large number of elements.'
      endif
!
   40 continue
!
!  ...skip evaluating partition
      goto 50
!
!  ...evaluate current partition
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if (RANK.eq.ROOT) write(*,200) '4. evaluating current partition...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call zoltan_w_eval
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if (RANK.eq.ROOT) write(*,300) end_time - start_time
!
   50 continue
!
!  ...run mesh verification routines
      call MPI_BARRIER (MPI_COMM_WORLD, ierr);
      if (RANK.eq.ROOT) write(*,200) '5. verify distributed mesh consistency...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call par_verify
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if (RANK.eq.ROOT) write(*,300) end_time - start_time
!
   60 continue
!
   enddo
!
!..update geometry degrees of freedom for current mesh
   call update_gdof
!
!..update Dirichlet DOFs for signal field
   if (RANK.eq.ROOT) write(*,4200) ' Updating signal Dirichlet DOFs...'
   NO_PROBLEM = 3
   call set_physAm(NO_PROBLEM, physNick,flag)
   call update_Ddof
!
!..update Dirichlet DOFs for pump field
   if (PLANE_PUMP.eq.0) then
      if (RANK.eq.ROOT) write(*,4200) ' Updating pump Dirichlet DOFs...'
      NO_PROBLEM = 4
      call set_physAm(NO_PROBLEM, physNick,flag)
      call update_Ddof
   endif
!
!..set stopping criterion for nonlinear iteration
   stopEpsilon = 1.0d-4
   L2NormDiff = 1.d0
   fieldNormQ = 1.d0
   L2NormDiffIter = 0.d0
!
   SignalRes = 0.d0
   PumpRes = 0.d0
!
!..set attributes for norm of signal difference
   currAttr = 5 ! signal field, current solution
   prevAttr = 7 ! signal field, previous iterate (auxiliary solution)
!
   if (RANK.eq.ROOT) write(*,200) '6. Beginning nonlinear iterations...'
 4200 format(/,A)
 4201 format(A,/)
!
   i = 0
!
!..do until stopping criterion is satisfied
   do
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); aux_time = MPI_Wtime()
!
!  ...check stopping criterion
      if (i .eq. 0) goto 405
      L2NormDiffIter(i) = L2NormDiff/fieldNormQ
      if((L2NormDiff/fieldNormQ).lt.stopEpsilon) then
         if (RANK.eq.ROOT) then
            write(*,4230) '   ', L2NormDiff/fieldNormQ, ' < ', stopEpsilon
            write(*,*) ' Stopping criterion satisfied.'
            write(*,*) '---------------------------------------------'
            write(*,4210) ' Ending nonlinear loop after ', i, ' iterations.'
            write(*,*) '---------------------------------------------'
     4210   format(A,I3,A)
         endif
         exit
      elseif (RANK.eq.ROOT) then
         write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,4220) ' Stopping criterion not yet satisfied. i = ', i
   4220  format(A,I3)
         write(*,4230) '   ', L2NormDiff/fieldNormQ, ' > ', stopEpsilon
   4230  format(A, F11.5,A,F11.5)
         write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,4201) ' Proceed with nonlinear iterations..'
      endif
!
      if (i.ge.MAX_ITER) then
         if (RANK.eq.ROOT) write(*,4210) ' Ending nonlinear loop after', &
                               MAX_ITER, ' iterations, no convergence.'
         exit
      endif
!
  405 continue
!
!  ...solve for signal first
      if (RANK.eq.ROOT) write(*,4200) '   Signal solve..'
      NO_PROBLEM = 3
      call set_physAm(NO_PROBLEM, physNick,flag)
#if HP3D_USE_INTEL_MKL
      if (NUM_PROCS .eq. 1) then
         call pardiso_sc('H')
      else
         call par_nested('H')
      endif
#else
      call par_mumps_sc('H')
#endif
!
!  ...compute signal residual
      if (ires) then
         QUIET_MODE = .true.; IPRINT_TIME = 0
         if (RANK.eq.ROOT) write(*,4200) '   Signal residual:'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call residual(SignalRes(i+1))
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
         if (RANK.eq.ROOT) write(*,300) end_time - start_time
         QUIET_MODE = .false.; IPRINT_TIME = 1
      endif
!
!  ...if assuming a pump plane wave, skip the pump field computation
      if (PLANE_PUMP.eq.1) then
         if (RANK.eq.ROOT) write(*,*) '   Assuming constant pump plane wave...'
         goto 410
      elseif (PLANE_PUMP.eq.2) then
         if (RANK.eq.ROOT) write(*,*) '   Computing pump plane wave solution by ODE model...'
         call pump_ode_solve
         goto 410
      endif
!
!  ...next solve for pump
      if (RANK.eq.ROOT) write(*,4200) '   Pump solve..'
      NO_PROBLEM = 4
      call set_physAm(NO_PROBLEM, physNick,flag)
#if HP3D_USE_INTEL_MKL
      if (NUM_PROCS .eq. 1) then
         call pardiso_sc('H')
      else
         call par_nested('H')
      endif
#else
      call par_mumps_sc('H')
#endif
!
!  ...compute pump residual
      if (ires) then
         QUIET_MODE = .true.; IPRINT_TIME = 0
         if (RANK.eq.ROOT) write(*,4200) '   Pump residual:'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call residual(PumpRes(i+1))
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
         if (RANK.eq.ROOT) write(*,300) end_time - start_time
         QUIET_MODE = .false.; IPRINT_TIME = 1
      endif
!
  410 continue
!
!  ...calculate norm corresponding to signal
      NO_PROBLEM = 3
      call set_physAm(NO_PROBLEM, physNick,flag)
      if (i .gt. 0) then
         call get_L2NormDiff(currAttr,prevAttr, L2NormDiff)
         call get_L2NormAttr(prevAttr, fieldNormQ)
      endif
      if (RANK.eq.ROOT) write(*,4240) '   L2NormDiff = ', L2NormDiff
      if (RANK.eq.ROOT) write(*,4240) '   FieldNormQ = ', fieldNormQ
 4240 format(A,F10.4)
!
!  ...copy current solution components of all fields into previous solution
      if (RANK.eq.ROOT) write(*,4200) 'copy_attr...'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call copy_attr(currAttr,prevAttr)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
      if (RANK.eq.ROOT) write(*,300) end_time - start_time
      if (RANK.eq.ROOT) write(*,*)
!
      if (RANK.eq.ROOT) write(*,4210) ' Iteration ',i,' total time:'
      if (RANK.eq.ROOT) write(*,300) end_time - aux_time
      i = i + 1
!  ...end nonlinear loop
   end do
!
!..compute final residual values (if not previously computed)
   if (.not. ires) then
      if (PLANE_PUMP.eq.0) then
         QUIET_MODE = .true.; IPRINT_TIME = 0
         if (RANK.eq.ROOT) write(*,4200) '   Pump residual:'
         NO_PROBLEM = 4
         call set_physAm(NO_PROBLEM, physNick,flag)
         call residual(PumpRes(i))
      endif
!
      if (RANK.eq.ROOT) write(*,4200) '   Signal residual:'
      NO_PROBLEM = 3
      call set_physAm(NO_PROBLEM, physNick,flag)
      call residual(SignalRes(i))
      QUIET_MODE = .false.; IPRINT_TIME = 1
   endif
!
!..display stats
   if (RANK.eq.ROOT) then
      write(*,*) 'L2NormDiff/FieldNormQ:'
      do j=1,i
         write(*,4241) L2NormDiffIter(j)
      enddo
      if (ires) then
         write(*,*) 'Signal Residuals:'
         do j=1,i
            write(*,4241) SignalRes(j)
         enddo
         if (PLANE_PUMP.eq.0) then
            write(*,*) 'Pump Residuals:'
            do j=1,i
               write(*,4241) PumpRes(j)
            enddo
         endif
      endif
 4241 format(es14.5)
   endif
!
!..compute power in fiber for signal and pump field
   if (RANK.eq.ROOT) write(*,200) '7. computing power...'
   numPts = 2**IMAX; fld = 2
   if (PLANE_PUMP.eq.1) fld = 1
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   !call get_power(fld,numPts,-1) ! write to stdout
   call get_power(fld,numPts,0) ! write to file
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
   if (RANK.eq.ROOT) write(*,300) end_time - start_time
!
!..write paraview output
   if (RANK.eq.ROOT) write(*,200) '8. writing paraview output...'
   iPvAttr = (/.false.,.false.,.false.,.false.,.true.,.false.,.false./)
   call paraview_select_attr(iPvAttr)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   call my_paraview_driver
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
   if (RANK .eq. ROOT) write(*,300) end_time - start_time
!
   if (PLANE_PUMP.eq.2) call pump_ode_dealloc
!
  100 format(/,'/////////////////////////////////////////////////////////////', &
             /,'             ',A,I2,/)
  200 format(/,'-------------------------------------------------------------', &
             /,'  --> ',A,/)
  300 format(' timer: ',f12.5,' seconds')
!
   if(RANK.eq.ROOT) then
      write(*,*)
      write(*,*) '====================='
      write(*,*) 'exec_job_nl: finished'
      write(*,*) '====================='
      write(*,*)
   endif
!
end subroutine exec_job_nl
