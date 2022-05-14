!----------------------------------------------------------------------
!
!     program name      - main
!
!----------------------------------------------------------------------
!
!     latest revision:  - April 2022
!
!     purpose:          - main driver for MPI test program
!                         sheathed hose
!
!----------------------------------------------------------------------
!
program main
!
  use environment
  use common_prob_data
  use data_structure3D
  use GMP
  use control
  use assembly
  use assembly_sc, only: IPRINT_TIME
  use stc        , only: STORE_STC,HERM_STC
!
  use MPI        , only: MPI_COMM_WORLD
  use mpi_param  , only: ROOT,RANK,NUM_PROCS
  use mpi_wrapper, only: mpi_w_init,mpi_w_finalize
!
  implicit none
!
!..auxiliary variables
   integer :: i, ierr
!
!..OMP variables
   integer :: num_threads, omp_get_num_threads
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
!
!-----------------------------------------------------------------------
!
!..Initialize MPI environment
   call mpi_w_init
!
!..Set common hp3D environment parameters (reads in options arguments)
   call begin_environment
   call set_environment
   call end_environment
!
   if (RANK .eq. ROOT) then
!  ...print header
      write(6,*)'                                                        '
      write(6,*)'//                                                    //'
      write(6,*)'// --            HYBRID DPG FORMULATION            -- //'
      write(6,*)'//                                                    //'
      write(6,*)'// -- SHEATHED HOSE PROBLEM FOR LINEAR ELASTICITY  -- //'
      write(6,*)'//                                                    //'
      write(6,*)'                                                        '
   endif
   flush(6)
!
   IPRINT_TIME = 1
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..initialize physics, geometry, etc.
   do i = 0, NUM_PROCS-1
      if ((RANK .eq. i) .and. (RANK .eq. ROOT)) then
         write(6,*)
         write(6,1020) "Master proc [", RANK, "], initialize.."
         QUIET_MODE = .FALSE.
      else if ((RANK .eq. i) .and. (RANK .ne. ROOT)) then
         write(6,1020) "Worker proc [", RANK, "], initialize.."
         QUIET_MODE = .TRUE.
      endif
   enddo
 1020 format (A,I3,A)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   call initialize
!
!   IBC_PROB : 0 - uniform traction ; 1 - clamped ends ; 2 - free ends ; 3 - periodic ends
    if ((IBC_PROB.eq.2).or.(IBC_PROB.eq.3)) then
      call remove_RBM(IBC_PROB)
      call update_Ddof
    endif
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
   if (RANK .eq. ROOT) write(6,1015) end_time-start_time
 1015 format(' initialize : ',f12.5,' seconds',/)
!
!..FLAGS
 ISTC_FLAG = .true.
 STORE_STC = .true.
 HERM_STC  = .false.
!
!..determine number of omp threads running
 if (RANK .eq. ROOT) then
    write(6,1025) ' Initial polynomial order: ',IP
!$OMP parallel
!$OMP single
    num_threads = omp_get_num_threads()
    write(6,1025) ' Number of OpenMP threads: ',num_threads
1025 format(A,I2)
!$OMP end single
!$OMP end parallel
 endif
!
 if (JOB .ne. 0) then
    call exec_job
 else
    if (RANK .eq. ROOT) then
       call master_main
    else
       call worker_main
    endif
 endif
!
 call finalize
 call mpi_w_finalize
!..END MPI
!
end program main
!
!
!----------------------------------------------------------------------
! master_main
!----------------------------------------------------------------------
subroutine master_main()
  !
     use environment
     use common_prob_data
     use data_structure3D
     use GMP
  !
     use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
     use mpi_param     , only: ROOT,RANK,NUM_PROCS
     use par_mesh      , only: DISTRIBUTED,HOST_MESH
     use zoltan_wrapper, only: zoltan_w_set_lb
  !
     implicit none
  !
  !..MPI variables
     integer :: ierr
  !
  !..auxiliary variables
     integer :: idec, r, lb, count, src
  !
  !----------------------------------------------------------------------
  !
     if (RANK .ne. ROOT) then
        write(*,*) 'master_main: RANK .ne. ROOT'
        stop
     endif
  !
  !..start user interface, with idec
  !..broadcast user command to workers
  !
  !..test accessing data structures
     write(6,8020) '[', RANK, '] : ', 'NRELIS,NRELES,NRNODS = ',NRELIS,NRELES,NRNODS
   8020 format(A,I3,A,A,I4,', ',I4,', ',I4)
  !
     flush(6)
     call MPI_BARRIER (MPI_COMM_WORLD, ierr)
  !
  #if DEBUG_MODE
     write(*,*) '========================='
     write(*,*) '  RUNNING in DEBUG_MODE  '
     write(*,*) '========================='
  #endif
  !
  !..display menu in infinite loop
     idec = 1
     do while(idec /= 0)
  !
        write(*,*)
        write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
        write(*,*) 'SELECT'
        write(*,*) 'QUIT....................................0'
        write(*,*) '                                         '
        write(*,*) '     ---- Visualization, I/O ----        '
        write(*,*) 'HP3D graphics (graphb)..................1'
        write(*,*) 'HP3D graphics (graphg)..................2'
        write(*,*) 'Paraview................................3'
        write(*,*) '                                         '
        write(*,*) '    ---- Print Data Structure ----       '
        write(*,*) 'Print arrays (interactive).............10'
        write(*,*) 'Print data structure arrays............11'
        write(*,*) 'Print current partition (elems)........15'
        write(*,*) 'Print current subdomains (nodes).......16'
        write(*,*) 'Print partition coordinates............17'
        write(*,*) '                                         '
        write(*,*) '        ---- Refinements ----            '
        write(*,*) 'Single uniform h-refinement............20'
        write(*,*) 'Single uniform p-refinement............21'
        write(*,*) 'Multiple uniform h-refs + solve........22'
        write(*,*) 'Single anisotropic h-refinement (z)....23'
        write(*,*) 'Refine a single element................26'
        write(*,*) '                                         '
        write(*,*) '        ---- MPI Routines ----           '
        write(*,*) 'Distribute mesh........................30'
        write(*,*) 'Collect dofs on ROOT...................31'
        write(*,*) 'Suggest mesh partition (Zoltan)........32'
        write(*,*) 'Evaluate mesh partition (Zoltan).......33'
        write(*,*) 'Run verification routines..............35'
        write(*,*) '                                         '
        write(*,*) '          ---- Solvers ----              '
        write(*,*) 'MUMPS (MPI)............................40'
        write(*,*) 'MUMPS (OpenMP).........................41'
        write(*,*) 'Pardiso (OpenMP).......................42'
        write(*,*) 'Frontal (Seq)..........................43'
        write(*,*) 'MUMPS (Nested Dissection)..............44'
        write(*,*) 'PETSc (MPI)............................45'
        write(*,*) '                                         '
        write(*,*) '     ---- Error and Residual ----        '
        write(*,*) 'Compute exact error....................50'
        write(*,*) '                                         '
        write(*,*) '          ---- TESTING ----              '
        write(*,*) 'Flush dof, update_gdof, update_Ddof....60'
        write(*,*) 'P-refine an element....................65'
        write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
  !
        read( *,*) idec
        write(6,8010) '[', RANK, '] : ','Broadcast: idec = ', idec
   8010 format(A,I3,A,A,I3)
        count = 1; src = ROOT
        call MPI_BCAST (idec,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
  !
        select case(idec)
  !     ...QUIT
           case(0) ; goto 89
  !
  !     ...HP3D graphics
           case(1)
              if (DISTRIBUTED .and. (.not. HOST_MESH)) then
                 write(*,*) 'Cannot use HP3D graphics with distributed mesh!'
              else
                 call graphb
              endif
           case(2)
              if (DISTRIBUTED .and. (.not. HOST_MESH)) then
                 write(*,*) 'Cannot use HP3D graphics with distributed mesh!'
              else
                 call graphg
              endif
  !
  !     ...Paraview graphics
           case(3) ; call exec_case(idec)
  !
  !     ...Print data structure
           case(10,11)
              r = ROOT
              if (NUM_PROCS .gt. 1) then
                 write(*,*) 'Select processor RANK: '
                 read (*,*) r
                 count = 1; src = ROOT
                 call MPI_BCAST (r,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
              endif
              if (r .eq. RANK) then
                 call exec_case(idec)
              endif
           case(15,16,17)
              call exec_case(idec)
  !
  !     ...Refinements
           case(20,21,22,23,26)
              call exec_case(idec)
  !
  !     ...MPI Routines
           case(31,33,35)
              call exec_case(idec)
  !
  !     ...MPI Routines (partitioner)
           case(30,32)
  !        ...Load balancing strategy
              if (DISTRIBUTED) then
                 write(*,*) 'Select load balancing strategy:'
                 write(*,*) '  0: nelcon'
                 write(*,*) '  1: BLOCK'
                 write(*,*) '  2: RANDOM'
                 write(*,*) '  3: RCB'
                 write(*,*) '  4: RIB'
                 write(*,*) '  5: HSFC'
                 write(*,*) '  6: GRAPH'
                 write(*,*) '  7: FIBER'
                 read (*,*) lb
                 count = 1; src = ROOT
                 call MPI_BCAST (lb,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
                 call zoltan_w_set_lb(lb)
              endif
              call exec_case(idec)
  !
  !     ...Solvers
           case(40,41,42,43,44,45)
              call exec_case(idec)
  !
  !     ...Error and Residual
           case(50)
              call exec_case(idec)
  !
  !     ...TODO testing
           case(60)
              call exec_case(idec)
  !
           case(65)
              call exec_case(idec)
  !
        end select
  !
        call MPI_BARRIER (MPI_COMM_WORLD, ierr)
  !
  !..end infinite loop
     enddo
  !
     call MPI_BARRIER (MPI_COMM_WORLD, ierr)
  !
  89 continue
     write(6,8030) '[', RANK, '] : ','master_main end.'
   8030 format(A,I3,A,A)
  !
  end subroutine master_main
  !
  !
!   display menu
  solved = .FALSE.
  idec=1
!   display menu in infinite loop
  do while(idec /= 0)
  !
    write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    write(*,*) 'SELECT'
    write(*,*) 'QUIT ...................................0'
    write(*,*) '                                         '
    write(*,*) 'Geometry graphics (X11) ................1'
    write(*,*) 'HP3D graphics (X11) ....................2'
    write(*,*) 'Paraview ...............................3'
    ! write(*,*) 'Compute BC Interpolation Error .........5'
    write(*,*) '                                         '
    write(*,*) 'Print Data Structure arrays ...........10'
    write(*,*) 'Dumpout Data Structure ................11'
    write(*,*) '                                         '
    write(*,*) ' -- Refinements --                       '
    write(*,*) 'Single Uniform H-refinement ...........20'
    ! write(*,*) 'Single Adaptive H-Refinements .........21'
    write(*,*) 'Single Custom H-Refinement ............22'
    write(*,*) 'Multi-step Uniform H-refinement .......23'
    ! write(*,*) 'Multi-step Adaptive H-refinement ......24'
    write(*,*) 'Multi-step Custom H-refinement ........25'
    write(*,*) '                                         '
    write(*,*) ' -- Solves --                            '
    write(*,*) 'FRONTAL SOLVE PROBLEM .................30'
    write(*,*) 'MUMPS SOLVE (SEQ.) ....................40'
    write(*,*) 'MUMPS SOLVE (OMP) .....................45'
    write(*,*) 'UHM SOLVE (PAR.) ......................50'
    write(*,*) '                                         '
    write(*,*) 'EXACT ERROR ...........................60'
    ! write(*,*) 'Write error to file ...................61'
    write(*,*) 'RESIDUAL ..............................70'
    write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    read( *,*) idec
  !
    select case(idec)
!  QUIT
    case( 0) ; call finalize ; stop
!-----------------------------------------------------------------------
!                                 GRAPHICS
!  GMP x11 graphics
    case(1) ; call graphg
!
!  hp3d x11 graphics
    case(2) ; call graphb
!
!  Paraview graphics
    case(3) ; call paraview_driver
              ! call paraview_stress_driver
!
!  BC Interpolation Error
    ! case(5) ; call compute_BC_interp_error
!-----------------------------------------------------------------------
!                                 ?
!  print data structure
    case(10) ; call result
!
!  dump out
    case(11) ; call dumpout
!-----------------------------------------------------------------------
!                               REFINEMENTS
!-----------------------------------------------------------------------
!
!   Single uniform refinement
    case(20)
      if (.not.solved) write(*,*) 'You have not solved.'
      if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
      call refine_DPG(IUNIFORM,1,0.25d0, nstop)
      solved=.FALSE.
! !   Single adaptive refinement
!     case(21)
!       if (.not.solved) then
!         write(*,*) 'You have not solved. Cannot adaptively refine.'
!       else
!         nreflag=0
!         do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
!           write(*,*) 'REFINEMENT FLAG:h-refine=1,p-refine=2,hp-refine=3'
!           write(*,*) '0.d0<FACTOR<1.d0'
!           write(*,*) 'Provide: REFINEMENT FLAG, FACTOR'
!           read(*,*) nreflag,factor
!         enddo
!         call refine_DPG(IADAPTIVE,nreflag,factor, nstop)
!         if (nstop.eq.1) write(*,*) 'No elements were refined.'
!         if (nstop.eq.0) solved=.FALSE.
!       endif
!   Single custom refinement
    case(22)
      if (.not.solved) write(*,*) 'You have not solved.'
      if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
      call refine_DPG(ICUSTOMREF,1,0.25d0, nstop)
      solved=.FALSE.
!   Multi-step uniform h refinement
    case(23)
      nsteps=0
      do while (nsteps.le.0)
        write(*,*) 'NUMBER OF REFINEMENTS>0'
        write(*,*) 'SOLVER:frontal=1,MUMPS=2,UHM=3,OMP MUMPS=4'
        write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
        read(*,*) nsteps,idec_solve
      enddo
      do i=0,nsteps
!    ...solve first if needed
        if (.not.solved) then
          select case(idec_solve)
          case(1)
            call solve1(NR_RHS_PROB)
          case(2)
            call mumps_solve_seq(NR_RHS_PROB)
          case(3)
            call uhm_solve
            call uhm_solver_flush(UHM_SOLVER_PTR)
          case(4)
            call mumps_interf(NR_RHS_PROB)
          end select
        endif
!    ...declare solved and save results to paraview file
        solved=.TRUE.
        ! call paraview_driver
        ! PREFIX_tmp=PREFIX
        ! PARAVIEW_DOMAIN=1
        ! PREFIX="Steel_"
        ! call paraview_custom_driver
        ! PARAVIEW_DOMAIN=2
        ! PREFIX="Rubber_"
        ! call paraview_custom_driver
        ! PREFIX=PREFIX_tmp
        call paraview_custom_dump
!    ...then display error and refine if necessary
        if (i.ne.nsteps) then
          call refine_DPG(IUNIFORM,1,0.25d0, nstop)
          if (nstop.eq.1) then
            write(*,*) 'No elements were refined.'
            write(*,7000) i
 7000       format('Exiting loop after ',i2,' refinements...')
            cycle
          else
            solved=.FALSE.
          endif
        else ! Last step only display (no refinement)
          call refine_DPG(INOREFINEMENT,1,0.25d0, nstop)
        endif
      enddo
! !   Multi-step adaptive h refinement
!     case(24)
!       nreflag=0
!       do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
!         write(*,*) 'NUMBER OF REFINEMENTS>0'
!         write(*,*) 'REFINEMENT FLAG:h-refine=1,p-refine=2,hp-refine=3'
!         write(*,*) '0.d0<FACTOR<1.d0'
!         write(*,*) 'SOLVER:frontal=1,MUMPS=2,UHM=3,OMP MUMPS=4'
!         write(*,*) 'Provide: NUMBER OF REFINEMENTS,',  &
!                    ' REFINEMENT FLAG, FACTOR, SOLVER'
!         read(*,*) nsteps,nreflag,factor,idec_solve
!       enddo
!       do i=0,nsteps
! !    ...solve first if needed
!         if (.not.solved) then
!           select case(idec_solve)
!           case(1)
!             call solve1(NR_RHS_PROB)
!           case(2)
!             call mumps_solve_seq(NR_RHS_PROB)
!           case(3)
!             call uhm_solve
!             call uhm_solver_flush(UHM_SOLVER_PTR)
!           case(4)
!             call mumps_interf(NR_RHS_PROB)
!             write(*,8000) MTime(1), MTime(2), MTime(3), MTime(4)
!  8000 format('-- SOLVE --',/,  &
!              'Step 1: determine the first dof offsets for active nodes',/,  &
!              'Time: ',f15.8,/,  &
!              'Step 2: compute element matrices, assemble global stiffness matrix',/,  &
!              'Time: ',f15.8,/,  &
!              'Step 3: MUMPS solve',/,  &
!              'Time: ',f15.8,/,  &
!              'Step 4: reconstruct global to local connectivities',/,  &
!              'Time: ',f15.8)
!           end select
!         endif
! !    ...say it has solved and save results to paraview file
!         solved=.TRUE.
!         call paraview_driver
!         call paraview_stress_driver
! !    ...then display error and refine if necessary
!         if (i.ne.nsteps) then
!           call refine_DPG(IADAPTIVE,nreflag,factor, nstop)
!           if (nstop.eq.1) then
!             write(*,*) 'No elements were refined.'
!             write(*,7000) i
!             cycle
!           else
!             solved=.FALSE.
!           endif
!         else ! Last step only display (no refinement)
!             call refine_DPG(INOREFINEMENT,nreflag,factor, nstop)
!         endif
!         write(*,8001) MTime(5)
!  8001 format('-- REFINE --',/,  &
!              'Calculate residual :',/,  &
!              'Time: ',f15.8)
!       enddo
!   Multi-step uniform h refinement
    case(25)
      nsteps=0
      do while (nsteps.le.0)
        write(*,*) 'NUMBER OF REFINEMENTS>0'
        write(*,*) 'SOLVER:frontal=1,MUMPS=2,UHM=3,OMP MUMPS=4'
        write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
        read(*,*) nsteps,idec_solve
      enddo
      do i=0,nsteps
!    ...solve first if needed
        if (.not.solved) then
          select case(idec_solve)
          case(1)
            call solve1(NR_RHS_PROB)
          case(2)
            call mumps_solve_seq(NR_RHS_PROB)
          case(3)
            call uhm_solve
            call uhm_solver_flush(UHM_SOLVER_PTR)
          case(4)
            call mumps_interf(NR_RHS_PROB)
          end select
        endif
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        ! call paraview_driver
        ! call paraview_stress_driver
!    ...then display error and refine if necessary
        if (i.ne.nsteps) then
          call refine_DPG(ICUSTOMREF,1,0.25d0, nstop)
          if (nstop.eq.1) then
            write(*,*) 'No elements were refined.'
            write(*,9000) i
 9000       format('Exiting loop after ',i2,' refinements...')
            cycle
          else
            solved=.FALSE.
          endif
        else ! Last step only display (no refinement)
          call refine_DPG(INOREFINEMENT,1,0.25d0, nstop)
        endif
      enddo
!-----------------------------------------------------------------------
! !  uniform global h-refinements
!     case(20)
!       call global_href
!       call update_gdof
!       call update_ddof
!       call verify_orient
!       call verify_neig
!       solved=.FALSE.
!-----------------------------------------------------------------------
!   global p-enrichments
    ! case(24) ; call global_pref
!-----------------------------------------------------------------------
!                                 SOLVERS
!  Frontal Solver
    case(30)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call solve1(NR_RHS_PROB)
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        call paraview_driver
        ! call paraview_stress_driver
      endif
!
!  MUMPS
    case(40)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call uhm_time_in
        call mumps_solve_seq(NR_RHS_PROB)
        call uhm_time_out(t)
        write(*,*) 'time mumps = ', t
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        call paraview_driver
        ! call paraview_stress_driver
      endif
!
!  OMP MUMPS
    case(45)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call uhm_time_in
        call mumps_interf(NR_RHS_PROB)
        call uhm_time_out(t)
        write(*,*) 'time mumps = ', t
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        call paraview_driver
        ! call paraview_stress_driver
      endif
!
!  UHM
    case(50)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call uhm_time_in
        call uhm_solve
        call uhm_time_out(t)
        write(*,*) 'time uhm = ', t
        call uhm_solver_flush(UHM_SOLVER_PTR)
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        call paraview_driver
        ! call paraview_stress_driver
      endif
!-----------------------------------------------------------------------
!                                  ERROR
    case(60)
!      select case(isolver)
!        case(1) ; call solve1(         NR_RHS_PROB)
!        case(2) ; call mumps_solve_seq(NR_RHS_PROB)
!        case(3) ; call uhm_solve
!                  call uhm_solver_flush(UHM_SOLVER_PTR)
!       endselect
      call exact_error
!    case(61) ; call dumpout_error_to_file(FILE_ERR)
    case(70) ; call compute_residual
    endselect
!
! end infinite loop
  enddo
!
! finalize library
  call finalize
!
end program main