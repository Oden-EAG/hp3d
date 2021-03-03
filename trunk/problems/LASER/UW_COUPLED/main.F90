!----------------------------------------------------------------------
!
!     program name      - main
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!     purpose:          - main driver for MPI Test Program
!                         Poisson Galerkin implementation
!
!----------------------------------------------------------------------
!
program main
!
   use control
   use data_structure3D
   use environment
   use GMP
   use physics
   use commonParam
   use laserParam
!
   use assembly
   use assembly_sc, only: IPRINT_TIME
   use stc        , only: STORE_STC,HERM_STC
!
   use MPI        , only: MPI_COMM_WORLD,MPI_BARRIER, &
                          MPI_GET_PROCESSOR_NAME,MPI_MAX_PROCESSOR_NAME
   use mpi_param  , only: ROOT,RANK,NUM_PROCS
   use mpi_wrapper, only: mpi_w_init,mpi_w_finalize
!
   implicit none
!
!..auxiliary variables
   integer :: i, iargs, ierr, req, ret, plen
!
   !integer :: flag(6)
   !integer :: physNick
!
!..OMP variables
   integer :: num_threads, omp_get_num_threads
!
!..MPI variables
   character(MPI_MAX_PROCESSOR_NAME) :: pname
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
!
   character      :: arg
   character(32)  :: args
   character(512) :: cmd
!
!----------------------------------------------------------------------
!
!..Initialize MPI environment
   call mpi_w_init
!
   if (RANK .eq. ROOT) then
      write(*,*) '    =========================    '
      write(*,*) '    Program command and args     '
      write(*,*) '    =========================    '
!  ...read command
      call get_command(cmd)
      write(*,*) cmd
!  ...read number of arguments
      iargs = command_argument_count()
!  ...read arguments
      do i=0,iargs
         call get_command_argument(i, args)
         write(*,*) args
      end do
      write(*,*) '    =========================    '
   endif
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..Set common hp3D environment parameters (reads in options arguments)
   call begin_environment
   call set_environment_laser
   call end_environment
!
!..add flags 6,7,8 to Dirichlet flag list
   call add_dirichlet_to_list(6)
   call add_dirichlet_to_list(7)
   call add_dirichlet_to_list(8)
!
!..PML
   call set_PML
!
   if (RANK .eq. ROOT) then
!  ...print header
      write(6,*)
      write(6,*) '//                              //'
      write(6,*) '//  -- MPI LASER AMPLIFIER  --  //'
      write(6,*) '//                              //'
      write(6,*)
#if DEBUG_MODE
      write(*,*) '    =========================    '
      write(*,*) '      RUNNING in DEBUG_MODE      '
      write(*,*) '    =========================    '
#endif
   endif
   flush(6)
!
   IPRINT_TIME = 1
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..initialize physics, geometry, etc.
   call MPI_GET_PROCESSOR_NAME (pname,plen,ierr);
   do i = 0, NUM_PROCS-1
      if ((RANK .eq. i) .and. (RANK .eq. ROOT)) then
         write(6,*)
         write(6,1020) "Master proc [",RANK,"] on node [",trim(pname),"]: initialize..."
         QUIET_MODE = .FALSE.
      else if ((RANK .eq. i) .and. (RANK .ne. ROOT)) then
         write(6,1020) "Worker proc [",RANK,"] on node [",trim(pname),"]: initialize..."
         QUIET_MODE = .TRUE.
      endif
   enddo
 1020 format (A,I4,A,A,A)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   call initialize
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
   if (RANK .eq. ROOT) write(6,1015) end_time-start_time
 1015 format(' initialize : ',f12.5,' seconds',/)
!
!..determine number of omp threads running
   if (RANK .ne. ROOT) goto 80
!
!..print problem parameters
   if ((GEOM_NO .eq. 1) .and. (ZL .gt. 1.d0)) then ! rectangular waveguide
      write(*,9001) ' Wavelengths/Unit length  = ', sqrt(OMEGA*OMEGA-PI*PI)/(2.d0*PI)
      write(*,9000) ' Waveguide length         = ', ZL
   endif
   if (GEOM_NO .eq. 5) then ! fiber waveguide
      write(*,9000) ' Fiber length             = ', ZL
      write(*,9000) ' Signal frequency         = ', OMEGA_SIGNAL
      write(*,9000) ' Wavelengths/Unit length  = ', 1.d0/(LAMBDA_SIGNAL/REF_INDEX_CORE)
      write(*,9000) ' Numerical Aperture       = ', NA
      write(*,9000) ' V-number                 = ', VNUM
      write(*,9000) ' Core ref index           = ', CORE_NX
      write(*,9000) ' Clad ref index           = ', CLAD_NX
   endif
   if (ANISO_REF_INDEX .eq. 1) then
      write(*,9010) ' ANISO_REF_INDEX          = ', ANISO_REF_INDEX
      write(*,9000) ' CORE_NY                  = ', CORE_NY
      write(*,9000) ' CLAD_NY                  = ', CLAD_NY
   endif
   if (ART_GRATING .eq. 1) then
      write(*,9010) ' ART_GRATING              = ', ART_GRATING
   endif
   if (USE_PML) then
      write(*,9000) ' PML_REGION               = ', PML_REGION
   endif
   if (NONLINEAR_FLAG .eq. 1) then
      write(*,9020) ' Raman gain               = ', RAMAN_GAIN
      write(*,9020) ' Active gain              = ', ACTIVE_GAIN
      write(*,9010) ' COPUMP                   = ', COPUMP
      write(*,9010) ' FAKE_PUMP                = ', FAKE_PUMP
   endif
   write(*,9030) ' Polynomial order (x,y,z) = ', ORDER_APPROX_X,ORDER_APPROX_Y,ORDER_APPROX_Z
   write(*,9010) ' ISOL                     = ', ISOL
   write(*,9010) ' NEXACT                   = ', NEXACT
   write(*,9010) ' FAST INTEGRATION         = ', FAST_INT
   write(*,9015) ' OUTPUT_DIR               = ', trim(OUTPUT_DIR)
   if (HEAT_FLAG .eq. 1) then
      write(*,9010) ' NSTEPS                   = ', NSTEPS
      write(*,9000) ' DELTA_T                  = ', DELTA_T
      if (ANISO_HEAT .eq. 1) then
         write(*,9020) ' ALPHA_Z                  = ', ALPHA_Z
      endif
   endif
 9000 format(A,F11.6)
 9001 format(A,F11.3)
 9010 format(A,I3)
 9015 format(A,A)
 9020 format(A,ES11.2)
 9030 format(A,' (',I1,',',I1,',',I1,') ')
!
!$OMP parallel
!$OMP single
   num_threads = omp_get_num_threads()
   write(6,1025) ' Number of OpenMP threads: ',num_threads
 1025 format(A,I2)
!$OMP end single
!$OMP end parallel
!
   80 continue
!
!..set interface variables
!  (1) - H1 field for heat (1 component)
!  (2) - Hcurl for Maxwell trace for signal (2 components)
!  (3) - Hcurl for Maxwell trace for pump   (2 components)
!  (4) - Hdiv trace for heat (1 component)
!  (5) - L2 field for Maxwell (signal, 6 components)
!  (6) - L2 field for Maxwell (pump  , 6 components)
   PHYSAi(1:6) = (/.false.,.true.,.true.,.true.,.false.,.false./)
!
!..set static condensation flags
   ISTC_FLAG = .true.
   STORE_STC = .true.
   HERM_STC = .true.
!
   if (HERM_STC) then
      arg = 'H'
   else
      arg = 'G'
   endif
!
   if (RANK .eq. ROOT) then
      write(*,*) 'FLAGS:'
      write(*,*) ' ISTC_FLAG: ', ISTC_FLAG
      write(*,*) ' STORE_STC: ', STORE_STC
      write(*,*) ' HERM_STC : ', HERM_STC
      write(*,*) ' PHYSAi   : ', PHYSAi
      write(*,*)
   endif
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..Maxwell signal solve
   NO_PROBLEM = 3
   ! physNick = 1; flag=0; flag(5)=1
   PHYSAm(1:6) = (/.false.,.true.,.false.,.false.,.true.,.false./)
!
   if (JOB .ne. 0) then
      if (NONLINEAR_FLAG .eq. 0) then
         if (HEAT_FLAG .eq. 0) then
            call exec_job           ! Linear Maxwell
         else
            call exec_job_heat      ! Linear Heat
         endif
      else
         if (HEAT_FLAG .eq. 0) then
            call exec_job_nl        ! Nonlinear gain fiber
         else
            call exec_job_coupled   ! Coupled Heat/Maxwell
         endif
      endif
   else
      if (RANK .eq. 0) then
         call master_main
      else
         call worker_main
      endif
   endif
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   if (RANK.eq.ROOT) write(*,*) 'Back in main. Just finishing up...'
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
   use commonParam
   use data_structure3D
   use GMP
!
   use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
   use mpi_param     , only: ROOT,RANK,NUM_PROCS
   use par_mesh      , only: DISTRIBUTED
   use zoltan_wrapper, only: zoltan_w_set_lb
!
   implicit none
!
!..MPI variables
   integer :: ierr
!
!..auxiliary variables
   integer :: idec, i, r, lb, count, src
   character(len=8) :: filename
   integer :: mdle,nr_elem_ref,kref
   real(8) :: res
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
      write(*,*) '         ----    I/O    ----             '
      write(*,*) 'Paraview ...............................3'
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
      write(*,*) 'Anisotropic h-refinements (z)..........23'
      write(*,*) '                                         '
      write(*,*) '        ---- MPI Routines ----           '
      write(*,*) 'Distribute mesh........................30'
      write(*,*) 'Collect dofs on ROOT...................31'
      write(*,*) 'Evaluate mesh partition (Zoltan).......33'
      write(*,*) 'Run verification routines..............35'
      write(*,*) '                                         '
      write(*,*) '          ---- Solvers ----              '
      write(*,*) 'MUMPS (MPI, Nested)....................40'
      write(*,*) 'MUMPS (MPI)............................41'
      write(*,*) 'MUMPS (OpenMP).........................42'
      write(*,*) 'Pardiso (OpenMP).......................43'
      write(*,*) 'Frontal (Seq)..........................44'
      write(*,*) 'PETSc (MPI)............................45'
      write(*,*) '                                         '
      write(*,*) '     ---- Error and Residual ----        '
      write(*,*) 'Compute exact error....................50'
      write(*,*) 'Compute residual.......................51'
      write(*,*) '                                         '
      write(*,*) '            ---- Misc ----               '
      write(*,*) 'Compute Power .........................60'
      write(*,*) 'Compute Temperature ...................61'
      write(*,*) '                                         '
      write(*,*) '          ---- Debugging ----            '
      write(*,*) 'Refine a single element................70'
      write(*,*) 'Random refinements.....................71'
      write(*,*) 'Read HIST file.........................72'
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
!     ...Paraview
         case(3)
            call exec_case(idec)
!
!     ...Print data structure
         case(10,11)
            r = ROOT
            if (NUM_PROCS > 1) then
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
         case(20,21,22,23)
            call exec_case(idec)
!
!     ...MPI Routines
         case(31,33,35)
            call exec_case(idec)
!
!     ...MPI Routines (partitioner)
         case(30)
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
         case(50,51)
            call exec_case(idec)
!
!     ...Miscellanea
         case(60,61)
            call exec_case(idec)
!
!     ...Debugging routines
         case(70,71)
            call exec_case(idec)
!
         case(72)
            if (NUM_PROCS > 1) then
               write(*,*) 'cannot use NHIST for MPI currently. returning...'
               cycle
            endif
            ! WRITE FOR DEBUGGING
            filename='HIST.dat'
            open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="UNKNOWN",ACTION="READ")
            do
               write(*,*) 'A. reading from file and refining'
               read(UNIT=9, FMT="(I6)") nr_elem_ref
               write(*,*) '   nr_elem_ref = ', nr_elem_ref
               do i=1,nr_elem_ref
                  read(UNIT=9, FMT="(I6)") mdle
                  select case (NODES(mdle)%type)
                     case('mdlb'); kref = 110
                     case('mdlp'); kref = 10
                     case default
                        write(*,*) 'READING UNEXPECTED ELEMENT TYPE (mdle): ',NODES(mdle)%type,' (',mdle,')'
                        call pause
                  end select
                  call refine(mdle,kref)
               enddo
               call pause
               write(*,*) 'B. calling close_mesh'
               call close_mesh
               call pause
               write(*,*) 'C. calling update_gdof'
               call update_gdof
!               call pause
               write(*,*) 'D. calling update_Ddof'
               call update_Ddof
               write(*,*) 'solve?  1 = YES; 0 = NO'
               read(*,*) i
               if (i .eq. 1) then
                  write(*,*) 'E. calling pardiso_sc'
                  call pardiso_sc('H')
                  write(*,*) '   computing residual'
                  call residual(res)
               endif
               write(*,*) 'continue?  1 = YES; 0 = NO'
               read(*,*) i
               if (i .eq. 0) exit
            enddo
            close(UNIT=9)
            ! END WRITE FOR DEBUGGING
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
!----------------------------------------------------------------------
! worker_main
!----------------------------------------------------------------------
subroutine worker_main()
!
   use environment
   use commonParam
   use data_structure3D
   use GMP
!
   use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
   use mpi_param     , only: ROOT,RANK,NUM_PROCS
   use par_mesh      , only: DISTRIBUTED
   use zoltan_wrapper, only: zoltan_w_set_lb
!
   implicit none
!
!..MPI variables
   integer :: ierr
!
!..auxiliary variables
   integer :: idec, i, r, lb, count, src
!
!----------------------------------------------------------------------
!
   if (RANK .eq. ROOT) then
      write(*,*) 'worker_main: RANK .eq. ROOT'
      stop
   endif
!
!..test accessing data structures
   write(6,9020) '[', RANK, '] : ', 'NRELIS,NRELES,NRNODS = ',NRELIS,NRELES,NRNODS
 9020 format(A,I3,A,A,I4,', ',I4,', ',I4)
!
   flush(6)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!
!..receive broadcast from master on how to proceed
!..do that in a loop, using idec
!..if receiving '0', end worker_main
   idec = 1
   do while(idec /= 0)
!
      write(6,9030) '[', RANK, '] : ','Waiting for broadcast from master...'
      count = 1; src = ROOT
      call MPI_BCAST (idec,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
      write(6,9010) '[', RANK, '] : ','Broadcast: idec = ', idec
 9010 format(A,I3,A,A,I3)
!
      select case(idec)
!     ...QUIT
         case(0) ; goto 99
!
!     ...Paraview
         case(3)
            call exec_case(idec)
!
!     ...Print data structure
         case(10,11)
            count = 1; src = ROOT
            call MPI_BCAST (r,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
            if (r .eq. RANK) then
               call exec_case(idec)
            endif
         case(15,16,17)
            call exec_case(idec)
!
!     ...Refinements
         case(20,21,22,23)
            call exec_case(idec)
!
!     ...MPI Routines
         case(31,33,35)
            call exec_case(idec)
!
!     ...MPI Routines (partitioner)
         case(30)
            if (DISTRIBUTED) then
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
         case(50,51)
            call exec_case(idec)
!
!     ...Miscellanea
         case(60,61)
            call exec_case(idec)
!
!     ...Debugging routines
         case(70,71)
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
99 continue
   write(6,9030) '[', RANK, '] : ','worker_main end.'
 9030 format(A,I4,A,A)
!
end subroutine worker_main
