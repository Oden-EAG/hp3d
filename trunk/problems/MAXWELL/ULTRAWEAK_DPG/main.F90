!----------------------------------------------------------------------
!
!     program name      - main
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 2023
!
!     purpose:          - main driver for UW Maxwell application
!
!----------------------------------------------------------------------
program main
!
   use control
   use data_structure3D
   use environment
   use GMP
   use physics
   use commonParam
   use parametersDPG
!
   use assembly
   use assembly_sc, only: IPRINT_TIME
   use stc        , only: STORE_STC,HERM_STC
!
   use MPI        , only: MPI_COMM_WORLD,MPI_BARRIER,MPI_ABORT, &
                          MPI_GET_PROCESSOR_NAME,MPI_MAX_PROCESSOR_NAME
   use mpi_param  , only: ROOT,RANK,NUM_PROCS
   use mpi_wrapper, only: mpi_w_init,mpi_w_finalize
!
   implicit none
!
!..auxiliary variables
   integer :: i, ierr, req, ret, plen
!
!..OMP variables
#if HP3D_USE_OPENMP
   integer :: num_threads, omp_get_num_threads
#endif
!
!..MPI variables
   character(MPI_MAX_PROCESSOR_NAME) :: pname
!
!..timer
   real(8) :: MPI_Wtime, start_time, end_time
!
   character       :: arg
   character(1024) :: cmd
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
      write(*,*) trim(cmd)
      write(*,*) '    =========================    '
   endif
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..Set common hp3D environment parameters (reads in options arguments)
   call begin_environment
   call set_environment_maxwell
   call end_environment
!
   if (RANK .eq. ROOT) then
!  ...print header
      write(6,*)
      write(6,*) '//                              //'
      write(6,*) '//  ------- UW MAXWELL -------  //'
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
   write(*,9010) ' Initial polynomial order = ', IP
   write(*,9010) ' NORD_ADD (Delta p)       = ', NORD_ADD
   write(*,9010) ' ISOL                     = ', ISOL
   write(*,9010) ' NEXACT                   = ', NEXACT
   write(*,9010) ' IBCFLAG                  = ', IBCFLAG
   write(*,9020) ' ALPHA_NORM               = ', ALPHA_NORM
   write(*,9015) ' OUTPUT_DIR               = ', trim(OUTPUT_DIR)
 9000 format(A,F13.6)
 9001 format(A,F13.3)
 9010 format(A,I3)
 9015 format(A,A)
 9020 format(A,ES13.2)
!
#if HP3D_USE_OPENMP
   !$OMP parallel
   !$OMP single
      num_threads = omp_get_num_threads()
      write(6,1025) ' Number of OpenMP threads: ',num_threads
 1025 format(A,I2)
   !$OMP end single
   !$OMP end parallel
#endif
!
   80 continue
!
!..set interface variables
!  (1) - Hcurl for Maxwell trace (2 components)
!  (2) - L2 field for Maxwell (6 components)
   PHYSAi(1:2) = (/.true., .false./)
!
!..set homogeneous Dirichlet flags
   if (NEXACT.eq.0) then
      PHYSAd(1:2) = (/.false.,.false./)
   endif
!
   PHYSAm(1:2) = (/.true.,.true./)
!
!..set static condensation flags
   ISTC_FLAG = .true. ! activate automatic static condensation
   STORE_STC = .true. ! store Schur complement factors
   HERM_STC = .true.  ! assume Hermitian element matrix
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
      write(*,*) ' PHYSAd   : ', PHYSAd
      write(*,*)
   endif
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
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
      !write(*,*) 'Multiple uniform h-refs + solve........22'
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
      write(*,*) '          ---- Debugging ----            '
      write(*,*) 'Refine a single element................70'
      write(*,*) 'Random refinements.....................71'
      write(*,*) 'Read HIST file.........................72'
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
!
      read (*,*) idec
      !write(6,8010) '[', RANK, '] : ','Broadcast: idec = ', idec
 !8010 format(A,I3,A,A,I3)
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
         case(20,21)
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
         case(40,41,42,43,44,45,46)
            call exec_case(idec)
!
!     ...Error and Residual
         case(50,51)
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
                  select case (NODES(mdle)%ntype)
                     case(MDLB); kref = 110
                     case(MDLP); kref = 10
                     case default
                        write(*,*) 'READING UNEXPECTED ELEMENT TYPE (mdle): ',s_type(NODES(mdle)%ntype),' (',mdle,')'
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
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..receive broadcast from master on how to proceed
!..do that in a loop, using idec
!..if receiving '0', end worker_main
   idec = 1
   do while(idec /= 0)
!
      !write(6,9030) '[', RANK, '] : ','Waiting for broadcast from master...'
      count = 1; src = ROOT
      call MPI_BCAST (idec,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
      !write(6,9010) '[', RANK, '] : ','Broadcast: idec = ', idec
 !9010 format(A,I3,A,A,I3)
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
         case(20,21)
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
         case(40,41,42,43,44,45,46)
            call exec_case(idec)
!
!     ...Error and Residual
         case(50,51)
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
