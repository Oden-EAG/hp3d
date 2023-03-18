!----------------------------------------------------------------------
!
!     program name      - main
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!     purpose:          - main driver for MPI Test Program
!                         Maxwell Galerkin implementation
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
   use MPI        , only: MPI_COMM_WORLD,MPI_BARRIER,MPI_Wtime, &
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
   real(8) :: start_time, end_time
!
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
   call set_environment
   call end_environment
!
   if (RANK .eq. ROOT) then
!  ...print header
      write(6,*)
      write(6,*) '//                              //'
      write(6,*) '// --  MPI MAXWELL GALERKIN  -- //'
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
!..FLAGS
   ISTC_FLAG = .true.
   STORE_STC = .true.
   HERM_STC  = .false.
!
   if (RANK .eq. ROOT) then
      write(*,*) 'FLAGS:'
      write(*,*) ' ISTC_FLAG: ', ISTC_FLAG
      write(*,*) ' STORE_STC: ', STORE_STC
      write(*,*) ' HERM_STC : ', HERM_STC
      write(*,*)
   endif
!
!..determine number of omp threads running
 1025 format(A,I2)
   if (RANK .eq. ROOT) then
      write(6,1025) ' Initial polynomial order: ',IP
   endif
#if HP3D_USE_OPENMP
   if (RANK .eq. ROOT) then
   !$OMP parallel
   !$OMP single
      num_threads = omp_get_num_threads()
      write(6,1025) ' Number of OpenMP threads: ',num_threads
   !$OMP end single
   !$OMP end parallel
   endif
#endif
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
   integer :: idec, i, r, lb, count, src
!
!----------------------------------------------------------------------
!
   if (RANK .ne. ROOT) then
      write(*,*) 'master_main: RANK .ne. ROOT'
      stop
   endif
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
!
!..start user interface, with idec
!..broadcast user command to workers
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
!
!  ...broadcasting to worker procs
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
   use common_prob_data
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
!  ...receiving from master proc
      count = 1; src = ROOT
      call MPI_BCAST (idec,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
!
      select case(idec)
!     ...QUIT
         case(0) ; goto 99
!
!     ...HP3D graphics (do not use with distributed mesh)
         case(1,2) ;
!
!     ...Paraview graphics
         case(3) ; call exec_case(idec)
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
         case(20,21,22,23,26)
            call exec_case(idec)
!
!     ...MPI Routines
         case(31,33,35)
            call exec_case(idec)
!
!     ...MPI Routines (partitioner)
         case(30,32)
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
99 continue
!
end subroutine worker_main
