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
   integer :: i, ierr, req, ret
!
!----------------------------------------------------------------------
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
      write(6,*)
      write(6,*) '//                          //'
      write(6,*) '// --  MPI TEST PROGRAM  -- //'
      write(6,*) '//                          //'
      write(6,*)
      IPRINT_TIME = 1
   endif
!
   flush(6)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
   write(6,1010) "Hello world! My MPI RANK is ", RANK, " out of ", NUM_PROCS, " processes."
 1010 format (A,I3,A,I3,A)
   flush(6)
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..printing in order
!..this does not guarantee printing in order
!..but it is "more likely" to print in order
!..how "fast" print statements are written to stdout depends
!..on the mpi implementation and runtime environment
   do i = 0, NUM_PROCS-1
      if (RANK == i .and. RANK == ROOT) then
         write(6,*)
         write(6,1020) "Master proc [", RANK, "], initialize.."
         QUIET_MODE = .FALSE.
         call initialize
      else if (RANK == i) then
         write(6,1020) "Worker proc [", RANK, "], initialize.."
         QUIET_MODE = .TRUE.
         call initialize
         QUIET_MODE = .FALSE.
      else
      endif
      flush(6)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   enddo
 1020 format (A,I3,A,/)
!
!..FLAGS
   ISTC_FLAG = .true.
   STORE_STC = .true.
   HERM_STC  = .false.
!
   if (RANK .eq. 0) then
      call master_main
   else
      call worker_main
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
   use par_mesh      , only: DISTRIBUTED
   use zoltan_wrapper, only: zoltan_w_set_lb
!
   implicit none
!
!..MPI variables
   integer :: ierr
!
!..OMP variables
   integer :: num_threads, omp_get_num_threads
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
!..start user interface, with idec
!..broadcast user command to workers
!
!..determine number of omp threads running
!$OMP parallel
!$OMP single
      num_threads = omp_get_num_threads()
      write(6,8010) '[', RANK, '] : ','Number of OpenMP threads: ',num_threads
 8010 format(A,I3,A,A,I3)
!$OMP end single
!$OMP end parallel
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
      write(*,*) '    ---- Print Data Structure ----       '
      write(*,*) 'Print arrays (interactive).............10'
      write(*,*) 'Print data structure arrays............11'
      write(*,*) 'Print current partition (elems)........15'
      write(*,*) 'Print current subdomains (nodes).......16'
      write(*,*) '                                         '
      write(*,*) '        ---- Refinements ----            '
      write(*,*) 'Single uniform h-refinement............20'
      write(*,*) 'Single uniform p-refinement............21'
      write(*,*) 'Multiple uniform h-refs + solve........22'
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
      write(*,*) '                                         '
      write(*,*) '     ---- Error and Residual ----        '
      write(*,*) 'Compute exact error....................50'
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
!
      read( *,*) idec
      write(6,8010) '[', RANK, '] : ','Broadcast: idec = ', idec
      count = 1; src = ROOT
      call MPI_BCAST (idec,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
!
      select case(idec)
!     ...QUIT
         case(0) ; goto 89
!
!     ...Print data structure
         case(10,11)
            write(*,*) 'Select processor RANK: '
            read (*,*) r
            count = 1; src = ROOT
            call MPI_BCAST (r,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
            if (r .eq. RANK) then
               call exec_case(idec)
            endif
         case(15,16)
            call exec_case(idec)
!
!     ...Refinements
         case(20,21,22)
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
               read (*,*) lb
               count = 1; src = ROOT
               call MPI_BCAST (lb,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
               call zoltan_w_set_lb(lb)
            endif
            call exec_case(idec)
!
!     ...Solvers
         case(40,41,42,43)
            call exec_case(idec)
!
!     ...Error and Residual
         case(50)
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
!..OMP variables
   integer :: num_threads, omp_get_num_threads
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
!..determine number of omp threads running
!$OMP parallel
!$OMP single
      num_threads = omp_get_num_threads()
      write(6,9010) '[', RANK, '] : ','Number of OpenMP threads: ',num_threads
 9010 format(A,I3,A,A,I3)
!$OMP end single
!$OMP end parallel
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
!
      select case(idec)
!     ...QUIT
         case(0) ; goto 99
!
!     ...Print data structure
         case(10,11)
            count = 1; src = ROOT
            call MPI_BCAST (r,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
            if (r .eq. RANK) then
               call exec_case(idec)
            endif
         case(15,16)
            call exec_case(idec)
!
!     ...Refinements
         case(20,21,22)
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
         case(40,41,42,43)
            call exec_case(idec)
!
!     ...Error and Residual
         case(50)
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
 9030 format(A,I3,A,A)
!
end subroutine worker_main
!
!
!----------------------------------------------------------------------
! exec_case
!----------------------------------------------------------------------
subroutine exec_case(idec)
!
   use data_structure3D
   use par_mesh
   use common_prob_data
   use zoltan_wrapper, only: zoltan_w_partition,zoltan_w_eval
!
   implicit none
!
   integer, intent(in) :: idec
!
   integer :: i,nsteps,nstop
   logical :: solved
   integer :: mdle_subd(NRELES)
!
!----------------------------------------------------------------------
!
   solved = .false.
!
   select case(idec)
!
!  ...print data structure (interactive)
      case(10); call result
!
!  ...print general data structure info
      case(11)
         write(*,110) NRELIS,NRELES,NRNODS
 110     format(' NRELIS,NRELES,NRNODS            = ',3I10)
         write(*,111) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
 111     format(' NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ = ',4I10)
         write(*,112) MAXNODS,NPNODS
 112     format(' MAXNODS,NPNODS                  = ',2I10)
!
!  ...print current partition (elems)
      case(15)
         write(*,*) 'printing current partition (elems)...'
         call print_partition
!
!  ...print current subdomains (nodes)
      case(16)
         write(*,*) 'printing current subdomains (nodes)...'
         call print_subd
!
!  ...single uniform h-refinement
      case(20)
         write(*,*) 'global h-refinement...'
         call global_href
         call close_mesh
         call update_gdof
         call update_Ddof
!
!  ...single uniform p-refinement
      case(21)
         write(*,*) 'global p-refinement...'
         call global_pref
         call close_mesh ! not needed?
         call update_gdof
         call update_Ddof
!
!  ...Multi-step uniform h refinement
      case(22)
         nsteps=0
         do while (nsteps.le.0)
            write(*,*) 'Provide: number of uniform h-refinements'
            read(*,*) nsteps
         enddo
         do i=0,nsteps
!        ...solve first if needed
            if (.not. solved) then
               call par_mumps_sc('G')
               solved = .true.
            endif
!        ...display error and refine if necessary
            if (i.ne.nsteps) then
               call uniform_href(IUNIFORM,1,0.25d0, nstop)
               if (nstop.eq.1) then
                  write(*,*) 'No elements were refined.'
                  write(*,220) i
 220              format('Exiting loop after ',i2,' refinements...')
                  cycle
               else
                  solved = .false.
               endif
            else ! Last step only display (no refinement)
               call uniform_href(INOREFINEMENT,1,0.25d0, nstop)
            endif
        enddo
!
!  ...distribute mesh
      case(30)
         write(*,*) 'distribute mesh...'
         call distr_mesh
!
!  ...collect dofs on ROOT processor
      case(31)
         write(*,*) 'collecting dofs on ROOT...'
         call collect_dofs
!
!  ...suggest new mesh partition (Zoltan)
      case(32)
         if (DISTRIBUTED) then
            write(*,*) 'computing new mesh partition (Zoltan)...'
            call zoltan_w_partition(mdle_subd)
         else
            write(*,*) 'distribute mesh first to use Zoltan...'
         endif
!
!  ...evaluate current partition
      case(33)
         if (DISTRIBUTED) then
            write(*,*) 'evaluating current partition...'
            call zoltan_w_eval
         else
            write(*,*) 'distribute mesh first to use Zoltan...'
         endif
!
!  ...run mesh verification routines
      case(35)
         write(*,*) 'verify distributed mesh consistency...'
         call par_verify
!
!  ...solve problem with omp_mumps (OpenMP MUMPS)
      case(40)
         write(*,*) 'calling MUMPS (MPI) solver...'
         call par_mumps_sc('G')
!
!  ...solve problem with par_mumps (MPI MUMPS)
      case(41)
         write(*,*) 'calling MUMPS (OpenMP) solver...'
         call mumps_sc('G')
!
!  ...solve problem with pardiso (OpenMP)
      case(42)
         write(*,*) 'calling Pardiso (OpenMP) solver...'
         call pardiso_sc('G')
!
!  ...solve problem with Frontal solver (sequential)
      case(43)
         write(*,*) 'calling Frontal (Seq) solver...'
         call solve1(1)
!
      case(50)
         write(*,*) 'computing error and residual...'
         call exact_error
!
      case default
         write(*,*) 'exec_case: unknown case...'
   end select
!
end subroutine exec_case

