!----------------------------------------------------------------------
!                                                                     
!     program name      - main
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 2019
!                                                                     
!     purpose:          - main driver for MPI Poisson Test Program
!                                                                    
!----------------------------------------------------------------------
!    
program main
!
   use environment
   use common_prob_data
   use data_structure3D
   use GMP
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
   call begin_environment  ! <-- found inside src/modules/environment.F90
!
!..This sets default values for: FILE_REFINE,FILE_VIS,VLEVEL,
!                                FILE_CONTROL,FILE_GEOM,FILE_ERR,
!                                FILE_HISTORY,FILE_PHYS
   call set_environment  ! <-- found inside ../common/set_environment.F90
!
!..Exit if this is a "dry run"
   call end_environment  ! <-- found inside src/modules/environment.F90
!
   if (RANK .eq. ROOT) then
!  ...print header
      write(6,*)
      write(6,*) '//                          //'
      write(6,*) '// --  MPI TEST PROGRAM  -- //'
      write(6,*) '//                          //'
      write(6,*)
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
   if (RANK .eq. ROOT) then
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
   use MPI      , only: MPI_COMM_WORLD,MPI_INTEGER
   use mpi_param, only: ROOT,RANK,NUM_PROCS
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
   integer :: idec, i, r, count, src
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
      write(*,*) 'Print data structure arrays ...........11'
      write(*,*) '                                         '
      write(*,*) '        ---- Refinements ----            '
      write(*,*) 'Single uniform h-refinement............20'
      write(*,*) 'Single uniform p-refinement............21'
      write(*,*) '                                         '
      write(*,*) '        ---- MPI Routines ----           '
      write(*,*) 'Distribute mesh........................30'
      write(*,*) 'Collect dofs on ROOT...................31'
      write(*,*) 'Run verification routines..............35'
      write(*,*) '                                         '
      write(*,*) '          ---- Solvers ----              '
      write(*,*) 'MUMPS (OpenMP).........................40'
      write(*,*) 'MUMPS (MPI)............................41'
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
!
      read( *,*) idec
      write(6,8010) '[', RANK, '] : ','Broadcast: idec = ', idec
      count = 1; src = ROOT
      call MPI_BCAST (idec,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
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
!
!     ...Refinements
         case(20,21)
            call exec_case(idec)
!
!     ...MPI Routines
         case(30,31,35)
            call exec_case(idec)
!
!     ...Solvers
         case(40,41)
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
   use MPI      , only: MPI_COMM_WORLD,MPI_INTEGER
   use mpi_param, only: ROOT,RANK,NUM_PROCS
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
   integer :: idec, i, r, count, src
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
      flush(6)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
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
!
!     ...Refinements
         case(20,21)
            call exec_case(idec)
!
!     ...MPI Routines
         case(30,31,35)
            call exec_case(idec)
!
!     ...Solvers
         case(40,41)
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
!
   implicit none
!
   integer, intent(in) :: idec
!
!----------------------------------------------------------------------
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
!  ...run mesh verification routines
      case(35)
         write(*,*) 'verify distributed mesh consistency...'
         call par_verify
!
!  ...solve problem with omp_mumps (OpenMP MUMPS)
      case(40)
         write(*,*) 'calling OpenMP MUMPS solver...'
         call mumps_sc('G')
!
!  ...solve problem with par_mumps (MPI MUMPS)
      case(41)
         write(*,*) 'calling MPI MUMPS solver...'
         call par_mumps_sc('G')
!
      case default
         write(*,*) 'exec_case: unknown case...'
   end select
!
end subroutine exec_case

