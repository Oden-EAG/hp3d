!----------------------------------------------------------------------
!                                                                     
!     program name      - main
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 2019
!                                                                     
!     purpose:          - main driver for MPI Test Program
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
   implicit none
!
include 'mpif.h'
!
   integer :: idec, i
!
!..OMP variables
   integer :: num_threads, omp_get_num_threads
!
!..MPI variables
   integer ierr, num_procs, rank
!
!----------------------------------------------------------------------
!
!..initialization
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
!..Before this line, everything is executed by every MPI proc
!
!..MPI
   call MPI_INIT ( ierr )
!..find out my rank (process id)
   call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
   write(*,1010) "Hello world! My MPI Rank is ", rank, " out of ", num_procs, " processes."
 1010 format (A,I2,A,I2,A)

   call MPI_Barrier (MPI_COMM_WORLD, ierr)
!
!..printing in order
!..this does not guarantee printing in order
!..but it is "more likely" to print in order
!..how "fast" print statements are written to stdout depends
!..on the mpi implementation and runtime environment
   do i = 0, num_procs-1
      if (rank == i .and. rank == 0) then
         write(*,1020) "Master proc [", rank, "]"
      else if (rank == i) then
         write(*,1020) "Slave  proc [", rank, "]"
      else
      endif
      call MPI_Barrier (MPI_COMM_WORLD, ierr)
   enddo
 1020 format (A,I2,A)
!
   if (rank .ne. 0) then
      goto 99
   endif
!
!..print fancy header
   write(*,*)'                      '
   write(*,*)'// --  MPI TEST PROGRAM  -- //'
   write(*,*)'                      '
!
!..Initialize common library 
!  (set common parameters, load solvers, and create initial mesh)
   call initialize  ! <-- found inside ../common/initialize.F90
!
!..determine number of omp threads running
!$OMP parallel
!$OMP single
      num_threads = omp_get_num_threads()
      write(*,1100) ' Number of OpenMP threads: ', num_threads
 1100 format(A,I2,/)
!$OMP end single
!$OMP end parallel
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
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      write(*,*) 'SELECT'
      write(*,*) 'QUIT ...................................0'
      write(*,*) '                                         '
      write(*,*) 'Print Data Structure arrays ...........10'
      write(*,*) '                                         '
      write(*,*) '         ---- Refinements ----           '
      write(*,*) 'Single Uniform h-refinement............20'
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='

      read( *,*) idec
!
      select case(idec)
!     ...QUIT
         case(0) ; goto 98
!
!-----------------------------------------------------------------------
!   ...DATA STRUCTURE
!-----------------------------------------------------------------------
!
!     ...print data structure
         case(10) ; call result
!         
!-----------------------------------------------------------------------
!  ...REFINEMENTS
!-----------------------------------------------------------------------
!
!     ...Single uniform refinement
         case(20)
            call global_href
            call close_mesh
            call update_gdof
            call update_Ddof
!
      end select
!
!..end infinite loop
   enddo
!
!..finalize library
   98 call finalize ! <-- found inside ../common/finalize.F90
!
   99 call MPI_FINALIZE ( ierr )
!..END MPI
!
   end program main
