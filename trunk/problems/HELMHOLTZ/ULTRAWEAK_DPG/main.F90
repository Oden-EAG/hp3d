!----------------------------------------------------------------------
!
!     program name      - main
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2020
!
!     purpose:          - main driver for MPI Test Program
!                         Poisson Galerkin implementation
!
!----------------------------------------------------------------------
!
program main
!
   use environment
   use common_prob_data_UW
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
!..OMP variables
   integer :: num_threads, omp_get_num_threads
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
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
      write(6,*) '//                                 //'
      write(6,*) '// --  MPI ACOUSTICS ULTRAWEAK  -- //'
      write(6,*) '//                                 //'
      write(6,*)
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
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time   = MPI_Wtime()
   if (RANK .eq. ROOT) write(6,1015) end_time-start_time
 1015 format(' initialize : ',f12.5,' seconds',/)
!
!..Set interface variables
!..(1)--trace
!..(2)--fluxv
!..(3)--L2 field (x4)
PHYSAi(1:3) = (/.true.,.true.,.false./)
PHYSAm(1:3) = (/.true.,.true.,.true./)
PHYSAd(1:3) = (/.false.,.false.,.false./)
!..FLAGS
   ISTC_FLAG = .true.
   STORE_STC = .true.
   HERM_STC  = .true.
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
   if (IBC_PROB.eq.3 .or. IBC_PROB.eq.4 .or. IBC_PROB.eq.6) call propagate_flag(2,3)
!
   if (JOB .ne. 0) then
      call exec_job
   else
      if (RANK .eq. 0) then
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
   use common_prob_data_UW
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
   use common_prob_data_UW
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

!..dump mesh to a txt file for MATLAB

subroutine matlab_dump_mesh(k)
!
   use parameters, only       : MAXbrickH
   use data_structure3D, only : NRELES, NODES

   implicit none
   integer, intent(in) :: k
   integer :: mdle, nordh,nordv,nord,i
   real*8, allocatable  :: xnod(:,:)
   character(len=20) :: str

!---------------------------------------------------------------------------------

   allocate(xnod(3,MAXbrickH)) ; xnod = 0.0d0
!
   open(12, file='output/MATLAB/meshk4unif_'//trim(str(k))//'.txt',       &
   status='replace',form='formatted',access='sequential')
   mdle=0;
   do i=1,NRELES
      call nelcon(mdle, mdle);
      call nodcor_vert(Mdle, xnod)
      write(12,*)                                                     &
      (/xnod(1,1),xnod(2,1),xnod(3,1), xnod(1,2),xnod(2,2),xnod(3,2), &
        xnod(1,3),xnod(2,3),xnod(3,3), xnod(1,4),xnod(2,4),xnod(3,4), &
        xnod(1,5),xnod(2,5),xnod(3,5), xnod(1,6),xnod(2,6),xnod(3,6), &
        xnod(1,7),xnod(2,7),xnod(3,7), xnod(1,8),xnod(2,8),xnod(3,8)/)
   enddo
!
!..dump out the order of approximation
   mdle = 0
   do i = 1, NRELES
      call nelcon(mdle,mdle)
      nord = NODES(mdle)%order
      call decode(nord, nordh,nordv)
      write(12,*) nordv
   enddo
   write(12,*) NRELES
   close(12)
   deallocate(xnod)
!
!
   end subroutine matlab_dump_mesh



!..Convert an integer to string

   character(len=20) function str(k)
!
   integer, intent(in) :: k
   write (str, *) k
   str = adjustl(str)
!
   end function str
