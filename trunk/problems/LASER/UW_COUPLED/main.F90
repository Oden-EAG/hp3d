!---------------------------------------------------------
!  The code is run through a bash script
!---------------------------------------------------------
program main
!
   use environment
   use paraview
   use control
   use parametersDPG
   use GMP
   use data_structure3D
   use physics
   use CommonParam
   use LaserParam
   use assembly
   use m_assembly
!
   use MPI
!
   implicit none
!
!..declare variables
   real*8  :: err,rnorm,rvoid,factor,geo_err,geo_rnorm
   integer :: mdle,i,kref,idec,nreflag,istep,nstop
   integer :: numRef,iref,numLaserIter,nord,idec_solve,naniso,numProb
   integer, dimension(6) :: flag
   integer, allocatable, dimension(:) :: list_elem
   integer :: ic, iel, iso_ans, iii, info, ref_xyz
   integer :: physNick,chooseProbType
   logical :: ires,solved
   real*8  :: FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!
!..variables for nonlinear loop
   integer :: No,No1,No2,nonlin_steps
   real*8  :: L2NormDiff,stopEpsilon,scale
!
!..MPI variables
   integer ierr, num_procs, rank
!
!----------------------------------------------------------------------
!
!..initialize environment
!..set the variables for laser problem
   call begin_environment
   call set_environment_laser
   call end_environment
!
!..add flags 6,7,8 to Dirichlet flag list
   call add_dirichlet_to_list(6)
   call add_dirichlet_to_list(7)
   call add_dirichlet_to_list(8)
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
!..how "fast" print statements are written to stdout
!..depends on the mpi implementation and runtime environment
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
!..display menu
   idec = 1
   flag = 0
!..flags
   solved = .false.
   ires = .true.
!
!..IPARATTR
!..this parameter specifies which fields we dump out to paraview:
!     - 1 h1 (heat)
!     - 2 hcurl (signal, E and H flux)
!     - 2 hcurl (pump,   E and H flux)
!     - 1 hdiv (heat flux)
!     - 6 l2 (signal, E and H field)
!     - 6 l2 (pump,   E and H Field)
!
   allocate(IPARATTR(6))
   IPARATTR = (/0,0,0,0,0,0/)
   !IPARATTR = (/0,0,0,0,6,0/)
   !IPARATTR = (/0,0,0,0,6,6/)
!
!..PML
   call set_PML
   !STORE_STC = .FALSE.
   STORE_STC = .TRUE.
!
!..print fancy header
   write(*,*)'                   '
   write(*,*)'                   '
   write(*,*)'                   '
   write(*,*)'// ------ THE LASER PROBLEM ------ //'
   write(*,*)'                   '
   write(*,*)'                   '
   write(*,*)'                   '
!
!..initialize problem
   call initialize
   write(*,9000) ' OMEGA_RATIO_SIGNAL       = ', OMEGA_RATIO_SIGNAL
   write(*,9000) ' OMEGA_RATIO_PUMP         = ', OMEGA_RATIO_PUMP
   write(*,9000) ' OMEGA_RATIO              = ', OMEGA_RATIO
   write(*,9000) ' Numerical Aperture       = ', NA
   write(*,9000) ' (REF_INDEX_CORE**2-1.d0) = ', (REF_INDEX_CORE**2-1.d0)
   write(*,9000) ' V-number                 = ', VNUM
   write(*,9000) ' PML_REGION               = ', PML_REGION
   write(*,9010) ' NEXACT                   = ', NEXACT
   write(*,9010) ' NSTEPS                   = ', NSTEPS
 9000 format(A,F10.6)
 9010 format(A,I3)
!
!..check if LASER_MODE and NONLINEAR_FLAG are compatible:
   if((LASER_MODE.ne.0).and.(NONLINEAR_FLAG.ne.1)) then
      write(*,*) 'error from main: ', &
         'NONLINEAR_FLAG must be 1 when running LASER_MODE. stop.'
      stop
   endif
   write(*,*)
!
!..display menu in infinite loop
 10 continue
!
!     -- INTERACTIVE MODE --
!
!     display menu in infinite loop
!
   write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   write(*,*) 'Quit ...................................0'
   write(*,*) '                                         '
   write(*,*) 'Geometry graphics (X11) ................1'
   write(*,*) 'HP3D graphics (X11) ....................2'
   write(*,*) 'Paraview ...............................3'
   write(*,*) 'Print Data Structure arrays ............5'
   write(*,*) 'DumpOut Data Structure .................6'
!   write(*,*) 'DumpIn Data Structure  .................7'
   write(*,*) '                                         '
   write(*,*) '    ----------  Geometry  ----------     '
   write(*,*) 'Geometry error ........................10'
!   write(*,*) 'BC interpolation error ................11 (disabled)'
!   write(*,*) 'BC interpolation error rates ..........12'
   write(*,*) '                                         '
   write(*,*) '            --  Refinements  --          '
   write(*,*) 'Single Uniform  p-refinement...........15'
   write(*,*) 'Single Uniform  h-refinement ..........16'
   write(*,*) 'Single Adaptive h-refinement ..........17'
   write(*,*) 'Multi-step Uniform  h-refinements .....18'
   write(*,*) 'Multi-step Adaptive hp-refinements ....19'
   write(*,*) 'Interactive h-refinement ..............20'
   write(*,*) 'M uniform, N anisotropic refs+solve ...21'
!   write(*,*) 'Refine face of prismatic core .........22'
   write(*,*) '0 uniform and N aniso refs in z .......23'
   write(*,*) '                                         '
   write(*,*) '           ---- Solvers ----             '
   write(*,*) 'Solve (frontal) .......................40'
   write(*,*) 'Solve (MUMPS)   .......................41'
   write(*,*) 'Solve (PARDISO) .......................42'
   write(*,*) '                                         '
   write(*,*) '            ---- Error ----              '
   write(*,*) 'Compute Field Error  ..................50'
   write(*,*) 'Compute norm of current solution.......51'
   write(*,*) '                                         '
   write(*,*) '            ---- Misc ----               '
!   write(*,*) 'Test Bessel Functions..................60'
   write(*,*) 'Compute Power .........................61'
   write(*,*) 'Test Copy Coms ........................63'
   write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   read (*,*) idec
   write(*,*) 'SIGMA is ', SIGMA
!
!----------------------------------------------------------------------
!
   select case(idec)
!  ...quit
      case(0) ; goto 98
!
!-----------------------------------------------------------------------
!                                 GRAPHICS
!-----------------------------------------------------------------------
!
!  ...GMP graphics
      case(1) ; call graphg
!
!  ...hp3d graphics
      case(2) ; call graphb
!
!  ...Paraview graphics
      case(3) ; call paraview_driver(IPARATTR)
!
!-----------------------------------------------------------------------
!                           DATA STRUCTURE
!-----------------------------------------------------------------------
!
      case(5) ; call result
!
!  ...dump out data structure
      case(6)
         write(*,*)'dumping out GMP...'
         call dumpout_GMP
         write(*,*)'dumping out physics...'
         call dumpout_physics
         write(*,*)'dumping out HP3D...'
         call dumpout_hp3d('./files/dumpc3Dhp')

!  ...dump in data structure
!     case(7)
!         write(*,*)'dumping in GMP...'
!         call dumpin_GMP('./files/dumpGMP')
!         write(*,*)'dumping in physics...'
!         call dumpin_physics_from_file('./files/dumpPHYS')
!         write(*,*)'dumping in HP3D...'
!         call dumpin_hp3d('./files/dumpc3D')
!
!-----------------------------------------------------------------------
!                            GEOMETRY ERROR
!-----------------------------------------------------------------------
!
!  ...geometry error
      case(10) ; call geom_error(err,rnorm)
      !case(11) ; call compute_BC_interp_error
!
!  ...rate test geometry error
      case(12)
         write(*,*) 'Testing h-convergence rates for geometry error'
!     ...Get #refinements to be done
         write(*,*) 'Enter number of uniform H-refinements:'
         read(*,*) numRef
         iref = 0
         do while(iref.lt.numRef)
!        ...Compute geometry error
            geo_rnorm = 0.d0
            geo_err = 0.d0
            call geometry_error(geo_err,geo_rnorm)
!        ...Do uniform h-refinements
            call global_href
            call close_mesh
            call update_gdof
            call update_Ddof
            iref = iref+1
         enddo
!
!-----------------------------------------------------------------------
!                            REFINEMENTS
!-----------------------------------------------------------------------
!
!  ...Uniform p refinement
      case(15)
         mdle=0
         do iel=1,NRELES
            call nelcon(mdle, mdle)
            nord = NODES(mdle)%order
            select case(NODES(mdle)%type)
               case('mdln','mdld'); nord = nord+1
               case('mdlp'); nord = nord+11
               case('mdlb'); nord = nord+111
            end select
            call nodmod(mdle, nord)
         enddo
!
!  ...Single uniform refinement
      case(16)
         if (.not.solved) write(*,*) 'You have not solved.'
         if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
         call select_phys_problem(NO_PROBLEM)
         select case(NO_PROBLEM)
            case(1,2)
               physNick = 1000
               flag=0;flag(1)=1
            case(3)
               physNick = 1
               flag=0;flag(5)=1
            case(4)
               physNick = 1
               flag=0;flag(6)=1
         end select
         call refine_DPG(IUNIFORM,1,0.25d0,flag,physNick,ires,0, nstop)
         solved=.FALSE.
!
!  ...Single adaptive refinement
      case(17)
         if (.not.solved) then
            write(*,*) 'You have not solved. Cannot adaptively refine.'
         else
            call select_phys_problem(NO_PROBLEM)
            select case(NO_PROBLEM)
               case(1,2)
                  physNick = 1000
                  flag=0;flag(1)=1
               case(3)
                  physNick = 1
                  flag=0;flag(5)=1
               case(4)
                  physNick = 1
                  flag=0;flag(6)=1
            end select
            nreflag=0
            do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
               write(*,*) 'REFINEMENT FLAG:h-refine=1,p-refine=2,hp-refine=3'
               write(*,*) '0.d0<FACTOR<1.d0'
               write(*,*) 'Provide: REFINEMENT FLAG, FACTOR, NO_PROBLEM'
               read(*,*) nreflag,factor
            enddo
            call refine_DPG(IADAPTIVE,nreflag,factor,flag,physNick,ires,0, nstop)
            if (nstop.eq.1) write(*,*) 'No elements were refined.'
            if (nstop.eq.0) solved=.FALSE.
         endif
!
!     ...Multi-step uniform h refinement
         case(18)
!        ...select number of refinements and the solver
            write(*,*) 'NUMBER OF REFINEMENTS > 0'
            write(*,*) 'SOLVER:  1 = frontal'
            write(*,*) '         2 = MUMPS'
            write(*,*) '         3 = PARDISO'
            write(*,*) '         4 = OMP MUMPS'
            write(*,*) '         5 = OMP MUMPS 64'
            write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
            read(*,*) nsteps,idec_solve
!        ...select problem
            call select_phys_problem(NO_PROBLEM)
            call set_physAm(NO_PROBLEM, physNick,flag)
!        ...begin refinements
            do i=0,nsteps
!           ...solve first if needed
               if (.not.solved) then
                  select case(idec_solve)
                     case(1)
                        call solve1(MY_NR_RHS)
                     case(2)
                        call mumps_solve_seq(MY_NR_RHS)
                     case(3)
                        IPRINT_TIME = 1
                        call pardiso_sc('H')
                        IPRINT_TIME = 0
                     case(4)
                        IPRINT_TIME = 1
                        call mumps_sc('H')
                        IPRINT_TIME = 0
                     case(5)
                        write(*,*) 'mumps_sc_64 not yet configured.'
                        !call mumps_sc_64('H')
                  end select
               endif
!           ...say it has solved and save results to paraview file
               solved=.TRUE.
!           ...display error and refine if necessary
               if (i.ne.nsteps) then
                  call refine_DPG(IUNIFORM,1,0.25d0,flag,physNick,ires,0, nstop)
                  if (nstop.eq.1) then
                     write(*,*) 'No elements were refined.'
                     write(*,1800) i
 1800                format('Exiting loop after ',i2,' refinements...')
                     cycle
                  else
                     solved=.FALSE.
                  endif
               else
!              ...Last step only display (no refinement)
                  call refine_DPG(INOREFINEMENT,1,0.25d0,flag,physNick,ires,0, nstop)
               endif
            enddo
!
!     ...Multi-step adaptive h refinement
         case(19)
            nreflag=0
            do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
               write(*,*) 'NUMBER OF REFINEMENTS>0'
               write(*,*) 'REFINEMENT FLAG:  1 =  h-refine'
               write(*,*) '                  2 =  p-refine'
               write(*,*) '                  3 = hp-refine'
               write(*,*) '0.d0 < FACTOR < 1.d0'
               write(*,*) 'SOLVER:  1 = frontal'
               write(*,*) '         2 = MUMPS'
               write(*,*) '         3 = PARDISO'
               write(*,*) '         4 = OMP MUMPS'
               write(*,*) 'Provide: NUMBER OF REFINEMENTS, ', &
                          'REFINEMENT FLAG, FACTOR, SOLVER'
               read(*,*) nsteps,nreflag,factor,idec_solve
!
               call select_phys_problem(NO_PROBLEM)
               call set_physAm(NO_PROBLEM, physNick,flag)
            enddo
!
            do i=0,nsteps
!           ...solve first if needed
               if (.not.solved) then
                  select case(idec_solve)
                     case(1)
                        call solve1(MY_NR_RHS)
                     case(2)
                        call mumps_solve_seq(MY_NR_RHS)
                     case(3)
                        call pardiso_sc('H')
                     case(4)
                        call mumps_sc('H')
                     case default
                        write(*,*) 'invalid SOLVER param. stop.'
                        stop
                  end select

               endif
!           ...say it has solved and save results to paraview file
               solved=.TRUE.
!
!           ...then display error and refine if necessary
               if (i.ne.nsteps) then
                  call refine_DPG(IADAPTIVE,nreflag,factor,flag,physNick,ires,0, nstop)
                  if (nstop.eq.1) then
                     write(*,*) 'No elements were refined.'
                     write(*,2100) i
                  cycle
               else
                  solved=.FALSE.
               endif
            else
!           ...Last step only display (no refinement)
               call refine_DPG(INOREFINEMENT,nreflag,factor,flag,physNick,ires,0, nstop)
         endif
      enddo
!  ...interactive refinements
      case(20)
         write(*,*)'Active elements:'
         mdle=0
         do i=1,NRELES
            call nelcon(mdle, mdle)
!
            select case(NODES(mdle)%type)
            case('mdlp') ; write(*,1900) mdle
            case('mdlb') ; write(*,1910) mdle
            case('mdln') ; write(*,1920) mdle
            case('mdld') ; write(*,1930) mdle
            endselect
 1900       format(' mdle = ',i6,' ; PRISM')
 1910       format(' mdle = ',i6,' ; BRICK')
 1920       format(' mdle = ',i6,' ; TET')
 1930       format(' mdle = ',i6,' ; PYRAMID')
!
         enddo
         call display_ref_kind
         write(*,1940)
 1940    format(' mdle,kref =')
         read(*,*) mdle,kref
!
!     ...refine element
         call refine(mdle,kref)
!
!     ...recover 1-irregular mesh, update geometry and Dirichlet dof's
         call close_mesh
         call update_gdof
         call update_ddof
!
!  ...M uniform, N anisotropic refinements
      case(21)
         nsteps=-1
         do while (nsteps.lt.0)
            write(*,*) 'M uniform, N anisotropic refinements'
            write(*,*) 'SOLVER:frontal=1,MUMPS=2,PARDISO=3,OMP MUMPS=4'
            write(*,*) 'Provide: NUMBER OF UNIFORM REFINEMENTS, SOLVER, NUMBER OF ANISOTROPIC REFINEMENTS and PROBLEM NUMBER'
            read (*,*) nsteps,idec_solve,naniso
            call select_phys_problem(NO_PROBLEM)
            call set_physAm(NO_PROBLEM, physNick,flag)
         enddo
!
         do i=0,nsteps
!        ...solve first if needed
            if (.not.solved) then
               select case(idec_solve)
                  case(1)
                     call solve1(MY_NR_RHS)
                  case(2)
                     call mumps_solve_seq(MY_NR_RHS)
                  case(3)
                     call pardiso_sc('H')
                  case(4)
                     call mumps_sc('H')
               end select
            endif
!        ...say it has solved and save results to paraview file
            solved=.TRUE.
!        ...display error and refine if necessary
            if (i.ne.nsteps) then
               call refine_DPG(IUNIFORM,1,0.25d0,flag,physNick,ires,0, nstop)
               if (nstop.eq.1) then
                  write(*,*) 'No elements were refined.'
                  write(*,2100) i
 2100             format('Exiting loop after ',i2,' refinements...')
                  cycle
               else
                  solved=.FALSE.
               endif
            else ! Last step only anisotropic refinements
               write(*,*) 'Last Uniform Refinement, will now to anisotropic refinement'
               do iii=1,naniso
                  call refine_DPG(IANISOTROPIC,1,0.25d0,flag,physNick,ires,1, nstop)
                  select case(idec_solve)
                     case(1)
                        call solve1(MY_NR_RHS)
                     case(2)
                        call mumps_solve_seq(MY_NR_RHS)
                     case(3)
                        call pardiso_sc('H')
                     case(4)
                        call mumps_sc('H')
                  end select
               enddo
               call refine_DPG(INOREFINEMENT,nreflag,0.25d0,flag,physNick,ires,0, nstop)
            endif
         enddo
!
!  ...anisotropic h-refinement of core face only
      case(22)
         write(*,*)'set number of uniform refinements in xy-direction'
         write(*,*)'-1 for infinite loop'
         read(*,*), iso_ans
         info=0
         if(iso_ans.ne.-1) then
!        ...ref_xyz = 1 for anisotropic refinements in xy-direction
            ref_xyz=1
            do iii=1,iso_ans
               call setAnisoRef(ref_xyz)
            enddo
            call update_gdof
            call update_Ddof
         else
!        ...ref_xyz = 1 for anisotropic refinements in xy-direction
            ref_xyz=1
            do while(1.gt.0)
               call setAnisoRef(ref_xyz)
            enddo
            call update_gdof
            call update_Ddof
         endif
!
!  ...anisotropic h-refinement in z direction
      case(23)
         write(*,*)'set number of anisotropic refinements'
         write(*,*)'-1 for infinite loop'
         read(*,*), iso_ans
         info=0
!     ...refine elements
         if(iso_ans.ne.-1) then
!        ...ref_xyz = 2 for anisotropic refinements in z-direction
            ref_xyz=2
            do iii=1,iso_ans
               call setAnisoRef(ref_xyz)
            enddo
!        ...update gdof and Ddof
            call update_gdof
            call update_Ddof
         else
!     ...ref_xyz = 2 for anisotropic refinements in z-direction
         ref_xyz=2
         do while(1.gt.0)
            call setAnisoRef(ref_xyz)
         enddo
!     ...update gdof and Ddof
         call update_gdof
         call update_Ddof
      endif
!
!-----------------------------------------------------------------------
!                            SOLVERS
!-----------------------------------------------------------------------
!
!  ...frontal solve
      case(40)
         call select_phys_problem(NO_PROBLEM)
         call set_physAm(NO_PROBLEM, physNick,flag)
         call solve1(MY_NR_RHS)
!
!  ...MUMPS solve
      case(41)
         write(*,*) 'Select:  1 = Linear  Problem'
         write(*,*) '         2 = Lasing  Problem'
         write(*,*) '         3 = Coupled Problem'
         read(*,*) chooseProbType
         select case(chooseProbType)
            case(1)
               if(NONLINEAR_FLAG.eq.1) then
                  write(*,*) 'NONLINEAR_FLAG must be 0 for linear problem'
                  stop 1
               endif
               call select_phys_problem(NO_PROBLEM)
               call set_physAm(NO_PROBLEM, physNick,flag)
!
               IPRINT_TIME = 1
               call mumps_sc('H')
               IPRINT_TIME = 0
               solved = .true.
!
            case(2)
               if(NONLINEAR_FLAG.ne.1) then
                  write(*,*) 'NONLINEAR_FLAG must be 1 for nonlinear problem. stop.'
                  stop
               endif
!           ...set components
               No1 = 1; No2 = 2
               if(NEXACT.ne.0) then
                  write(*,*) 'Error in main: NEXACT must be 0 for Lasing problem'
                  stop
               endif
               write(*,*) 'Enter tolerance for Lasing problem stopping criterion'
               read(*,*) stopEpsilon
               physNick = 1
               L2NormDiff = 1.d0
               FieldNormQ = 1.d0
               do
                  if((L2NormDiff/FieldNormQ).lt.(stopEpsilon)) exit
!
!              ...solve for signal first
                  PHYSAm(1:6) = (/.false.,.true.,.false.,.false.,.true.,.false./)
                  NO_PROBLEM = 3
                  call update_Ddof
                  IPRINT_TIME = 1
                  call mumps_sc('H')
                  IPRINT_TIME = 0
!
!              ...next solve for pump
                  PHYSAm(1:6) = (/.false.,.false.,.true.,.false.,.false.,.true./)
                  NO_PROBLEM = 4
                  call update_Ddof
!
                  IPRINT_TIME = 1
                  call mumps_sc('H')
                  IPRINT_TIME = 0
!
!              ...copy components and calculate norm corresponding to signal
                  flag = 0; flag(5)= 1
                  NO_PROBLEM = 3
                  call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
                  call get_Norm(flag,No2, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
                  write(*,*) 'from main: L2NormDiff = ', L2NormDiff
                  write(*,*) 'from main: FieldNormQ = ', FieldNormQ
                  call copy_coms(No1,No2)
!           ...end nonlinear loop
               end do
               solved = .true.
!           ...calculating power
               write(*,*) 'Computing power..'
               call get_power
!
            case(3)
               write(*,*) 'MUMPS: case 3 not configured yet'
!
            case default
               write(*,*) 'Invalid problem selected. stop.'
               stop
         end select
!
!  ...PARDISO solve
      case(42)
!     ...select problem
         write(*,*) 'Select:  1 = Linear  Problem'
         write(*,*) '         2 = Lasing  Problem'
         write(*,*) '         3 = Coupled Problem'
         read(*,*) chooseProbType
         select case(chooseProbType)
!        ...linear problem
            case(1)
               if(NONLINEAR_FLAG.eq.1) then
                  write(*,*) 'NONLINEAR_FLAG must be 0 for linear problem. stop.'
                  stop
               endif
               call select_phys_problem(NO_PROBLEM)
               call set_physAm(NO_PROBLEM, physNick,flag)
!           ...solve with pardiso solver
               IPRINT_TIME = 1
               call pardiso_sc('H')
               IPRINT_TIME = 0
               solved = .true.
!        ...lasing problem (nonlinear, no heat)
            case(2)
               if(NONLINEAR_FLAG.ne.1) then
                  write(*,*) 'NONLINEAR_FLAG must be 1 for nonlinear problem. stop.'
                  stop
               endif
!           ...set components
               No1 = 1; No2 = 2
               if(NEXACT.ne.0) then
                  write(*,*) 'NEXACT must be 0 for Lasing problem. stop.'
                  stop
               endif
               write(*,4200) 'Enter tolerance for Lasing problem stopping criterion (e.g., 0.001d0)'
               read(*,*) stopEpsilon
               physNick = 1
               L2NormDiff = 1.d0
               FieldNormQ = 1.d0
!           ...do until stopping criterion is satisfied
               write(*,4200) ' Beginning nonlinear iterations..'
           4200 format(/,A)
           4201 format(A,/)
               i = 0
               do
!              ...check stopping criterion
                  if((L2NormDiff/FieldNormQ).lt.(stopEpsilon)) then
                     write(*,*) ' Stopping criterion satisfied.'
                     write(*,4210) ' Ending nonlinear loop after ', i, 'iterations.'
                 4210 format(A,I3,A)
                     exit
                  else
                     write(*,4220) ' Stopping criterion not yet satisfied. i = ', i
                 4220 format(A,I3)
                     write(*,4230) '   ', L2NormDiff/FieldNormQ, ' > ', stopEpsilon
                 4230 format(A, F10.4,A,F10.4)
                     write(*,*) '---------------------------------------------'
                     write(*,4201) ' Proceed with nonlinear iterations..'
                  endif
!              ...solve for signal first
                  write(*,*) '   Signal solve..'
                  NO_PROBLEM = 3
                  PHYSAm(1:6) = (/.false.,.true.,.false.,.false.,.true.,.false./)
                  call update_Ddof
                  IPRINT_TIME = 1
                  call pardiso_sc('H')
                  IPRINT_TIME = 0
!
!              ...next solve for pump
                  write(*,*) '   Pump solve..'
                  NO_PROBLEM = 4
                  PHYSAm(1:6) = (/.false.,.false.,.true.,.false.,.false.,.true./)
                  call update_Ddof
                  IPRINT_TIME = 1
                  call pardiso_sc('H')
                  IPRINT_TIME = 0
!
!              ...copy components and calculate norm corresponding to signal
                  flag = 0; flag(5) = 1
                  NO_PROBLEM = 3
                  call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
                  call get_Norm(flag,No2, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
                  write(*,4240) '   L2NormDiff = ', L2NormDiff
                  write(*,4240) '   FieldNormQ = ', FieldNormQ
              4240 format(A,F10.4)
                  write(*,*)
                  call copy_coms(No1,No2)
                  i = i + 1
!           ...end nonlinear loop
               end do
               solved = .true.
!           ...calculating power
               write(*,*) 'Computing power..'
               call get_power
!        ...coupled problem (nonlinear Maxwell weakly coupled with heat equation)
            case(3)
               if((LASER_MODE.ne.1).or.(NONLINEAR_FLAG.ne.1).or.(NEXACT.ne.0)) then
                  write(*,*) 'NONLINEAR_FLAG must be 1, LASER_MODE must be 1, and ..'
                  write(*,*) ' ..NEXACT must be 0 for coupled nonlinear problem. stop.'
                  stop
               endif
               write(*,*) 'Enter tolerance for Lasing problem stopping criterion (e.g., 0.001d0)'
               read(*,*) stopEpsilon
!           ...start time stepping
               write(*,*) 'Begin time stepping..'
               do nonlin_steps = 1,NSTEPS
                  write(*,*) '  time step = ', nonlin_steps
!              ...first solve Nonlinear Maxwell loop for SIGNAL and PUMP
!              ...set components
                  No1 = 1; No2 = 2
                  physNick = 1
                  L2NormDiff = 1.d0
                  FieldNormQ = 1.d0
                  do
                     if((L2NormDiff/FieldNormQ).lt.(stopEpsilon)) exit
!                 ...solve for signal first
                     PHYSAm(1:6) = (/.false.,.true.,.false.,.false.,.true.,.false./)
                     NO_PROBLEM = 3
                     call update_Ddof
!
                     IPRINT_TIME = 1
                     call pardiso_sc('H')
                     IPRINT_TIME = 0
!
!                 ...next solve for pump
                     PHYSAm(1:6) = (/.false.,.false.,.true.,.false.,.false.,.true./)
                     NO_PROBLEM = 4
                     call update_Ddof
!
                     IPRINT_TIME = 1
                     call pardiso_sc('H')
                     IPRINT_TIME = 0
!
!                 ...copy components and calculate norm corresponding to signal
                     flag = 0; flag(5)= 1
                     NO_PROBLEM = 3
                     call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
                     call get_Norm(flag,No2, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
                     write(*,*) 'from main: L2NormDiff = ', L2NormDiff
                     write(*,*) 'from main: FieldNormQ = ', FieldNormQ
                     call copy_coms(No1,No2)
                  end do
                  call get_power
!              ...now solve heat equation
                  PHYSAm(1:6) = (/.true.,.false.,.false.,.true.,.false.,.false./)
                  NO_PROBLEM = 2
                  call update_Ddof
!
                  IPRINT_TIME = 1
                  call pardiso_sc('H')
                  IPRINT_TIME = 0
!
!           ...end loop for time stepping
               enddo
               solved = .true.
!
            case default
               write(*,*) 'Invalid problem selected. stop.'
               stop
         end select
!
      case(50)
         call select_phys_problem(NO_PROBLEM)
         select case(NO_PROBLEM)
            case(1,2)
               flag= 0; flag(1)=1
               PHYSAm(1:6) = (/.true.,.false.,.false.,.true.,.false.,.false./)
            case(3)
               flag= 0; flag(5)=1
               PHYSAm(1:6) = (/.false.,.true.,.false.,.false.,.true.,.false./)
            case(4)
               flag= 0; flag(6)=1
               PHYSAm(1:6) = (/.false.,.false.,.true.,.false.,.false.,.true./)
         end select
         call compute_error(flag,1)
!
      case(51)
         call select_phys_problem(NO_PROBLEM)
         read(*,*) numProb
         NO_PROBLEM=numProb
         select case(NO_PROBLEM)
            case(1,2)
               flag= 0; flag(1)=1
               PHYSAm(1:6) = (/.true.,.false.,.false.,.true.,.false.,.false./)
            case(3)
               flag= 0; flag(5)=1
               PHYSAm(1:6) = (/.false.,.true.,.false.,.false.,.true.,.false./)
            case(4)
               flag= 0; flag(6)=1
               PHYSAm(1:6) = (/.false.,.false.,.true.,.false.,.false.,.true./)
         end select
         No = 1
         call get_Norm(flag,No, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
         write(*,*) 'The Field Norms are :'
         write(*,*) 'FieldNormH = ', FieldNormH
         write(*,*) 'FieldNormE = ', FieldNormE
         write(*,*) 'FieldNormV = ', FieldNormV
         write(*,*) 'FieldNormQ = ', FieldNormQ
!
!-----------------------------------------------------------------------
!                            TESTS
!-----------------------------------------------------------------------
!
      case(60)
         call test_Bessel
      case(61)
         call get_power
      case(62)
         call my_sizetest
      case(63)
         No1 = 1; No2 = 2
         call copy_coms(No1,No2)
         flag = 0; flag(5)= 1
         !NO_PROBLEM = 3
         L2NormDiff = 1.d0
         call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
         write(*,*) 'from test copy coms: L2NormDiff = ', L2NormDiff
!
!..end select menu options
   end select
!
!..go back to menu
   goto 10
98 call finalize
!
99 call MPI_FINALIZE ( ierr )
!..END MPI
!
   stop
!
end program main
!
!
!-----------------------------------------------------------------------
!                   subroutine select_phys_problem
!-----------------------------------------------------------------------
!
subroutine select_phys_problem(NO_PROBLEM)
!
   implicit none
!
   integer, intent(out) :: NO_PROBLEM
   integer              :: numProb
!
   write(*,*) 'Select problem:  1 = Heat Step'
   write(*,*) '                 2 = Multi-step Heat equation'
   write(*,*) '                 3 = Time harmonic Maxwell (signal)'
   write(*,*) '                 4 = Time harmonic Maxwell (pump)'
   read(*,*) numProb
   NO_PROBLEM = numProb
!
end subroutine select_phys_problem
!
!
!-----------------------------------------------------------------------
!                       subroutine set_physAm
!-----------------------------------------------------------------------
!
subroutine set_physAm(NO_PROBLEM, PhysNick,Flag)
!
   use physics, only: PHYSAm
!
   implicit none
!
   integer,               intent(in)  :: NO_PROBLEM
   integer,               intent(out) :: PhysNick
   integer, dimension(6), intent(out) :: Flag
!
   select case(NO_PROBLEM)
      case(1,2)
         PhysNick = 1000
         Flag=0; Flag(1)=1
         PHYSAm(1:6) = (/.true.,.false.,.false.,.true.,.false.,.false./)
      case(3)
         PhysNick = 1
         Flag=0; Flag(5)=1
         PHYSAm(1:6) = (/.false.,.true.,.false.,.false.,.true.,.false./)
      case(4)
         PhysNick = 1
         Flag=0; Flag(6)=1
         PHYSAm(1:6) = (/.false.,.false.,.true.,.false.,.false.,.true./)
      case default
         write(*,*) 'invalid NO_PROBLEM param. stop.'
         stop
   end select
!
end subroutine set_physAm












