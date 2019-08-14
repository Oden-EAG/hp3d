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
   use commonParam
   use laserParam
   use assembly
   use assembly_sc, only: IPRINT_TIME
   use stc,         only: STORE_STC, HERM_STC
!
   implicit none
!
!----------------------------------------------------------------------
!
!..declare variables
   real*8  :: err,rnorm,rvoid,factor,geo_err,geo_rnorm
   integer :: mdle,i,j,kref,idec,nreflag,istep,nstop
   integer :: numRef,iref,numLaserIter,nord,idec_solve,numProb
   integer :: ic,iel,iso_ans,iii,ref_xyz,nuniform,naniso
   integer :: physNick,chooseProbType,fld,numPts
   logical :: ires,solved
   real*8  :: FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!
   integer              :: flag(6)
   integer, allocatable :: list_elem(:)
!
!..variables for nonlinear loop
   integer :: No,No1,No2,time_step
   real*8  :: L2NormDiff,stopEpsilon,scale
   real*8  :: L2NormDiffIter(100)
!
!..Paraview info array
   integer :: iParAttr(6)
   character(len=2) :: vis_level
!
!..timer variables
   integer*8 :: t1,t2,clock_rate,clock_max
!
!..OpenMP variables
   integer :: num_threads, omp_get_num_threads
!
   character :: arg
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
!..display menu
   idec = 1
   solved = .false.
   ires = .true.
!
   flag(1:6) = 0
!
!..PML
   call set_PML
!
!..print header
   write(*,*)'                   '
   write(*,*)'                   '
   write(*,*)'// ------ THE LASER PROBLEM ------ //'
   write(*,*)'                   '
   write(*,*)'                   '
!
!..initialize problem
   call initialize
!
!..print problem parameters
   if (GEOM_NO .eq. 5) then
      write(*,9000) ' Fiber length             = ', ZL
      write(*,9000) ' Signal frequency         = ', OMEGA_SIGNAL
      write(*,9000) ' Wavelengths/Unit length  = ', 1.d0/(LAMBDA_SIGNAL/REF_INDEX_CORE)
      write(*,9000) ' Numerical Aperture       = ', NA
      write(*,9000) ' V-number                 = ', VNUM
   endif
   if (ANISO_REF_INDEX .eq. 1) then
      write(*,9010) ' ANISO_REF_INDEX          = ', ANISO_REF_INDEX
      write(*,9000) ' CORE_NY                  = ', CORE_NY
      write(*,9000) ' CLAD_NY                  = ', CLAD_NY
   endif
   if (USE_PML) then
      write(*,9000) ' PML_REGION               = ', PML_REGION
   endif
   if (NONLINEAR_FLAG .eq. 1) then
      write(*,9020) ' Raman gain               = ', RAMAN_GAIN
      write(*,9020) ' Active gain              = ', ACTIVE_GAIN
   endif
   write(*,9030) ' Polynomial order (x,y,z) = ', ORDER_APPROX_X,ORDER_APPROX_Y,ORDER_APPROX_Z
   write(*,9010) ' NEXACT                   = ', NEXACT
   write(*,9010) ' FAST INTEGRATION         = ', FAST_INT
   if (HEAT_FLAG .eq. 1) then
      write(*,9010) ' NSTEPS                   = ', NSTEPS
      write(*,9000) ' DELTA_T                  = ', DELTA_T
   endif
   if (ANISO_HEAT .eq. 1) then
      write(*,9020) ' ALPHA_Z                  = ', ALPHA_Z
   endif
 9000 format(A,F10.6)
 9010 format(A,I3)
 9020 format(A,ES10.2)
 9030 format(A,' (',I1,',',I1,',',I1,') ')
!
!..check if HEAT_FLAG and NONLINEAR_FLAG are compatible:
   if((HEAT_FLAG.eq.1).and.(NONLINEAR_FLAG.ne.1)) then
      write(*,*) 'error from main: ', &
         'NONLINEAR_FLAG must be 1 when running HEAT_FLAG=1. stop.'
      stop
   endif
   write(*,*)
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
   write(*,*) 'FLAGS:'
   write(*,*) ' ISTC_FLAG: ', ISTC_FLAG
   write(*,*) ' STORE_STC: ', STORE_STC
   write(*,*) ' HERM_STC : ', HERM_STC
   write(*,*) ' PHYSAi   : ', PHYSAi
   write(*,*)
!
   IPRINT_TIME = 1
!
#if DEBUG_MODE
   write(*,*) '========================='
   write(*,*) '  RUNNING in DEBUG_MODE  '
   write(*,*) '========================='
#endif
!
!..display menu in infinite loop
 10 continue
!
!     -- INTERACTIVE MODE --
!
   write(*,*)
   write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   write(*,*) 'Quit ...................................0'
   write(*,*) '                                         '
   write(*,*) '     ---------     I/O    ----------     '
   write(*,*) 'Geometry graphics (X11) ................1'
   write(*,*) 'HP3D graphics (X11) ....................2'
   write(*,*) 'Paraview ...............................3'
   write(*,*) 'Print Data Structure arrays ............5'
   write(*,*) 'DumpOut Data Structure .................6'
   write(*,*) '                                         '
   write(*,*) '    ----------  Geometry  ----------     '
   write(*,*) 'Geometry error ........................10'
   write(*,*) '                                         '
   write(*,*) '            --  Refinements  --          '
   write(*,*) 'Single Uniform  p-refinement...........15'
   write(*,*) 'Single Uniform  h-refinement ..........16'
   write(*,*) 'Single Adaptive h-refinement ..........17'
   write(*,*) 'Multi-step Uniform  h-refinements .....18'
   write(*,*) 'Multi-step Adaptive hp-refinements ....19'
   write(*,*) 'Interactive h-refinement ..............20'
   write(*,*) 'M uniform, N z-refs, and solve ........21'
   write(*,*) 'M xz-refs, N z-refs, and solve ........22'
   write(*,*) '0 uniform, N anisotropic refinements ..23'
   write(*,*) 'Fiber core refinements ................24'
   write(*,*) 'Fiber input refinements ...............25'
   write(*,*) '                                         '
   write(*,*) '      --------  Solvers  --------        '
   write(*,*) 'Solve (Frontal) .......................40'
   write(*,*) 'Solve (MUMPS)   .......................41'
   write(*,*) 'Solve (PARDISO) .......................42'
   write(*,*) '                                         '
   write(*,*) '      ---------  Error  ---------        '
   write(*,*) 'Compute Field Error  ..................50'
   write(*,*) 'Compute norm of current solution.......51'
   write(*,*) '                                         '
   write(*,*) '      ---------  Misc   ---------        '
   write(*,*) 'Compute Power .........................61'
   write(*,*) 'Compute Temperature ...................62'
   write(*,*) 'Test Copy Coms ........................63'
   write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   read (*,*) idec
!
!----------------------------------------------------------------------
!
   select case(idec)
!  ...quit
      case(0) ; call finalize ; stop
!
!-----------------------------------------------------------------------
!                              GRAPHICS
!-----------------------------------------------------------------------
!
!  ...GMP graphics
      case(1) ; call graphg
!
!  ...hp3d graphics
      case(2) ; call graphb
!
!  ...Paraview graphics
!  ...iParAttr parameter specifies which fields we dump out to paraview:
!     - 1 H1 (heat)
!     - 2 H(curl) (signal, E and H flux)
!     - 2 H(curl) (pump,   E and H flux)
!     - 1 H(div) (heat flux)
!     - 6 L2 (signal, E and H field)
!     - 6 L2 (pump,   E and H Field)
!
      case(3)
         iParAttr = (/1,2,2,1,6,6/)
         write(*,300) ' paraview output: select fields...'   , &
                      '  - 1 H1 (heat)'                      , &
                      '  - 2 H(curl) (signal, E and H flux)' , &
                      '  - 2 H(curl) (pump,   E and H flux)' , &
                      '  - 1 H(div) (heat flux)'             , &
                      '  - 6 L2 (signal, E and H field)'     , &
                      '  - 6 L2 (pump,   E and H Field)'
         read (*,*) iParAttr(1),iParAttr(2),iParAttr(3), &
                    iParAttr(4),iParAttr(5),iParAttr(6)
     300 format(A,/,A,/,A,/,A,/,A,/,A,/,A)
!
         write(*,*) 'paraview output: select VLEVEL (0-4)...'
         read (*,*) vis_level
         select case(vis_level)
            case('0','1','2','3','4')
               VLEVEL = vis_level
            case default
               write(*,*) ' invalid VLEVEL. setting VLEVEL=3 (default).'
               VLEVEL = '3'
         end select
!
         call my_paraview_driver(iParAttr)
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
!
!-----------------------------------------------------------------------
!                            GEOMETRY ERROR
!-----------------------------------------------------------------------
!
!  ...geometry error
      case(10) ; call geom_error(err,rnorm)
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
         !if (.not.solved) write(*,*) 'You have not solved.'
         !if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
         !call select_phys_problem(NO_PROBLEM)
         !call set_physAm(NO_PROBLEM, physNick,flag)
         !call refine_DPG(IUNIFORM,1,0.25d0,flag,physNick,ires,0, nstop)
         call global_href
         call close_mesh
         call update_gdof
         call update_Ddof
         solved=.false.
!
!  ...Single adaptive refinement
      case(17)
         if (.not.solved) then
            write(*,*) 'You have not solved. Cannot adaptively refine.'
         else
            call select_phys_problem(NO_PROBLEM)
            call set_physAm(NO_PROBLEM, physNick,flag)
            nreflag=0
            do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
               write(*,*) 'REFINEMENT FLAG:h-refine=1,p-refine=2,hp-refine=3'
               write(*,*) '0.d0<FACTOR<1.d0'
               write(*,*) 'Provide: REFINEMENT FLAG, FACTOR, NO_PROBLEM'
               read(*,*) nreflag,factor
            enddo
            call refine_DPG(IADAPTIVE,nreflag,factor,flag,physNick,ires,0, nstop)
            if (nstop.eq.1) write(*,*) 'No elements were refined.'
            if (nstop.eq.0) solved=.false.
         endif
!
!     ...Multi-step uniform h refinement
         case(18)
!        ...select number of refinements and the solver
            write(*,*) 'NUMBER OF REFINEMENTS > 0'
            write(*,*) 'SOLVER:  1 = Frontal'
            write(*,*) '         2 = OMP MUMPS'
            write(*,*) '         3 = PARDISO'
            write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
            read(*,*) nuniform,idec_solve
!        ...select problem
            call select_phys_problem(NO_PROBLEM)
            call set_physAm(NO_PROBLEM, physNick,flag)
!        ...begin refinements
            do i=0,nuniform
!           ...solve first if needed
               if (.not.solved) then
                  select case(idec_solve)
                     case(1)
                        call solve1(MY_NR_RHS)
                     case(2)
                        call mumps_sc(arg)
                     case(3)
                        call pardiso_sc(arg)
                  end select
               endif
!           ...say it has solved and save results to paraview file
               solved=.true.
!           ...display error and refine if necessary
               if (i.ne.nuniform) then
                  call refine_DPG(IUNIFORM,1,0.25d0,flag,physNick,ires,0, nstop)
                  if (nstop.eq.1) then
                     write(*,*) 'No elements were refined.'
                     write(*,1800) i
 1800                format('Exiting loop after ',i2,' refinements...')
                     cycle
                  else
                     solved=.false.
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
               write(*,*) 'SOLVER:  1 = Frontal'
               write(*,*) '         2 = OMP MUMPS'
               write(*,*) '         3 = PARDISO'
               write(*,*) 'Provide: NUMBER OF REFINEMENTS, ', &
                          'REFINEMENT FLAG, FACTOR, SOLVER'
               read(*,*) nuniform,nreflag,factor,idec_solve
!
               call select_phys_problem(NO_PROBLEM)
               call set_physAm(NO_PROBLEM, physNick,flag)
            enddo
!
            do i=0,nuniform
!           ...solve first if needed
               if (.not.solved) then
                  select case(idec_solve)
                     case(1)
                        call solve1(MY_NR_RHS)
                     case(2)
                        call mumps_sc(arg)
                     case(3)
                        call pardiso_sc(arg)
                     case default
                        write(*,*) 'invalid SOLVER param. stop.'
                        stop
                  end select

               endif
!           ...say it has solved
               solved=.true.
!
!           ...then display error and refine if necessary
               if (i.ne.nuniform) then
                  call refine_DPG(IADAPTIVE,nreflag,factor,flag,physNick,ires,0, nstop)
                  if (nstop.eq.1) then
                     write(*,*) 'No elements were refined.'
                     write(*,2100) i
                  cycle
               else
                  solved=.false.
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
            end select
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
         write(*,*) 'M uniform, N anisotropic refinements.'
         write(*,*) 'Solver: Frontal=1, OMP MUMPS=2, PARDISO=3.'
         write(*,*) 'Provide: Uniform refinements, Anisotropic refinements, Solver.'
         read (*,*) nuniform,naniso,idec_solve
         call select_phys_problem(NO_PROBLEM)
         call set_physAm(NO_PROBLEM, physNick,flag)
!
!         if (naniso .gt. 0) then
!            write(*,*) 'Set direction of anisotropic refinements:'
!            write(*,*) '  x = 1, y = 2, z = 3, xy = 4, xz = 5, yz = 6'
!            read(*,*), ref_xyz
!            call set_ANISO_FLAG(ref_xyz)
!         endif
!
         if (nuniform .lt. 0 .or. naniso .lt. 0) then
            write(*,*) 'main: invalid param nuniform or naniso. stop.'
            stop
         endif
!
         do i=0,nuniform
!        ...solve first if needed
            if (.not.solved) then
               select case(idec_solve)
                  case(1)
                     call solve1(MY_NR_RHS)
                  case(2)
                     call mumps_sc(arg)
                  case(3)
                     call pardiso_sc(arg)
               end select
            endif
!        ...say it has solved
            solved=.true.
!        ...display error and refine if necessary
            if (i.lt.nuniform) then
               call refine_DPG(IUNIFORM,1,0.25d0,flag,physNick,ires,0, nstop)
               if (nstop.eq.1) then
                  write(*,*) 'No elements were refined.'
                  write(*,2100) i
 2100             format('Exiting loop after ',i2,' refinements...')
                  cycle
               else
                  solved=.false.
               endif
            else ! Last step only anisotropic refinements
               write(*,*) 'Finished uniform refinements; proceeding with anisotropic refinements.'
               ANISO_FLAG = 1 ! z
               do iii=1,naniso
                  call refine_DPG(IANISOTROPIC,1,0.25d0,flag,physNick,ires,1, nstop)
                  select case(idec_solve)
                     case(1)
                        call solve1(MY_NR_RHS)
                     case(2)
                        call mumps_sc(arg)
                     case(3)
                        call pardiso_sc(arg)
                  end select
               enddo
               call refine_DPG(INOREFINEMENT,nreflag,0.25d0,flag,physNick,ires,0, nstop)
            endif
         enddo
!
!  ...anisotropic h-refinement in x,y, or z-direction
      case(23)
         !write(*,*) 'Set direction of anisotropic refinements:'
         !write(*,*) '  x = 1, y = 2, z = 3, xy = 4, xz = 5, yz = 6.'
         !read(*,*), ref_xyz
         ref_xyz = 3
         call set_ANISO_FLAG(ref_xyz)
!
         write(*,*)'Set number of anisotropic refinements.'
         read(*,*), iso_ans
!     ...refine elements
         if(iso_ans.lt.1) then
            write(*,*) 'main: invalid iso_ans param. stop.'
            stop
         endif
!     ...do refinements
         do iii=1,iso_ans
            call setAnisoRef(ANISO_FLAG)
         enddo
!     ...update gdof and Ddof
         call close_mesh
         call update_gdof
         call update_Ddof
!
!  ...fiber core refinement
      case(24)
         call ref_core()
!
!     ...recover 1-irregular mesh, update geometry and Dirichlet dof's
         call close_mesh
         call update_gdof
         call update_ddof
!
!  ...fiber input refinement
      case(25)
         kref = 110 !xy
         call ref_layer(kref)
!
!     ...recover 1-irregular mesh, update geometry and Dirichlet dof's
         call close_mesh
         call update_gdof
         call update_ddof
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
               if(NONLINEAR_FLAG.ne.0) then
                  write(*,*) 'NONLINEAR_FLAG must be 0 for linear problem. stop.'
                  stop
               endif
               call select_phys_problem(NO_PROBLEM)
               call set_physAm(NO_PROBLEM, physNick,flag)
!
               call mumps_sc(arg)
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
                  write(*,*) 'Error in main: NEXACT must be 0 for Lasing problem. stop.'
                  stop
               endif
               write(*,*) 'Enter tolerance for Lasing problem stopping criterion.'
               read(*,*) stopEpsilon
               physNick = 1
               L2NormDiff = 1.d0
               FieldNormQ = 1.d0
               do
                  if((L2NormDiff/FieldNormQ).lt.(stopEpsilon)) exit
!
!              ...solve for signal first
                  NO_PROBLEM = 3
                  call set_physAm(NO_PROBLEM, physNick,flag)
                  call update_Ddof
                  call mumps_sc(arg)
!
!              ...next solve for pump
                  NO_PROBLEM = 4
                  call set_physAm(NO_PROBLEM, physNick,flag)
                  call update_Ddof
                  call mumps_sc(arg)
!
!              ...copy components and calculate norm corresponding to signal
                  NO_PROBLEM = 3
                  call set_physAm(NO_PROBLEM, physNick,flag)
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
               fld=2;numPts=0
               call get_power(fld,numPts,-1)
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
         write(*,*) 'Select:  1 = Linear    Maxwell Problem'
         write(*,*) '         2 = Nonlinear Maxwell Problem'
         write(*,*) '         3 = Coupled   Maxwell/Heat Problem'
         write(*,*) '         4 = Linear    Heat Problem'
         read(*,*) chooseProbType
         select case(chooseProbType)
!        ...linear problem
            case(1)
               if(NONLINEAR_FLAG.ne.0) then
                  write(*,*) 'NONLINEAR_FLAG must be 0 for linear problem. stop.'
                  stop
               endif
               call select_phys_problem_maxwell(NO_PROBLEM)
               if(NO_PROBLEM.ne.3 .and. NO_PROBLEM.ne.4) then
                  write(*,*) 'Linear Maxwell problem requires NO_PROBLEM=3,4. stop.'
                  stop
               endif
               call set_physAm(NO_PROBLEM, physNick,flag)
!           ...solve with pardiso solver
               call pardiso_sc(arg)
               call refine_DPG(INOREFINEMENT,nreflag,0.25d0,flag,physNick,ires,0, nstop)
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
                  write(*,*) 'NEXACT must be 0 for nonlinear Maxwell problem. stop.'
                  stop
               endif
               !write(*,4200) 'Enter tolerance for nonlinear Maxwell problem stopping criterion (e.g., 0.001d0)'
               !read(*,*) stopEpsilon
               stopEpsilon = 1.0d-4
               physNick = 1
               L2NormDiff = 1.d0
               FieldNormQ = 1.d0
!           ...do until stopping criterion is satisfied
               write(*,4200) ' Beginning nonlinear iterations..'
           4200 format(/,A)
           4201 format(A,/)
               L2NormDiffIter = 0.d0
               i = 0
               do
!              ...check stopping criterion
                  L2NormDiffIter(i+1) = L2NormDiff/FieldNormQ
                  if((L2NormDiff/FieldNormQ).lt.(stopEpsilon)) then
                     write(*,4230) '   ', L2NormDiff/FieldNormQ, ' < ', stopEpsilon
                     write(*,*) ' Stopping criterion satisfied.'
                     write(*,*) '---------------------------------------------'
                     write(*,4210) ' Ending nonlinear loop after ', i, ' iterations.'
                     write(*,*) '---------------------------------------------'
                 4210 format(A,I3,A)
                     exit
                  else
                     write(*,4220) ' Stopping criterion not yet satisfied. i = ', i
                 4220 format(A,I3)
                     write(*,4230) '   ', L2NormDiff/FieldNormQ, ' > ', stopEpsilon
                 4230 format(A, F11.5,A,F11.5)
                     write(*,*) '---------------------------------------------'
                     write(*,4201) ' Proceed with nonlinear iterations..'
                  endif
                  if (i .ge. MAX_ITER) then
                     write(*,4210) ' Ending nonlinear loop after', MAX_ITER,' iterations, no convergence.'
                     exit
                  endif
!              ...solve for signal first
                  write(*,*) '   Signal solve..'
                  NO_PROBLEM = 3
                  call set_physAm(NO_PROBLEM, physNick,flag)
                  call update_Ddof
                  call pardiso_sc(arg)
                  QUIET_MODE = .true.; IPRINT_TIME = 0
                  write(*,*) ''; write(*,*) '   Signal residual:'
                  call refine_DPG(INOREFINEMENT,1,0.25d0,flag,physNick,ires,0, nstop)
                  write(*,*); write(*,*)
                  QUIET_MODE = .false.; IPRINT_TIME = 1
!
!              ...next solve for pump
                  write(*,*) '   Pump solve..'
                  NO_PROBLEM = 4
                  call set_physAm(NO_PROBLEM, physNick,flag)
                  call update_Ddof
                  call pardiso_sc(arg)
                  QUIET_MODE = .true.; IPRINT_TIME = 0
                  write(*,*) ''; write(*,*) '   Pump residual:'
                  call refine_DPG(INOREFINEMENT,1,0.25d0,flag,physNick,ires,0, nstop)
                  QUIET_MODE = .false.; IPRINT_TIME = 1
!
!              ...copy components and calculate norm corresponding to signal
                  NO_PROBLEM = 3
                  call set_physAm(NO_PROBLEM, physNick,flag)
                  if (i .gt. 0) then
                     call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
                     call get_Norm(flag,No2, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
                  endif
                  write(*,4240) '   L2NormDiff = ', L2NormDiff
                  write(*,4240) '   FieldNormQ = ', FieldNormQ
              4240 format(A,F10.4)
                  write(*,*)
                  call copy_coms(No1,No2)
                  i = i + 1
!           ...end nonlinear loop
               end do
!
!!           ...Last step only display (no refinement)
!               QUIET_MODE = .true.; IPRINT_TIME = 0
!               write(*,*) ''; write(*,*) '   Pump residual:'
!               NO_PROBLEM = 4
!               call set_physAm(NO_PROBLEM, physNick,flag)
!               call refine_DPG(INOREFINEMENT,1,0.25d0,flag,physNick,ires,0, nstop)
!!
!               write(*,*) ''; write(*,*) '   Signal residual:'
!               NO_PROBLEM = 3
!               call set_physAm(NO_PROBLEM, physNick,flag)
!               call update_Ddof
!               call pardiso_sc(arg)
!               call refine_DPG(INOREFINEMENT,1,0.25d0,flag,physNick,ires,0, nstop)
!               QUIET_MODE = .false.; IPRINT_TIME = 1
!
               write(*,*) 'L2NormDiff/FieldNormQ:'
               do j=1,i+1
                  write(*,4241) L2NormDiffIter(j)
             4241 format(es14.5)
               enddo 
               solved = .true.
!
!        ...coupled problem (nonlinear Maxwell weakly coupled with heat equation)
            case(3)
               if((HEAT_FLAG.ne.1).or.(NONLINEAR_FLAG.ne.1).or.(NEXACT.ne.0)) then
                  write(*,*) 'NONLINEAR_FLAG must be 1, HEAT_FLAG must be 1, and NEXACT must be 0 for coupled nonlinear problem. stop.'
                  stop
               endif
               !write(*,4200) 'Enter tolerance for nonlinear Maxwell problem stopping criterion (e.g., 0.001d0)'
               !read(*,*) stopEpsilon
               stopEpsilon = 1.0d-4
!           ...set components
               No1 = 1; No2 = 2
!           ...start time stepping
               write(*,*) 'Begin time stepping..'
               do time_step = 1,NSTEPS
                  TIMESTEP = time_step !setting global variable
                  write(*,*) '---------------------------------------------'
                  write(*,4220) ' Proceeding with time step = ', time_step
                  !..activate to compute maxwell only in first time step
                  if(time_step .gt. 1 .and. time_step .le. 100) goto 420
!              ...first solve Nonlinear Maxwell loop for SIGNAL and PUMP
                  physNick = 1
                  L2NormDiff = 1.d0
                  FieldNormQ = 1.d0
!              ...do until stopping criterion is satisfied
                  write(*,4200) ' Beginning nonlinear iterations..'
                  i = 0
                  do
!                 ...check stopping criterion
                     if((L2NormDiff/FieldNormQ).lt.(stopEpsilon)) then
                        write(*,4230) '   ', L2NormDiff/FieldNormQ, ' < ', stopEpsilon
                        write(*,*) ' Stopping criterion satisfied.'
                        write(*,*) '---------------------------------------------'
                        write(*,4210) ' Ending nonlinear loop after ', i, ' iterations.'
                        write(*,*) '---------------------------------------------'
                        exit
                     else
                        write(*,4220) ' Stopping criterion not yet satisfied. i = ', i
                        write(*,4230) '   ', L2NormDiff/FieldNormQ, ' > ', stopEpsilon
                        write(*,*) '---------------------------------------------'
                        write(*,4201) ' Proceed with nonlinear iterations..'
                     endif
                     if (i .ge. MAX_ITER) then
                        write(*,4210) ' Ending nonlinear loop after', MAX_ITER,' iterations, no convergence.'
                     exit
                  endif
!
!                 ...next solve for pump (reordered for increasing pump power experiment)
                     write(*,*) '   Pump solve..'
                     NO_PROBLEM = 4
                     call set_physAm(NO_PROBLEM, physNick,flag)
                     call update_Ddof
                     call pardiso_sc(arg)
!
!                 ...solve for signal first
                     write(*,*) '   Signal solve..'
                     NO_PROBLEM = 3
                     call set_physAm(NO_PROBLEM, physNick,flag)
                     call update_Ddof
                     call pardiso_sc(arg)
!
!                 ...copy components and calculate norm corresponding to signal
                     NO_PROBLEM = 3
                     call set_physAm(NO_PROBLEM, physNick,flag)
                     if (i .gt. 0 .or. time_step .gt. 1) then
                        call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
                        call get_Norm(flag,No2, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
                     endif
                     write(*,4240) '   L2NormDiff = ', L2NormDiff
                     write(*,4240) '   FieldNormQ = ', FieldNormQ
                     call copy_coms(No1,No2)
                     i = i + 1
                  end do
!              ...calculating power
                  !if (time_step .eq. 1) then
                     write(*,*) 'Computing power..'
                     numPts=64
                     call get_power(2,numPts,time_step-1)
                     write(*,*) ''
                  !endif
                  !write(*,*) 'pausing before continuing with the heat solve'
                  !pause
 420              continue
!              ...now solve heat equation
                  write(*,*) ''
                  write(*,4210) ' Computing time-step ', time_step, ' for heat equation..'
                  NO_PROBLEM = 2
                  call set_physAm(NO_PROBLEM, physNick,flag)
                  call update_Ddof
                  call pardiso_sc(arg)
!              ...calculating temperature
                  call get_avgTemp(numPts,time_step-1)
!
!              ...write paraview output files
                  ! if(time_step .gt. 1 .and. time_step .le. 100) then
                  !    iParAttr = (/1,0,0,0,0,0/)
                  ! else
                  !    iParAttr = (/1,0,0,0,6,6/)
                  ! endif
                  iParAttr = (/1,0,0,0,0,0/)
                  call system_clock( t1, clock_rate, clock_max )
                  call my_paraview_driver(iParAttr)
                  if (IPRINT_TIME .eq. 1) then
                    call system_clock( t2, clock_rate, clock_max )
                    write(*,1420) real(t2 - t1,8)/real(clock_rate,8)
 1420                format(' Finished writing paraview output: ',f12.5,'  seconds',/)
                  endif
!
!           ...end loop for time stepping
               enddo
               solved = .true.
!
!        ...linear Heat problem
            case(4)
               if(NONLINEAR_FLAG.ne.0) then
                  write(*,*) 'NONLINEAR_FLAG must be 0 for linear problem. stop.'
                  stop
               endif
               call select_phys_problem_heat(NO_PROBLEM)
               if(NO_PROBLEM.ne.1 .and. NO_PROBLEM.ne.2) then
                  write(*,*) 'Linear Heat problem requires NO_PROBLEM=1,2. stop.'
                  stop
               endif
               call set_physAm(NO_PROBLEM, physNick,flag)
!
               if (NO_PROBLEM.eq.1) then
!              ...solve with pardiso solver
                  call pardiso_sc(arg)
                  call refine_DPG(INOREFINEMENT,nreflag,0.25d0,flag,physNick,ires,0, nstop)
                  solved = .true.
               else
!              ...start time stepping
                  write(*,*) 'Warning: Time stepping for linear heat problem with manufactured solution.'
                  write(*,*) 'Begin time stepping..'
                  do time_step = 1,NSTEPS
                     write(*,*) '  time step = ', time_step
!                 ...solve heat equation
                     call set_physAm(NO_PROBLEM, physNick,flag)
                     call update_Ddof
                     call pardiso_sc(arg)
                     solved = .true.
!              ...end loop for time stepping
                  enddo
               endif
!
            case default
               write(*,*) 'Invalid problem selected. stop.'
               stop
         end select
!
      case(50)
         call select_phys_problem_heat(NO_PROBLEM)
         call set_physAm(NO_PROBLEM, physNick,flag)
         call compute_error(flag,1)
!
      case(51)
         call select_phys_problem(NO_PROBLEM)
         call set_physAm(NO_PROBLEM, physNick,flag)
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
      case(61)
         write(*,*) 'Computing power..'
         write(*,*) 'Select field:  0 = pump'
         write(*,*) '               1 = signal'
         write(*,*) '               2 = signal and pump'
         read(*,*) fld
         write(*,*) 'Choose number sample points: 0 = default'
         read(*,*) numPts
         call get_power(fld,numPts,-1)
      case(62)
         write(*,*) 'Computing temperature..'
         write(*,*) 'Choose number sample points: 0 = default'
         read(*,*) numPts
         call get_avgTemp(numPts,-1)
      case(63)
         No1 = 1; No2 = 2
         call copy_coms(No1,No2)
         NO_PROBLEM = 3
         call set_physAm(NO_PROBLEM, physNick,flag)
         L2NormDiff = 1.d0
         call get_L2NormCOMS(flag,No1,No2, L2NormDiff)
         write(*,*) 'from test copy coms: L2NormDiff = ', L2NormDiff
!
!..end select menu options
   end select
!
!..go back to menu
   goto 10
   call finalize
!
!
end program main
!
!
!-----------------------------------------------------------------------
!                   subroutine select_phys_problem
!-----------------------------------------------------------------------
!
subroutine set_ANISO_FLAG(ref_xyz)
!
   use commonParam
!
   implicit none
!
   integer, intent(in) :: ref_xyz
!
   select case(ref_xyz)
      case(1); ANISO_FLAG = IREFINE_X
      case(2); ANISO_FLAG = IREFINE_Y
      case(3); ANISO_FLAG = IREFINE_Z
      case(4); ANISO_FLAG = IREFINE_XY
      case(5); ANISO_FLAG = IREFINE_XZ
      case(6); ANISO_FLAG = IREFINE_YZ
      case default
         write(*,*) 'set_ANISO_FLAG: invalid ref_xyz param. stop.'
         stop
   end select
!
end subroutine set_ANISO_FLAG
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
!                   subroutine select_phys_problem_maxwell
!-----------------------------------------------------------------------
!
subroutine select_phys_problem_maxwell(NO_PROBLEM)
!
   implicit none
!
   integer, intent(out) :: NO_PROBLEM
   integer              :: numProb
!
   write(*,*) 'Select problem:  3 = Time harmonic Maxwell (signal)'
   write(*,*) '                 4 = Time harmonic Maxwell (pump)'
   read(*,*) numProb
   NO_PROBLEM = numProb
!
end subroutine select_phys_problem_maxwell
!
!
!-----------------------------------------------------------------------
!                   subroutine select_phys_problem_heat
!-----------------------------------------------------------------------
!
subroutine select_phys_problem_heat(NO_PROBLEM)
!
   implicit none
!
   integer, intent(out) :: NO_PROBLEM
   integer              :: numProb
!
   write(*,*) 'Select problem:  1 = Heat Step'
   write(*,*) '                 2 = Multi-step Heat equation'
   read(*,*) numProb
   NO_PROBLEM = numProb
!
end subroutine select_phys_problem_heat
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
      write(*,*) 'set_physAm: invalid NO_PROBLEM param. stop.'
         stop
   end select
!
end subroutine set_physAm
!
!
!-------------------------
! subroutine ref_core
! refine fiber core in xy
!-------------------------
subroutine ref_core()
   use data_structure3D
   use assembly
!
   implicit none
!
   integer, parameter :: kref_mdlb = 110
   integer, parameter :: kref_mdlp = 10
!
   integer :: mdle, iel, ndom, nr_elem_ref
!
!..geometry dof (work space for nodcor)
   integer :: n_elem(NRELES)
!
   mdle=0; nr_elem_ref = 0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
!  ...check if mdle is in the fiber core 
      call find_domain(mdle, ndom)
      select case(ndom)
         case(1,2)
            nr_elem_ref = nr_elem_ref + 1
            n_elem(nr_elem_ref) = mdle
         case default
            ! do nothing
      end select
   enddo
!
   write(*,2410) nr_elem_ref
 2410 format(' found ', i6, ' elements to refine.')
!
   do iel=1,nr_elem_ref
      mdle = n_elem(iel)
      write(*,2400) mdle
!  ...refine element mdle
      select case(NODES(mdle)%type)
         case('mdlb')
            call refine(mdle,kref_mdlb)
         case('mdlp')
            call refine(mdle,kref_mdlp)
         case default
         write(*,*) 'ref_core: unsupported node type. stop.'
         stop
      end select
   enddo
!
 2400 format(' refining mdle = ', i6)
!
end subroutine ref_core
!
!
!-----------------------------------------------------------------------
!                       subroutine ref_layer
!-----------------------------------------------------------------------
!
subroutine ref_layer(Kref)
!
   use data_structure3D
   use assembly
!
   implicit none
!
   integer, intent(in) :: Kref
!
   integer :: mdle, iel, nr_elem_ref
!
!..geometry dof (work space for nodcor)
   real*8, dimension(3,MAXbrickH) :: xnod
   integer :: n_elem(NRELES)
!
   mdle=0; nr_elem_ref = 0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
!  ...check if mdle is at the fiber input
      call nodcor(mdle, xnod)
      ! For HEXA only
      if (minval(xnod(3,1:8)) .lt. 0.5d0) then
         nr_elem_ref = nr_elem_ref + 1
         n_elem(nr_elem_ref) = mdle
      endif
   enddo
!
   write(*,2510) nr_elem_ref
 2510 format(' found ', i6, ' elements to refine.')
!
   do iel=1,nr_elem_ref
      mdle = n_elem(iel)
! ...refine element mdle
      write(*,2500) mdle
      call refine(mdle,Kref)
   enddo
!
 2500 format(' refining mdle = ', i6)
!
end subroutine ref_layer

