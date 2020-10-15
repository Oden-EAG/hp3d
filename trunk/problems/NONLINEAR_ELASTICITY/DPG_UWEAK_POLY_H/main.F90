!-----------------------------------------------------------------------
!> Purpose : a code for elasticity
!-----------------------------------------------------------------------
!
program main
!
  use environment
  use common_prob_data
  use parameters, only : NSTD_OUT
  use data_structure3D_poly , only : NRELES, NRDOFSH
  ! use metis_tests
!
!-----------------------------------------------------------------------
  implicit none
  real*8  :: t,greedy,err,rnorm,rvoid,factor
  integer :: i,ic,idec,iter,nvoid,istep,nsteps,nreflag,nstop,idec_solve,idec_test
  logical :: solved
  integer(kind=8) :: t1,t2,clock_rate,clock_max
!-----------------------------------------------------------------------
!                             INITIALIZATION
!
! Set common hp3d environment parameters (reads in options arguments)
  call begin_environment  ! <-- found inside src/modules/environment.F90
!
! Set environment parameters specific to LINEAR_ELASTICITY
! This sets default values for: FILE_REFINE,FILE_VIS,VLEVEL,
!                               FILE_CONTROL,FILE_GEOM,FILE_ERR,
!                               FILE_HISTORY,FILE_PHYS
  call set_environment_poly  ! <-- found inside ../common/set_environment.F90
!
! Exit if this is a "dry run".
  call end_environment  ! <-- found inside src/modules/environment.F90







  ! idec = 0

  ! write(*,*) '***POLYDPG 3D***'
  ! write(*,*) 'main: Run code verification tests? Yes (1)  No (any other integer) '
  ! read(*,*) idec

  ! if (idec.eq.1) then
  !   call initialize_verify

  !   idec_test=1
  !   do while(idec_test /= 0)
  !     write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
  !     write(*,*) 'POLYDPG VERIFICATION UNIT TESTS'
  !     write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
  !     write(*,*) 'SELECT'
  !     write(*,*) 'QUIT ................................... 0'
  !     write(*,*) 'ALL .................................... 1'
  !     write(*,*) 'DIRICHLET BOUNDARY LIFT ................ 2'
  !     write(*,*) 'POLYHEDRAL ELEMENT CONNECTIVITY......... 3'
  !     write(*,*) 'GEOMETRIC UTILITIES FOR POLYHEDRA....... 4'
  !     write(*,*) 'INTEGRATION IN 2D....................... 5'
  !     write(*,*) 'INTEGRATION IN 3D....................... 6'
  !     write(*,*) 'SHAPE FUNCTIONS......................... 7'
  !     write(*,*) 'REDUCED ELEMENT DOFS.................... 8'
  !     write(*,*) 'STATIC CONDENSATION..................... 9'
  !     write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='

  !     read(*,*) idec_test


  !     select case(idec_test)

  !     case(0)
  !       write(*,*) 'Leaving verification'
  !       write(*,*) '******'
  !       ! goto 1111
  !     case(1)
  !       write(*,*) 'Starting ALL verification tests'
  !       write(*,*) '******'
  !       write(*,*) 'Starting polyhedral geometric utilities verification'
  !       call poly_geom_util_verify
  !       write(*,*) '******'
  !     case(2)
  !       write(*,*) 'Starting Dirichlet boundary data projection verification'

  !     case(3)
  !       write(*,*) 'Starting polyhedral connectivity verification'

  !     case(4)
  !       write(*,*) 'Starting polyhedral geometric utilities verification'
  !       call poly_geom_util_verify
  !       write(*,*) '******'

  !     case(5)
  !       write(*,*) 'Starting 2D-integration verification'

  !     case(6)
  !       write(*,*) 'Starting 3D-integration verification'

  !     case(7)
  !       write(*,*) 'Starting shape function verification'

  !     case(8)
  !       write(*,*) 'Starting reduced element dofs (celem logic) verification'

  !     case(9)
  !       write(*,*) 'Starting static condensation verification'


  !     end select
    
  !   enddo
  !   call finalize_verify
  ! else 
  !   write(*,*) 'Skipping verification'
  !   write(*,*) '******'
  ! endif
! 
! 
! 
! 
! 
! 
!   
!
!     print fancy header
  write(*,*)'                      '
  write(*,*)'// ------  POLYDPG FOR LINEAR ELASTICITY  ------ //'
  write(*,*)'// ------      WITH CONTINUOUS TRACES     ------ //'
  write(*,*)'                      '
!
! Initialize common library (set common parameters, load solvers, 
!                                              and create initial mesh)
  call initialize  ! <-- found inside ../common/initialize.F90
!-----------------------------------------------------------------------
!                            INTERACTIVE MODE
!
!   display menu
  solved = .FALSE.
  idec=1
!   display menu in infinite loop
 do while(idec /= 0)
!
    write(*,*) '                                         '
    write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    write(*,*) 'SELECT'
    write(*,*) 'QUIT ................................... 0'
    write(*,*) '                                         '
    ! write(*,*) 'Geometry graphics (X11) ................1'
    write(*,*) 'PolyHP3D graphics (X11) ................ 2'
    write(*,*) 'Paraview ............................... 3'
    write(*,*) 'Generate MATLAB visualization data ..... 8'
    ! write(*,*) '                                         '
    ! write(*,*) 'Print Data Structure arrays ...........10'
    write(*,*) 'Dumpout Data Structure .................11'
    ! write(*,*) '                                         '
    ! write(*,*) ' -- Refinements --                       '
    ! write(*,*) 'Single Uniform H-refinement ...........20'
    ! write(*,*) 'Single Adaptive H-Refinements .........21'
    ! write(*,*) 'Multi-step Uniform H-refinement .......22'
    ! write(*,*) 'Multi-step Adaptive H-refinement ......23'
    ! write(*,*) 'Interactive H-refinements .............21'
    ! write(*,*) 'Interactive P-enrichments .............22'
    ! write(*,*) 'Uniform     H-refinements .............23'
    ! write(*,*) 'Uniform     P-enrichments .............24'
    write(*,*) '                                         '
    write(*,*) ' -- Solver --                            '
    ! write(*,*) 'FRONTAL SOLVE PROBLEM .................30'
    write(*,*) 'MUMPS SOLVE (OMP) CONTINUOUS TRACES.....30'
    write(*,*) 'MUMPS SOLVE (OMP) DISCONTINUOUS TRACES..40'
    ! write(*,*) 'UHM SOLVE (PAR.) ......................50'
    write(*,*) ' -- Results --                            '
    write(*,*) '                                         '
    write(*,*) 'RESIDUAL ...............................50'
    write(*,*) 'EXACT ERROR ............................60'
!    write(*,*) 'Write error to file ...................61'
    write(*,*) 'METIS testing.......................... 90'
    write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    read( *,*) idec
!
    select case(idec)
!  QUIT
    case( 0) 
      ! call finalize ! <-- found inside ../common/finalize.F90
      stop
!
!-----------------------------------------------------------------------
!                                 GRAPHICS
!
! !  GMP x11 graphics
!     case(1) ; call graphg
! !
!  hp3d x11 graphics
    case(2) ; call graphb
! !
!  Paraview graphics
    case(3) ; call paraview_driver
!
!   MATLAB
    case(8)
      call generate_visual_poly_elem_hni
      call generate_visual_poly_face_ct
!-----------------------------------------------------------------------
!                                 DATA STRUCTURE 
!
!  print data structure
    ! case(10) ; call result
!
!  dump out
    case(11) ; call dumpout
!-----------------------------------------------------------------------
!                               REFINEMENTS
!
!
! !   Single uniform refinement
!     case(20)
!       if (.not.solved) write(*,*) 'You have not solved.'
!       if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
!       call refine_DPG(IUNIFORM,1,0.25d0, nstop)
!       solved=.FALSE.
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
! !   Multi-step uniform h refinement
!     case(22)  
!       nsteps=0    
!       do while (nsteps.le.0)
!         write(*,*) 'NUMBER OF REFINEMENTS>0'
!         write(*,*) 'SOLVER:frontal=1,MUMPS=2,UHM=3'
!         write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
!         read(*,*) nsteps,idec_solve
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
!           end select
!         endif
! !    ...say it has solved and save results to paraview file
!         solved=.TRUE.
!         call paraview_driver
! !    ...then display error and refine if necessary
!         if (i.ne.nsteps) then
!           call refine_DPG(IUNIFORM,1,0.25d0, nstop)
!           if (nstop.eq.1) then
!             write(*,*) 'No elements were refined.'
!             write(*,7000) i
!  7000       format('Exiting loop after ',i2,' refinements...')
!             cycle
!           else
!             solved=.FALSE.
!           endif
!         else ! Last step only display (no refinement)
!           call refine_DPG(INOREFINEMENT,1,0.25d0, nstop)
!         endif
!       enddo
! !   Multi-step adaptive h refinement
!     case(23)  
!       nreflag=0
!       do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
!         write(*,*) 'NUMBER OF REFINEMENTS>0'
!         write(*,*) 'REFINEMENT FLAG:h-refine=1,p-refine=2,hp-refine=3'
!         write(*,*) '0.d0<FACTOR<1.d0'
!         write(*,*) 'SOLVER:frontal=1,MUMPS=2,UHM=3'
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
!           end select
!         endif
! !    ...say it has solved and save results to paraview file
!         solved=.TRUE.
!         call paraview_driver
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
!       enddo
!
!-----------------------------------------------------------------------
!                               O L D
!-----------------------------------------------------------------------
! !   Adaptive refinements
!     case(20)      
!       nreflag=0
!       do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
!         write(*,7011)
!  7011   format('main: SET INITIAL, FINAL STEP,',  &
!                ' REFINEMENT FLAG (1-href,2-pref,3-hpref), FACTOR, idec_solve')
!         read(*,*) ibeg,iend,nreflag,factor,idec_solve
!       enddo
!       istep=ibeg
!       do while(istep.lt.iend)
!         ! nrdof_old = NRDOFSH
!         call adapt_DPG(idec_solve,istep,nreflag,factor, nstop)
!         if (nstop.eq.1) exit
!         ! This is already checked inside adaptDPG
!  !        nrdof_new = NRDOFSH !this changes within data_structure3D after being refined inside adapt
!  !        if (nrdof_new.eq.nrdof_old) then
!  !          istep=istep-1
!  !          factor = factor*0.25d0
!  !          write(*,7023) factor
!  ! 7023       format('main: NEW factor = ',e12.5)
!  !          if (factor.lt.0.000001d0) exit
!  !        endif
!         istep=istep+1
!       enddo
! !  ...solve the problem on the current mesh last time (so that you can plot right after)
!       select case(idec_solve)
!       case(1)
!         call solve1(NR_RHS)
!       case(2)
!         call mumps_solve_seq(NR_RHS)
!       case(3)
!         call uhm_solve
!         call uhm_solver_flush(UHM_SOLVER_PTR)
!       end select
!  !    !-------------------------------------------------------------------
!  !    ! interactive p-enrichments
!  !    case(22)
!  !      !
!  !      ! print active elements
!  !      call display_act_elem
!  !      !
!  !      ! select element, refinement kind
!  !      write(*,*) 'Select: mdle,ip,iq,ir ='
!  !      read( *,*) mdle,ip,iq,ir
!  !      write(*,7012) mdle,ip,iq,ir
!  ! 7012 format(' Refining mdle,ip,iq,ir = ',i4,2x,3(i2,2x),' ...')
!  !      !
!  !      ! refine element
!  !      call p_refine(mdle,ip,iq,ir)
!  !      call enforce_min_rule
!  !      call update_gdof
!  !      call update_ddof
!  !      call verify_orient
!  !      call verify_neig
! !-----------------------------------------------------------------------
! !  uniform global h-refinements
!     case(23)
!       ! call global_href
!       ! call update_gdof
!       ! call update_ddof
!       ! call verify_orient
!       ! call verify_neig
!       call refine_DPG(IUNIFORM,1,0.25, nstop)
! !-----------------------------------------------------------------------
! !   global p-enrichments
!     ! case(24) ; call global_pref
!
!-----------------------------------------------------------------------
!                                 SOLVERS
!
! !  Frontal Solver
!     case(30)
!       i=1
!       if (solved) then
!         write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
!         read(*,*) i
!       endif
!       if (i.eq.1) then
!         call solve1(NR_RHS_PROB)
! !    ...say it has solved and save results to paraview file
!         solved=.TRUE.
!         call paraview_driver
!       endif
!
!  MUMPS
    case(30)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call system_clock ( t1, clock_rate, clock_max )
        ! call mumps_sc_poly('H')
        call complete_Ddof_poly_ct
        ! call mumps_sc_symm_poly_ct('H')
        call mumps_sc_poly_ct('H')
        ! call mumps_symm_poly_ct('H')
        call system_clock ( t2, clock_rate, clock_max )
        t =  real(t2 - t1,8)/real(clock_rate,8)
        write(*,*) 'time mumps = ', t
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        ! call paraview_driver
      endif

!  MUMPS
    case(40)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call system_clock ( t1, clock_rate, clock_max )
        ! call mumps_sc_poly('H')
        call mumps_sc_symm_poly('H')
        call system_clock ( t2, clock_rate, clock_max )
        t =  real(t2 - t1,8)/real(clock_rate,8)
        write(*,*) 'time mumps = ', t
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        ! call paraview_driver
      endif


!
! !  UHM
!     case(50)
!       i=1
!       if (solved) then
!         write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
!         read(*,*) i
!       endif
!       if (i.eq.1) then
!         call uhm_time_in
!         call uhm_solve
!         call uhm_time_out(t)
!         write(*,*) 'time uhm = ', t
!         call uhm_solver_flush(UHM_SOLVER_PTR)
! !    ...say it has solved and save results to paraview file
!         solved=.TRUE.
!         call paraview_driver
!       endif
!
!-----------------------------------------------------------------------
!                                RESIDUAL
!
!  Compute residuals locally and globally
    case(50)
      call compute_residual
!
!-----------------------------------------------------------------------
!                                 ERROR
!

!  Compute errors of different variables
    case(60)
      call exact_error!(err,rnorm)
!
!-----------------------------------------------------------------------
!  METIS testing of mesh partitioning
    case(90)
      call all_tests()
    endselect
!
! end infinite loop
  enddo
!
! ! finalize library
!   call finalize ! <-- found inside ../common/finalize.F90
!
end program main