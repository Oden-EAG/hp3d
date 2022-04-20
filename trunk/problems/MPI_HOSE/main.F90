!-----------------------------------------------------------------------
!> Purpose : a code for elasticity
!-----------------------------------------------------------------------
!
program main
!
  use environment
  use common_prob_data
  use parameters,       only : NSTD_OUT
  use data_structure3D, only : NRELES, NRDOFSH
  use uhm2
  use m_assembly,       only : MTime
  use paraview,         only : PARAVIEW_DOMAIN
  ! use environment,      only : PREFIX
!
!-----------------------------------------------------------------------
  implicit none
  real*8  :: t,greedy,err,rnorm,rvoid,factor
  integer :: i,ic,idec,iter,nvoid,istep,nsteps,nreflag,nstop,idec_solve
  logical :: solved
  character(len=128) :: PREFIX_tmp
!-----------------------------------------------------------------------
!                             INITIALIZATION
!
! Set common hp3d environment parameters (reads in options arguments)
  call begin_environment  ! <---- found inside src/modules/environment.F90
!
! Set environment parameters specific to LINEAR_ELASTICITY
! This sets default values for: FILE_REFINE,FILE_VIS,VLEVEL,FILE_CONTROL,FILE_GEOM,
!                               FILE_ERR,FILE_HISTORY,FILE_PHYS
  call set_environment  ! <---- found inside ../common/set_environment.F90
!
! Exit if this is a "dry run".
  call end_environment  ! <---- found inside src/modules/environment.F90
!
!     print fancy header
      write(*,*)'                      '
      write(*,*)'//                                             //'
      write(*,*)'// --         HYBRID DPG FORMULATION        -- //'
      write(*,*)'//                                             //'
      write(*,*)'// -- SHEATHE PROBLEM FOR LINEAR ELASTICITY -- //'
      write(*,*)'//                                             //'
      write(*,*)'                                                 '
! Initialize common library (set common parameters, load solvers, and create initial mesh)
  call initialize  ! <---- found inside ../common/initialize.F90
!
! Mod out by infinitesimal rigid body motions
! IBC_PROB : 0 - uniform traction ; 1 - clamped ends ; 2 - free ends ; 3 - periodic ends
  if ((IBC_PROB.eq.2).or.(IBC_PROB.eq.3)) then
  ! if ((IBC_PROB.eq.2)) then
    call remove_RBM(IBC_PROB)
    ! call update_gdof
    call update_Ddof
  endif
!-----------------------------------------------------------------------
!                            INTERACTIVE MODE
!   display menu
  solved = .FALSE.
  idec=1
!   display menu in infinite loop
  do while(idec /= 0)
  !
    write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    write(*,*) 'SELECT'
    write(*,*) 'QUIT ...................................0'
    write(*,*) '                                         '
    write(*,*) 'Geometry graphics (X11) ................1'
    write(*,*) 'HP3D graphics (X11) ....................2'
    write(*,*) 'Paraview ...............................3'
    ! write(*,*) 'Compute BC Interpolation Error .........5'
    write(*,*) '                                         '
    write(*,*) 'Print Data Structure arrays ...........10'
    write(*,*) 'Dumpout Data Structure ................11'
    write(*,*) '                                         '
    write(*,*) ' -- Refinements --                       '
    write(*,*) 'Single Uniform H-refinement ...........20'
    ! write(*,*) 'Single Adaptive H-Refinements .........21'
    write(*,*) 'Single Custom H-Refinement ............22'
    write(*,*) 'Multi-step Uniform H-refinement .......23'
    ! write(*,*) 'Multi-step Adaptive H-refinement ......24'
    write(*,*) 'Multi-step Custom H-refinement ........25'
    write(*,*) '                                         '
    write(*,*) ' -- Solves --                            '
    write(*,*) 'FRONTAL SOLVE PROBLEM .................30'
    write(*,*) 'MUMPS SOLVE (SEQ.) ....................40'
    write(*,*) 'MUMPS SOLVE (OMP) .....................45'
    write(*,*) 'UHM SOLVE (PAR.) ......................50'
    write(*,*) '                                         '
    write(*,*) 'EXACT ERROR ...........................60'
    ! write(*,*) 'Write error to file ...................61'
    write(*,*) 'RESIDUAL ..............................70'
    write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    read( *,*) idec
  !
    select case(idec)
!  QUIT
    case( 0) ; call finalize ; stop
!-----------------------------------------------------------------------
!                                 GRAPHICS
!  GMP x11 graphics
    case(1) ; call graphg
!
!  hp3d x11 graphics
    case(2) ; call graphb
!
!  Paraview graphics
    case(3) ; call paraview_driver
              ! call paraview_stress_driver
!
!  BC Interpolation Error
    ! case(5) ; call compute_BC_interp_error
!-----------------------------------------------------------------------
!                                 ?
!  print data structure
    case(10) ; call result
!
!  dump out
    case(11) ; call dumpout
!-----------------------------------------------------------------------
!                               REFINEMENTS
!-----------------------------------------------------------------------
!
!   Single uniform refinement
    case(20)
      if (.not.solved) write(*,*) 'You have not solved.'
      if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
      call refine_DPG(IUNIFORM,1,0.25d0, nstop)
      solved=.FALSE.
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
!   Single custom refinement
    case(22)
      if (.not.solved) write(*,*) 'You have not solved.'
      if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
      call refine_DPG(ICUSTOMREF,1,0.25d0, nstop)
      solved=.FALSE.
!   Multi-step uniform h refinement
    case(23)
      nsteps=0
      do while (nsteps.le.0)
        write(*,*) 'NUMBER OF REFINEMENTS>0'
        write(*,*) 'SOLVER:frontal=1,MUMPS=2,UHM=3,OMP MUMPS=4'
        write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
        read(*,*) nsteps,idec_solve
      enddo
      do i=0,nsteps
!    ...solve first if needed
        if (.not.solved) then
          select case(idec_solve)
          case(1)
            call solve1(NR_RHS_PROB)
          case(2)
            call mumps_solve_seq(NR_RHS_PROB)
          case(3)
            call uhm_solve
            call uhm_solver_flush(UHM_SOLVER_PTR)
          case(4)
            call mumps_interf(NR_RHS_PROB)
          end select
        endif
!    ...declare solved and save results to paraview file
        solved=.TRUE.
        ! call paraview_driver
        ! PREFIX_tmp=PREFIX
        ! PARAVIEW_DOMAIN=1
        ! PREFIX="Steel_"
        ! call paraview_custom_driver
        ! PARAVIEW_DOMAIN=2
        ! PREFIX="Rubber_"
        ! call paraview_custom_driver
        ! PREFIX=PREFIX_tmp
        call paraview_custom_dump
!    ...then display error and refine if necessary
        if (i.ne.nsteps) then
          call refine_DPG(IUNIFORM,1,0.25d0, nstop)
          if (nstop.eq.1) then
            write(*,*) 'No elements were refined.'
            write(*,7000) i
 7000       format('Exiting loop after ',i2,' refinements...')
            cycle
          else
            solved=.FALSE.
          endif
        else ! Last step only display (no refinement)
          call refine_DPG(INOREFINEMENT,1,0.25d0, nstop)
        endif
      enddo
! !   Multi-step adaptive h refinement
!     case(24)
!       nreflag=0
!       do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
!         write(*,*) 'NUMBER OF REFINEMENTS>0'
!         write(*,*) 'REFINEMENT FLAG:h-refine=1,p-refine=2,hp-refine=3'
!         write(*,*) '0.d0<FACTOR<1.d0'
!         write(*,*) 'SOLVER:frontal=1,MUMPS=2,UHM=3,OMP MUMPS=4'
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
!           case(4)
!             call mumps_interf(NR_RHS_PROB)
!             write(*,8000) MTime(1), MTime(2), MTime(3), MTime(4)
!  8000 format('-- SOLVE --',/,  &
!              'Step 1: determine the first dof offsets for active nodes',/,  &
!              'Time: ',f15.8,/,  &
!              'Step 2: compute element matrices, assemble global stiffness matrix',/,  &
!              'Time: ',f15.8,/,  &
!              'Step 3: MUMPS solve',/,  &
!              'Time: ',f15.8,/,  &
!              'Step 4: reconstruct global to local connectivities',/,  &
!              'Time: ',f15.8)
!           end select
!         endif
! !    ...say it has solved and save results to paraview file
!         solved=.TRUE.
!         call paraview_driver
!         call paraview_stress_driver
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
!         write(*,8001) MTime(5)
!  8001 format('-- REFINE --',/,  &
!              'Calculate residual :',/,  &
!              'Time: ',f15.8)
!       enddo
!   Multi-step uniform h refinement
    case(25)
      nsteps=0
      do while (nsteps.le.0)
        write(*,*) 'NUMBER OF REFINEMENTS>0'
        write(*,*) 'SOLVER:frontal=1,MUMPS=2,UHM=3,OMP MUMPS=4'
        write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
        read(*,*) nsteps,idec_solve
      enddo
      do i=0,nsteps
!    ...solve first if needed
        if (.not.solved) then
          select case(idec_solve)
          case(1)
            call solve1(NR_RHS_PROB)
          case(2)
            call mumps_solve_seq(NR_RHS_PROB)
          case(3)
            call uhm_solve
            call uhm_solver_flush(UHM_SOLVER_PTR)
          case(4)
            call mumps_interf(NR_RHS_PROB)
          end select
        endif
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        ! call paraview_driver
        ! call paraview_stress_driver
!    ...then display error and refine if necessary
        if (i.ne.nsteps) then
          call refine_DPG(ICUSTOMREF,1,0.25d0, nstop)
          if (nstop.eq.1) then
            write(*,*) 'No elements were refined.'
            write(*,9000) i
 9000       format('Exiting loop after ',i2,' refinements...')
            cycle
          else
            solved=.FALSE.
          endif
        else ! Last step only display (no refinement)
          call refine_DPG(INOREFINEMENT,1,0.25d0, nstop)
        endif
      enddo
!-----------------------------------------------------------------------
! !  uniform global h-refinements
!     case(20)
!       call global_href
!       call update_gdof
!       call update_ddof
!       call verify_orient
!       call verify_neig
!       solved=.FALSE.
!-----------------------------------------------------------------------
!   global p-enrichments
    ! case(24) ; call global_pref
!-----------------------------------------------------------------------
!                                 SOLVERS
!  Frontal Solver
    case(30)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call solve1(NR_RHS_PROB)
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        call paraview_driver
        ! call paraview_stress_driver
      endif
!
!  MUMPS
    case(40)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call uhm_time_in
        call mumps_solve_seq(NR_RHS_PROB)
        call uhm_time_out(t)
        write(*,*) 'time mumps = ', t
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        call paraview_driver
        ! call paraview_stress_driver
      endif
!
!  OMP MUMPS
    case(45)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call uhm_time_in
        call mumps_interf(NR_RHS_PROB)
        call uhm_time_out(t)
        write(*,*) 'time mumps = ', t
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        call paraview_driver
        ! call paraview_stress_driver
      endif
!
!  UHM
    case(50)
      i=1
      if (solved) then
        write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
        read(*,*) i
      endif
      if (i.eq.1) then
        call uhm_time_in
        call uhm_solve
        call uhm_time_out(t)
        write(*,*) 'time uhm = ', t
        call uhm_solver_flush(UHM_SOLVER_PTR)
!    ...say it has solved and save results to paraview file
        solved=.TRUE.
        call paraview_driver
        ! call paraview_stress_driver
      endif
!-----------------------------------------------------------------------
!                                  ERROR
    case(60)
!      select case(isolver)
!        case(1) ; call solve1(         NR_RHS_PROB)
!        case(2) ; call mumps_solve_seq(NR_RHS_PROB)
!        case(3) ; call uhm_solve
!                  call uhm_solver_flush(UHM_SOLVER_PTR)
!       endselect
      call exact_error
!    case(61) ; call dumpout_error_to_file(FILE_ERR)
    case(70) ; call compute_residual
    endselect
!
! end infinite loop
  enddo
!
! finalize library
  call finalize
!
end program main