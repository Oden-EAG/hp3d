!-----------------------------------------------------------------------
!> Purpose : a code for elasticity
!-----------------------------------------------------------------------
!
program main
!
  use environment
  use testvars
  use common_prob_data
  use parameters!, only : NSTD_OUT
  use data_structure3D , only : NRELES, NRDOFSH
  use uhm2
!
!-----------------------------------------------------------------------
  implicit none
  real*8  :: t,greedy,err,rnorm,rvoid,factor
  integer :: i,ic,idec,iter,nvoid,ibeg,istep,iend,nreflag,nstop,  &
             nrdof_old,nrdof_new
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
      write(*,*)'// --  DPG METHOD  -- //'
      write(*,*)'                      '
!
! Initialize common library (set common parameters, load solvers, and create initial mesh)
  call initialize  ! <---- found inside ../common/initialize.F90
!
! 
!-----------------------------------------------------------------------
!                            INTERACTIVE MODE
!   display menu
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
    ! write(*,*) 'Paraview ...............................3'
    write(*,*) '                                         '
    write(*,*) 'Print Data Structure arrays ...........10'
    write(*,*) 'Dumpout Data Structure ................11'
    write(*,*) '                                         '
    write(*,*) ' -- Refinements --                       '
    write(*,*) 'Adaptive    H-Refinements (DPG) .......20 (DEBUGGING)'
    ! write(*,*) 'Interactive H-refinements .............21'
    ! write(*,*) 'Interactive P-enrichments .............22'
    write(*,*) 'Uniform     H-refinements .............23'
    ! write(*,*) 'Uniform     P-enrichments .............24'
    write(*,*) '                                         '
    write(*,*) ' -- Solves --                            '
    write(*,*) 'FRONTAL SOLVE PROBLEM .................30'
    write(*,*) 'MUMPS SOLVE (SEQ.) ....................40'
    write(*,*) 'UHM SOLVE (PAR.) ......................50'
    write(*,*) '                                         '
    write(*,*) 'EXACT ERROR ...........................60'
    write(*,*) 'Write error to file ...................61'
    write(*,*) 'RESIDUAL ..............................70'
    write(*,*) '                                         '
    write(*,*) 'teststuff .............................80'
    write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    read( *,*) idec
  !
    select case(idec)
!  QUIT
    case( 0) ; call finalize ; stop  ! <---- found inside ../common/finalize.F90
!-----------------------------------------------------------------------
!                                 GRAPHICS
!  GMP x11 graphics
    case(1) ; call graphg
!
!  hp3d x11 graphics
    case(2) ; call graphb
!
!  Paraview graphics
    ! case(3) ; call paraview_driver
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
!  H-adaptive refinements
    case(20)
 333  write(*,7011)
 7011 format('main: SET INITIAL, FINAL STEP,',  &
             ' REFINEMENT FLAG (1-href,2-pref,3-hpref), FACTOR')
      read(*,*) ibeg,iend,nreflag,factor
      if (nreflag.ne.1.and.nreflag.ne.2.and.nreflag.ne.3) go to 333
      istep=ibeg-1
      do while(istep.lt.iend)
        istep=istep+1
        nrdof_old = NRDOFSH
        call adapt_DPG(istep,nreflag,factor,nstop)
        if (nstop.eq.1) exit
        nrdof_new = NRDOFSH
        if (nrdof_new.eq.nrdof_old) then
          istep=istep-1
          factor = factor*0.25d0
          write(*,7023) factor
 7023       format('main: NEW factor = ',e12.5)
          if (factor.lt.0.000001d0) exit
        endif
      enddo
 !    !-------------------------------------------------------------------
 !    ! interactive p-enrichments
 !    case(22)
 !      !
 !      ! print active elements
 !      call display_act_elem
 !      !
 !      ! select element, refinement kind
 !      write(*,*) 'Select: mdle,ip,iq,ir ='
 !      read( *,*) mdle,ip,iq,ir
 !      write(*,7012) mdle,ip,iq,ir
 ! 7012 format(' Refining mdle,ip,iq,ir = ',i4,2x,3(i2,2x),' ...')
 !      !
 !      ! refine element
 !      call p_refine(mdle,ip,iq,ir)
 !      call enforce_min_rule
 !      call update_gdof
 !      call update_ddof
 !      call verify_orient
 !      call verify_neig
!-----------------------------------------------------------------------
!  uniform global h-refinements
    case(23)
      call global_href
      call update_gdof
      call update_ddof
      call verify_orient
      call verify_neig
!-----------------------------------------------------------------------
!   global p-enrichments
    ! case(24) ; call global_pref
!-----------------------------------------------------------------------
!                                 SOLVERS
!  Frontal Solver
    case(30) ; call solve1(NR_RHS_PROB)
!
!  MUMPS
    case(40)
      call uhm_time_in
      call mumps_solve_seq(NR_RHS_PROB)
      call uhm_time_out(t)
      write(*,*) 'time mumps = ', t
!
!  UHM (2)
    case(50)
      call uhm_time_in
      call uhm_solve
      call uhm_time_out(t)
      write(*,*) 'time uhm = ', t
      call uhm_solver_flush(UHM_SOLVER_PTR)
!-----------------------------------------------------------------------
!                                  ERROR
    case(60)
!      select case(isolver)
!        case(1) ; call solve1(         NR_RHS_PROB)
!        case(2) ; call mumps_solve_seq(NR_RHS_PROB)
!        case(3) ; call uhm_solve
!                  call uhm_solver_flush(UHM_SOLVER_PTR)
!       endselect
      call exact_error(err,rnorm)
    case(61)
      call dumpout_error_to_file(FILE_ERR)
!
    case(80) !teststuff  
      call teststuff
    endselect  
!
! end infinite loop
  enddo
!
! finalize library
  call finalize ! <---- found inside ../common/finalize.F90
!
end program main
