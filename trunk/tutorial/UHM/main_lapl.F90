!------------------------------------------------------------------------
!> Purpose : a code for Poisson equation
!------------------------------------------------------------------------
!
program main
  !
  use environment
  use lapl
  use parameters, only : NSTD_OUT
  use data_structure3D , only : NRELES
  use adaptivity_local

  !------------------------------------------------------------------------
  implicit none
  character(len=1024) :: argv
  integer :: i, idec, iter, mdle, kref, nfile, nr_vert
  real*8, dimension(MAXEQN_PROB) :: err, rnorm
  real*8 :: t

  !------------------------------------------------------------------------
  ! initialization
  call begin_environment
  call set_problem
  call get_option_real( &
       '-greedy', 'H-Adaptivity greedy factor', &
       0.80d0, GREEDY)
  call get_option_int( &
       '-ipara', 'Adaptivity include paraview output(-1, do not export)', &
       -1, IPARA)
  call end_environment
  call initialize

  ! Fortran Mystery or Our code mystery... 
  call get_command(argv)

  !
  ! set defaul values for graphics
  IEXACT_DISP   = 2
  ICHOOSE_DISP  = 1

  !
  ! diplay menu
  idec = 1
  do while (idec>0)
     call menu
     read(*,*) idec
     select case(idec)
     case(0);  exit 
     case(1);  call graphg
     case(2);  call graphb
     case(3);
        ISEL_PARAVIEW(2) = 0
        ISEL_PARAVIEW(4) = 0
        call paraview('h1_sc', 'H1_SC', 1,1)
        
        !ISEL_PARAVIEW(2) = 1
        !ISEL_PARAVIEW(4) = 10
        !call paraview('h1_gr', 'H1_GR', 1,1)

     case(10); call result
     case(11)
        call global_href
        call close
        !call update_gdof
        !call update_ddof
     case(12)
        mdle=0
        do i=1,NRELES
           call nelcon(mdle, mdle)
           write(*,7000)mdle,NODES(mdle)%type
7000       format(' mdle,type = ',i4,2x,a4)        
        enddo
        call display_ref_kind
        write(*,*) 'SET mdle,kref'
        read(*,*) mdle,kref
        call refine(mdle,kref)
        call close
        ! call update_gdof
        ! call update_ddof

     case(20); call solve1(NR_RHS_PROB)
     case(30); call mumps_solve_seq(NR_RHS_PROB)
     case(50);
        call exact_error(err, rnorm)
        write(*,7001) &
             0.d0,0.d0, -1, &
             err, rnorm, err/rnorm

     end select
  enddo

  call finalize

7001 format('max |u-u_hp| =', e12.5,2x, 'min |u-u_hp| =', e12.5,2x, &
          '# of elements =', i5 &
          '|u-u_hp| = ', e12.5,2x,' |u| = ', e12.5,2x, ' |u-u_hp|/|u| = ', e12.5)
end program main
!
!
!----------------------------------------------------------------------
!> Purpose : display menu
!----------------------------------------------------------------------
!
subroutine menu
  write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
  write(*,*) 'SELECT'
  write(*,*) 'QUIT ...................................0'
  write(*,*) '   '
  write(*,*) 'GMP GRAPHICS (X11) .....................1'
  write(*,*) 'hp3d GRAPHICS (X11) ....................2'
  write(*,*) 'PARAVIEW ...............................3'
  write(*,*) '   '
  write(*,*) 'PRINT DATA STRUCTURE ..................10'
  write(*,*) 'UNIFORM H-REFINEMENT ..................11'
  write(*,*) 'INTERACTIVE H-REFINEMENT ..............12'
  write(*,*) '   '
  write(*,*) 'FRONTAL SOLVE PROBLEM .................20'
  write(*,*) 'MUMPS SOLVE (SEQ.) ....................30'
  write(*,*) '   '
  write(*,*) 'EXACT ERROR ...........................50'
  write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='

end subroutine menu
