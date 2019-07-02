!----------------------------------------------------------------------
!
!   program name       - main 
!
!----------------------------------------------------------------------
!
!   latest revision    - Nov 10
!
!   purpose            - main program
!
!----------------------------------------------------------------------
!
program main
  !
  use uhm
  use parameters, only : NSTD_OUT
  use data_structure3D
  use thermo_elast
  implicit none
  !
  character(len=64), dimension(3) :: args
  integer :: i, iel, icnt=0, idec, mdle, kref

  if (iargc().ne.3) stop 100
  do i=1,3
     call getarg(i, args(i))
  enddo

  write(*,*)'==================================================='
  write(*,*)'|              main Thermo Elast                  |'
  write(*,*)'==================================================='
  !
  call initialize(args)
  !
  !----------------------------------------------------------------------
  !  ...display options      
10 continue
  write(*,*) 'SELECT'
  write(*,*) 'QUIT ...................................0'
  write(*,*) '   '
  write(*,*) 'GEOMETRY GRAPHICS (X11) ................1'
  write(*,*) 'hp3d GRAPHICS (X11) ....................2'
  write(*,*) '   '
  write(*,*) 'PRINT DATA STRUCTURE ..................10'
  write(*,*) 'UNIFORM H-ADAPTIVITY ..................11'
  write(*,*) 'INTERACTIVE H-ADAPTIVITY ..............12'
  write(*,*) '   '
  write(*,*) 'FRONTAL SOLVE..........................20'
  write(*,*) 'MUMPS SOLVE SEQ........................21'
  write(*,*) 'UHM SOLVE PAR .........................30'
  !
  read(*,*) idec
  !----------------------------------------------------------------------
  !
  select case(idec)
  case(0);     call finalize
  case(1);     call graphg
  case(2);     call graphb
  case(10);    call result  
  case(11)
     call global_href
     call close
     call update_gdof

     icnt=icnt+1
     write(*,*) 'icnt ', icnt, 'nrnods', NRNODS 

  case(12)
     mdle=0
     do i=1,NRELES
        call nelcon(mdle, mdle)
        write(*,7000)mdle,NODES(mdle)%type
7000    format(' mdle,type = ',i4,2x,a4)
     enddo
     call display_ref_kind
     write(*,*)'SET mdle,kref'
     read(*,*)mdle,kref
     call refine(mdle,kref)
     call close
     call update_gdof
  case(20); call solve1(NR_RHS_PROB)
  case(21); call mumps_solve_seq(NR_RHS_PROB)
  case(30);
     call uhm_method_disp_all(NSTD_OUT)
     write(*,*) 'choose a decomposition = '
     read(*,*) i

     call uhm_initialize( 4, 192, i )

     call uhm_set_nr_rhs(NR_RHS_PROB)
     call uhm_solve_problem
     call uhm_finalize
  end select

  !
  !  ...go back to menu
  goto 10      
  !
  !
end program main
!     
!----------------------------------------------------------------------
!  ...fake subroutines for linking      
!
subroutine customize_geometry
end subroutine customize_geometry
