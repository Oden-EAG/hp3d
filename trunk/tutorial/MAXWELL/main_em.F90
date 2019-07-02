!----------------------------------------------------------------------
!> Purpose : a code for Maxwell equations
!----------------------------------------------------------------------
!
program main_maxwell
  !
  use uhm
  use em
  use data_structure3D , only : NRELES
  !----------------------------------------------------------------------
  implicit none
  character(len=64), dimension(3) :: args
  integer :: i, idec, mdle, kref
  !----------------------------------------------------------------------
  !
  if (iargc().ne.3) stop 100
  do i=1,3
     call getarg(i, args(i))
  enddo
  !
  ! initialization
  call initialize(args)
  ! 
  ! set default values for graphics
  IEXACT_DISP   = 2 ; ITANGENT_DISP = 1 
  ICOMPLEX_DISP = 0 ; ICHOOSE_DISP  = 1
  !
  ! display menu  
  idec = 1
  do while (idec>0)
     call menu
     read(*,*) idec
     select case(idec)
     case(0);  call finalize
     case(1);  call graphg
     case(2);  call graphb
     case(3);  call plotgeomVTK
     case(10); call result
     case(11)
        call global_href
        call close
        !!call update_gdof
     case(12)
        mdle=0
        do i=1,NRELES
           call nelcon(mdle, mdle)
           write(*,7000)mdle,NODES(mdle)%type
7000       format(' mdle,type = ',i4,2x,a4)        
        enddo
        call display_ref_kind
        write(*,*)'SET mdle,kref'
        read(*,*)mdle,kref
        call refine(mdle,kref)
        call close
        !!call update_gdof

     case(20); call solve1(NR_RHS_PROB)
     case(21); call mumps_solve_seq(NR_RHS_PROB)
     case(30,31,32,33);

        select case (idec-30)
        case(0);       call uhm_initialize( 4, 192, UHM_LDL_INCPIV )
        case(1);       call uhm_initialize( 4, 192, UHM_INTERF_PARDISO )
        case(2);       call uhm_initialize( 4, 192, UHM_INTERF_WSMP )
        case(3);       call uhm_initialize( 4, 192, UHM_INTERF_MUMPS )
        end select

        call uhm_set_nr_rhs(NR_RHS_PROB)
        call uhm_solve_problem
        call uhm_finalize
     end select
  enddo
  !  
end program main_maxwell
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
  write(*,*) 'GEOMETRY GRAPHICS (X11) ................1'
  write(*,*) 'hp3d GRAPHICS (X11) ....................2'
  write(*,*) 'hp3d GRAPHICS (VTK) ....................3'
  write(*,*) '   '
  write(*,*) 'PRINT DATA STRUCTURE ..................10'
  write(*,*) 'UNIFORM H-ADAPTIVITY ..................11'
  write(*,*) 'INTERACTIVE H-ADAPTIVITY ..............12'
  write(*,*) '   '
  write(*,*) 'FRONTAL SOLVE PROBLEM .................20'
  write(*,*) 'MUMPS SOLVE (SEQ.) ....................21'
  write(*,*) 'UHM SOLVE (PAR.) ......................30'
  write(*,*) 'UHM + PARDISO SOLVE (PAR.) ............31'
  write(*,*) 'UHM + WSMP SOLVE (PAR.) ...............32'
  write(*,*) 'UHM + MUMPS SOLVE (SEQ.) ..............33'
  write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
end subroutine menu
