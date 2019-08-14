!---------------------------------------------------------------------------
!> Purpose : a code for elasticity
!---------------------------------------------------------------------------
!
program main_dome
!
  use dome
  use data_structure3D , only : NRELES
!---------------------------------------------------------------------------
  implicit none
  character(len=64), dimension(3) :: args
  integer :: i, idec, mdle, kref
!---------------------------------------------------------------------------
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
  IEXACT_DISP   = 2 ; ITANGENT_DISP = 0
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
      call update_gdof
    case(12)
      mdle=0
      do i=1,NRELES
        call nelcon(mdle, mdle)
        write(*,7000)mdle,NODES(mdle)%type
 7000   format(' mdle,type = ',i4,2x,a4)        
      enddo
      call display_ref_kind
      write(*,*)'SET mdle,kref'
      read(*,*)mdle,kref
      call refine(mdle,kref)
      call close
      call update_gdof

    case(20); call solve1(NR_RHS_PROB)
    case(21); call mumps_solve_seq(NR_RHS_PROB)
    
    end select
  enddo
!  
end program main_dome
!
!
!---------------------------------------------------------------------------
!> Purpose : display menu
!---------------------------------------------------------------------------
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
  write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
end subroutine menu
