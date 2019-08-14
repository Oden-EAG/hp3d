!------------------------------------------------------------------------
!> Purpose : a code for Poisson equation
!------------------------------------------------------------------------
!
program main_poisson
!
  use lapl
  use parameters, only : NSTD_OUT
  use data_structure3D , only : NRELES
!------------------------------------------------------------------------
  implicit none
  character(len=64), dimension(3) :: args
  integer :: i, idec, mdle, kref, nfile, nr_vert
!------------------------------------------------------------------------
!
  if (iargc().ne.3) stop 100
  do i=1,3
     call getarg(i, args(i))
  enddo
!  
! initialization
  call initialize(args)
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
    case(0);  call finalize
    case(1);  call graphg
    case(2);  call graphb
    case(3);  call plotgeomVTK
    case(4);
        nfile = 100
        open(unit=nfile, file='plane.vtk', form='formatted')
        !call write_mesh2vtk(nfile, nr_vert)
        !call write_vtk_data(nfile, nr_vert, 1, 1, 0, 0)
        close(nfile)
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
end program main_poisson
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
  write(*,*) 'hp3d GRAPHICS (PARAVIEW)................4'
  write(*,*) '   '
  write(*,*) 'PRINT DATA STRUCTURE ..................10'
  write(*,*) 'UNIFORM H-ADAPTIVITY ..................11'
  write(*,*) 'INTERACTIVE H-ADAPTIVITY ..............12'
  write(*,*) '   '
  write(*,*) 'FRONTAL SOLVE PROBLEM .................20'
  write(*,*) 'MUMPS SOLVE (SEQ.) ....................21'
  write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
end subroutine menu
