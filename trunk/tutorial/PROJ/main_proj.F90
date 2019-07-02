!------------------------------------------------------------------------
!> Purpose : a code for Mesh modification and its verification
!------------------------------------------------------------------------
!
program main
  !
  use environment
  use parameters, only : NSTD_OUT
  use data_structure3D , only : NRELES
  use uhm2

  use proj
  !------------------------------------------------------------------------
  implicit none
  character(len=1024) :: argv
  integer :: i, idec, iter, mdle, kref, ipx,ipy,ipz
  real*8 :: err, rnorm, t
  !------------------------------------------------------------------------
  ! default 
  IEXACT_DISP   = 2
  ICHOOSE_DISP  = 1
  H1_PROJ       = 0
  IPX = 1; IPY =0; IPZ = 0;

  call begin_environment
  call set_problem
  call end_environment
  call initialize
  
  UHM_SOLVER_REPORT_SHOW = .TRUE.
  call get_command(argv)
  call uhm_initialize(argv)
  call uhm_option_begin
  call uhm_option_end

  call uhm_direct_lu_piv_initialize( &
       UHM_DOUBLE, NR_RHS_PROB, 256, UHM_SOLVER_PTR)

  ! diplay menu
  idec = 1
  do while (idec>0)
     call menu
     read(*,*) idec
     select case(idec)
     case(0);  exit
     case(1);  call graphg
     case(2);  call graphb
     case(10); call result
     case(11); 
        call global_href
        call close
     case(12)
        mdle=0
        do i=1,NRELES
           call nelcon(mdle, mdle)
           write(*,7000) mdle,NODES(mdle)%type
7000       format(' mdle,type = ',i4,2x,a4)        
        enddo
        call display_ref_kind
        write(*,*) 'SET mdle,kref'
        read(*,*) mdle,kref
        call refine(mdle,kref)
     case(13)
        write(*,*) 'SET % to refine and # of iteration'
        read(*,*) per, niter
        call random_refine(per,niter)
     case(14)      
        call close
        call update_gdof
     case(15)
        write(*,*) 'SET NPX, NPY, NPZ'
        read(*,*) IPX, IPY, IPZ

     case(20); call solve1(NR_RHS_PROB)
     case(30); call mumps_solve_seq(NR_RHS_PROB)
     case(40); call uhm_solve; call uhm_solver_flush(UHM_SOLVER_PTR);
     case(50); call geometry_error(err, rnorm)
        write(*,7001) err, rnorm, err/rnorm
     case(51); call exact_error(err, rnorm)
        write(*,7001) err, rnorm, err/rnorm
     case(52); call calculate_volume(vol, smin, smax)
        write(*,7002)  vol, smin, smax
     case(53); call verify_orient
     case(54); call verify_neig
     end select
  enddo
7001 format('main: |u-u_hp| = ', e12.5,2x,' |u| = ', e12.5,2x, ' |u-u_hp|/|u| = ', e12.5)
7002 format('main: vol ', e12.5,2x, ' min size= ', e12.5,2x, ' max size= ', e12.5)  
  !  
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
  write(*,*) 'GEOMETRY GRAPHICS (X11) ................1'
  write(*,*) 'hp3d GRAPHICS (X11) ....................2'
  write(*,*) '   '
  write(*,*) 'PRINT DATA STRUCTURE ..................10'
  write(*,*) 'UNIFORM H-REFINEMENT ..................11'
  write(*,*) 'INTERACTIVE H-REFINEMENT ..............12'
  write(*,*) 'RANDOM H-REFINEMENT ...................13'
  write(*,*) 'CLOSE MESH WITH UPDATE TO DOFS ........14'
  write(*,*) '   '
  write(*,*) 'SET PROJ PROBLEM ......................15'
  write(*,*) '   '
  write(*,*) 'FRONTAL SOLVE PROBLEM .................20'
  write(*,*) 'MUMPS SOLVE (SEQ.) ....................30'
  write(*,*) 'UHM SOLVE (PAR.) ......................40'
  write(*,*) '   '
  write(*,*) 'GEOMETRY ERROR ........................50'
  write(*,*) 'PROJECTION ERROR ......................51'
  write(*,*) 'VOLUME OF DOMAIN ......................52'
  write(*,*) 'ORIENTATION CHECK .....................53'
  write(*,*) 'NEIGHBOR CHECK ........................54'
  write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
end subroutine menu
