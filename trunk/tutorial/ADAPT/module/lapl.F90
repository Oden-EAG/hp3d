!---------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!---------------------------------------------------------------------------
!
module lapl
  !
  use error
  use data_structure3D
  use environment
  save
  !
  ! problem dependent parameters
  !-------------------------------------------------------------------------
  ! E Q U A T I O N    S E T T I N G
  integer, parameter :: NRCOMS_PROB  = 1     ! number of component
  integer, parameter :: NR_RHS_PROB  = 1     ! rhs
  integer, parameter :: MAXNRHS_PROB = 1     ! MAX nr rhs
  integer, parameter :: MAXEQNH_PROB = 1     ! MAX H1 var
  integer, parameter :: MAXEQNE_PROB = 1     ! MAX Hcurl var
  integer, parameter :: MAXEQNV_PROB = 1     ! MAX Hdiv var
  integer, parameter :: MAXEQNQ_PROB = 1     ! MAX L2 var
  integer, parameter :: MAXEQN_PROB = &
       max(MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
  !     
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB    = 3     ! Dimension
  integer, parameter :: MANDIM_PROB  = 3     ! Manifold dimension

  integer, parameter :: MAXSU_PROB   = 100   ! surfaces
  integer, parameter :: MAXNP_PROB   = 1000  ! point
  integer, parameter :: MAXNC_PROB   = 1000  ! curve
  integer, parameter :: MAXTR_PROB   = 1000  ! triangle
  integer, parameter :: MAXRE_PROB   = 1000  ! rectangle
  integer, parameter :: MAXBT_PROB   = 1000  ! prism
  integer, parameter :: MAXHE_PROB   = 1000  ! hexa
  integer, parameter :: MAXTE_PROB   = 1000  ! tetra
  integer, parameter :: MAXPY_PROB   = 1000  ! pyramid
  ! 
  ! I N P U T   A R G U M E N T S
  character(len=128) :: &
       FILE_CONTROL, FILE_HISTORY, FILE_REFINE, FILE_VIS, FILE_GEOM, FILE_PHYS,&
       VLEVEL, PREFIX
  !
  ! display setting
  !-------------------------------------------------------------------------
  integer :: IEXACT_DISP, ICHOOSE_DISP
  !
  ! problem data
  integer, parameter :: IERROR_L2 = 1
  integer, parameter :: IERROR_H1 = 2
  integer :: IERROR_PROB

  !
  integer ::  IORDER_TETR_PROB, IORDER_PRIS_PROB
  !
  integer, parameter :: BC_NONE      = 0
  integer, parameter :: BC_DIRICHLET = 1
  integer, parameter :: BC_NEUMANN   = 2
  integer, parameter :: BC_ROBIN     = 3
  integer :: IBC_PROB

  integer, parameter :: IEXACT_HARMONIC_PROB = 1
  integer, parameter :: IEXACT_SHOCK_PROB    = 2
  integer :: IEXACT_PROB

  ! shock problem data
  real*8, parameter :: X0(3) = (/0.25d0, 0.25d0, 0.25d0/)
  real*8, parameter :: ALPHA = 20.d0
  real*8, parameter :: R0 = 0.1732050807568877d1
  
  ! PARAVIEW
  integer :: ISEL_PARAVIEW(10)
  !
contains
  !-------------------------------------------------------------------------
  subroutine disp_soldis(Nstream)
    implicit none
    integer, intent(in) :: Nstream

    write(Nstream,100)
    write(Nstream,310)
    select case (IEXACT_DISP)
    case(1); write(Nstream,101)
    case(2); write(Nstream,102)
    end select
    select case (ICHOOSE_DISP)
    case(1); write(Nstream,301)
    case(2); write(Nstream,302)
    case(3); write(Nstream,303)
    end select
    write(Nstream,310)

100 format('DISPLAY SETUP')
101 format('EXACT SOLUTION IS CHOSEN')
102 format('APPROX SOLUTION IS CHOSEN')
301 format('u IS CHOSEN')
302 format('u_x IS CHOSEN')
303 format('u_y IS CHOSEN')
304 format('u_z IS CHOSEN')
310 format('-----------------------')
  end subroutine disp_soldis


  subroutine mdle2vtk_data_upscale(Nout, Isel, Mdle, Z_harmonic)
    use data_structure3D
    use upscale
    implicit none
    !----------------------------------------------------------------------
    integer, intent(in) :: Nout, Isel(10), Mdle
    real*8, intent(in) :: Z_harmonic
    character(len=4) :: type
    integer :: &
         i, iv, ndom, &
         nodesl(27), norder(19), &
         norientl(27), nedge_orient(12), nface_orient(6)
    real*8 :: &
         t, xi(3), x(3), xnod(3,MAXbrickH), dxdxi(3,3)
    real*8 :: &
         zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE), &
         zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
    real*8 :: &
         zsolH(MAXEQNH),  zgradH(MAXEQNH,3), &
         zsolE(3,MAXEQNE),zcurlE(3,MAXEQNE), &
         zsolV(3,MAXEQNV),zdivV(MAXEQNV), &
         zsolQ(MAXEQNQ)
    real*8 :: &
         zvalH(MAXEQNH), &
         zdvalH(MAXEQNH,3),zd2valH(MAXEQNH,3,3), &
         zvalE(3,MAXEQNE), &
         zdvalE(3,MAXEQNE,3),zd2valE(3,MAXEQNE,3,3), &
         zvalV(3,MAXEQNV), &
         zdvalV(3,MAXEQNV,3),zd2valV(3,MAXEQNV,3,3), &
         zvalQ(MAXEQNQ), &
         zdvalQ(MAXEQNQ,3),zd2valQ(MAXEQNQ,3,3), &
         zval_curlE(3), z_scalar, z_vector(3)
    type(vis) :: vis_obj
    nodesl = 0; norientl = 0; norder = 0;
    nedge_orient = 0; nface_orient = 0;
    call find_domain(Mdle, ndom)
    if (Isel(ndom+4).eq.0) then
       return
    end if

    call elem_nodes(Mdle, nodesl, norientl)
    call find_orient_from_list(type, norientl, nedge_orient,nface_orient)
    call find_order_from_list(type, nodesl, norder)
    xnod = 0.d0; zdofH = ZERO; zdofE = ZERO;

    call nodcor(Mdle, xnod)
    call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)

    type = NODES(Mdle)%type
    vis_obj = vis_on_type(type)
    do iv=1, vis_obj%nr_vert
       call get_vis_point(vis_obj, iv-1, xi)
       call soleval( &
            mdle, xi, &
            nedge_orient, nface_orient, norder, &
            xnod, zdofH,zdofE,zdofV,zdofQ, &
            1, x, dxdxi, &
            zsolH,zgradH,zsolE,zcurlE,zsolv,zdivV,zsolQ)
       call exact(x, 1, &
            zvalH,zdvalH,zd2valH, &
            zvalE,zdvalE,zd2valE, &
            zvalV,zdvalV,zd2valV, &
            zvalQ,zdvalQ,zd2valQ)
       zval_curlE(1) = zdvalE(2,1,3) - zdvalE(3,1,2)
       zval_curlE(2) = zdvalE(3,1,1) - zdvalE(1,1,3)
       zval_curlE(3) = zdvalE(1,1,2) - zdvalE(2,1,1)

       ! Isel 1,2 are reserved
       select case(Isel(3))
       case(0); ! Approx solution
          select case(Isel(4)/10)
          case(0); z_scalar = zsolH(1)        ! H1 real
          case(1); z_vector = zgradH(1, 1:3)  ! H1 grad real
          end select
       case(1); ! Exact solution
          select case(Isel(4)/10)
          case(0); z_scalar = zvalH(1)        ! H1 real
          case(1); z_vector = zdvalH(1, 1:3)  ! H1 grad real
          end select
       end select

       select case(Isel(4)/10)
       case(0);
          z_scalar = z_scalar * Z_harmonic
          write(Nout,9001) real(z_scalar)
       case(1,2,3);
          z_vector = z_vector * Z_harmonic
          write(Nout,9002) real(z_vector)
       end select
    end do
9001 format(f20.15)
9002 format(3(f20.15,2x))
  end subroutine mdle2vtk_data_upscale

  subroutine set_problem

    call get_option_string( &
         '-file-refinement', 'Refinement lookup file location', &
         '../../files/ref', FILE_REFINE)
    call get_option_string( &
         '-file-vis-upscale', 'Visualization upscale file location', &
         '../../files/vis', FILE_VIS)
    call get_option_string( &
         '-vis-level', 'Visualization upscale level (0-3)', &
         '1', VLEVEL)

    call get_option_string( &
         '-file-control', 'Control file', &
         './files/control', FILE_CONTROL)
    call get_option_string( &
         '-file-geometry', 'Geometry file', &
         './files/cube', FILE_GEOM)
    call get_option_string( &
         '-file-phys', 'Physics file', &
         './files/physics', FILE_PHYS)

    call get_option_string( &
         '-file-history', 'History file', &
         './files/history', FILE_HISTORY)

    call get_option_int( &
         '-disp-exact', 'Soldis display exact value 1) exact, 2) approx, 3) exact error 4\
) geometry error', &
         1, IEXACT_DISP)
    call get_option_int( &
         '-disp-vector', 'Soldis display vector 1) sol vector, 2) tangent vector', &
         1, ITANGENT_DISP)
    call get_option_int( &
         '-disp-component', 'Soldis display component of vector 1) u, 2) u_x, 3) u_y 4) u_z', &
         1, ICHOOSE_DISP)

    call get_option_string( &
         '-prefix', 'Prefix to outputfile', &
         'lapl_', PREFIX)

    call get_option_int( &
         '-error', 'Error definition 1) L2 error 2) H1 error', &
         IERROR_L2, IERROR_PROB)

    call get_option_int( &
         '-order-tetr', 'Approximation order for TETR', &
         2, IORDER_TETR_PROB)

    call get_option_int( &
         '-order-pris', 'Approximation order for TETR', &
         2, IORDER_PRIS_PROB)

    call get_option_int( &
         '-bc', 'Boundary 1) dirichlet, 2) neumann', &
         BC_NONE, IBC_PROB)

    call get_option_int( &
         '-exact', 'Exact problem 1) harmonic, 2) shock', &
         IEXACT_HARMONIC_PROB, IEXACT_PROB)

    ISEL_PARAVIEW(1:10) =  0
    ISEL_PARAVIEW(2)    =  0         ! vector value if 1
    ISEL_PARAVIEW(3)    =  0         ! exact value if 1
    ISEL_PARAVIEW(5:10) =  1         ! domain display     

  end subroutine set_problem
  !
  !
end module lapl
