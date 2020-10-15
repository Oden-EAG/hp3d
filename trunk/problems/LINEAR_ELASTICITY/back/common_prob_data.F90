!> Purpose : define problem of plane wave such as exact solutions
!! @rev - Feb.2012
!------------------------------------------------------------------------------
module common_prob_data
  save
!------------------------------------------------------------------------------
! MISCELLANEOUS
  integer :: NR_RHS_PROB = 1
  ! integer :: ISOLVER

!------------------------------------------------------------------------------
! DPG
  integer :: TEST_NORM
  integer, parameter :: ADJOINTGRAPH = 1
  integer, parameter :: H1 = 2
  integer, parameter :: H1SYMM = 3
!------------------------------------------------------------------------------
! ORDER OF APPROXIMATION (set_initial_mesh.F90)
  integer :: IP

!------------------------------------------------------------------------------
! BOUNDARY CONDITION     (set_initial_mesh.F90)
  integer :: IBC_PROB
  integer, parameter :: BC_NONE      = 0
  integer, parameter :: BC_DIRICHLET = 1
  integer, parameter :: BC_NEUMANN   = 2
  ! integer, parameter :: BC_ROBIN     = 3
  integer, parameter :: BC_MIXED     = 8

!------------------------------------------------------------------------------
! EXACT SOLUTION (exact.F90)
  integer :: IEXACT_PROB
  integer, parameter :: IEXACT_TRILINEAR   = 1
  integer, parameter :: IEXACT_POLYNOMIAL  = 2
  integer, parameter :: IEXACT_EXPONENTIAL = 3
  integer, parameter :: IEXACT_SINUSOIDAL  = 4
  integer, parameter :: IEXACT_SINGULAR    = 5

! exponents defining exact polynomial solution (IEXACT_PROB = IEXACT_POLYNOMIAL)
  integer :: NP1 = 3, NP2 = 2, NP3 = 2
! pi for sinusiondal solution
  real*8,  parameter :: PI = 4.d0*datan(1.d0)
! tolerance for singular solution
  real*8,  parameter :: EPS = 1.d-10

!------------------------------------------------------------------------------
! ERROR (exact_error_element.F90)
  integer :: IERROR_PROB
  integer, parameter :: IERROR_L2     = 1
  integer, parameter :: IERROR_H1     = 2
  integer, parameter :: IERROR_ENERGY = 3
  ! integer, parameter :: IERROR_STRAIN = 4

  integer :: IVIS = 0
  real*8  :: RWORK(100,10)

!------------------------------------------------------------------------------
! DISPLAY SETTINGS (soldis.F90 ...)
  integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP
  integer :: ISEL_PARAVIEW(10), IDOMAIN_SMOOTHE = 0

! control variables/parameters of geometrical transformation on
! point coordinates
  logical :: COORD_TRANS        = .FALSE.
  logical :: COORD_TRANS_TRAS   = .FALSE.
  logical :: COORD_TRANS_ROT    = .TRUE.
  logical :: COORD_TRANS_SCAL   = .FALSE.
  double precision, dimension(3):: TRAS_VECTOR = (/1.d0,1.d0,1.d0/)
  double precision, dimension(3):: ROT_ANGLES = (/10.d0,45.d0,60.d0/)

!------------------------------------------------------------------------------
! ELEMENT CALCULATIONS (elem.F90)
  real*8,  parameter :: SYMMETRY_TOL = 1.d-9



contains

! NEEDED FOR PARAVIEW
  subroutine mdle2vtk_data_upscale(Nout, Isel, Mdle, Z_harmonic)
    use data_structure3D
    use upscale
    implicit none
    ! ** Arguments
    !----------------------------------------------------------------------
    integer, intent(in) :: Nout, Isel(10), Mdle
    real*8, intent(in) :: Z_harmonic
    !
    character(len=4) :: type
    !
    integer :: &
         i, iv, ndom, &
         norder(19), nedge_orient(12), nface_orient(6)
    real*8 :: &
         t, xi(3), x(3), xnod(3,MAXbrickH), dxdxi(3,3)
    !
    real*8 :: &
         zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE), &
         zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
    !
    real*8 :: &
         zsolH(MAXEQNH),  zgradH(MAXEQNH,3), &
         zsolE(3,MAXEQNE),zcurlE(3,MAXEQNE), &
         zsolV(3,MAXEQNV),zdivV(MAXEQNV), &
         zsolQ(MAXEQNQ)
    !
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
    !
    type(vis) :: vis_obj
    !
    norder = 0; nedge_orient = 0; nface_orient = 0;
    !
    call find_domain(Mdle, ndom)
    if (Isel(ndom+4).eq.0) then
       return
    end if

    call find_elem_nodes(Mdle, norder, nedge_orient, nface_orient)
    !
    xnod = 0.d0; zdofH = ZERO; zdofE = ZERO;
    call nodcor(Mdle, xnod)
    call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
    !
    type = NODES(Mdle)%type
    vis_obj = vis_on_type(type)
    !
    do iv=1, vis_obj%nr_vert
       call get_vis_point(vis_obj, iv-1, xi)
       !
       call soleval( &
            mdle, xi, &
            nedge_orient, nface_orient, norder, &
            xnod, zdofH,zdofE,zdofV,zdofQ, &
            1, x, dxdxi, &
            zsolH,zgradH,zsolE,zcurlE,zsolv,zdivV,zsolQ)
       !
       call exact(x, 1, &
            zvalH,zdvalH,zd2valH, &
            zvalE,zdvalE,zd2valE, &
            zvalV,zdvalV,zd2valV, &
            zvalQ,zdvalQ,zd2valQ)
       !
       zval_curlE(1) = zdvalE(2,1,3) - zdvalE(3,1,2)
       zval_curlE(2) = zdvalE(3,1,1) - zdvalE(1,1,3)
       zval_curlE(3) = zdvalE(1,1,2) - zdvalE(2,1,1)
       !
       ! Isel 1,2 are reserved
       select case(Isel(3))
       case(0); ! Approx solution
          select case(Isel(4)/10)
          case(0); z_scalar = zsolH(1)        ! H1 real
          case(1); z_vector = zgradH(1, 1:3)  ! H1 grad real
          case(2); z_vector = zsolE(1:3,1)    ! Hcurl grad real
          case(3); z_vector = zcurlE(1:3, 1)  ! Hcurl curl real
          end select
       case(1); ! Exact solution
          select case(Isel(4)/10)
          case(0); z_scalar = zvalH(1)        ! H1 real
          case(1); z_vector = zdvalH(1, 1:3)  ! H1 grad real
          case(2); z_vector = zvalE(1:3, 1)   ! Hcurl grad real
          case(3); z_vector = zval_curlE(1:3) ! Hcurl curl real
          end select
       end select
       !
       select case(Isel(4)/10)
       case(0);
          z_scalar = z_scalar * Z_harmonic
          select case(mod(Isel(4),10))
          case(0); write(Nout,9001) real(z_scalar)
          end select
       case(1,2,3);
          z_vector = z_vector * Z_harmonic
          select case(mod(Isel(4),10))
          case(0); write(Nout,9002) real(z_vector(1:3));
          end select
       end select
    end do
    !
9001 format(f20.15)
9002 format(3(f20.15,2x))
    !
  end subroutine mdle2vtk_data_upscale
end module common_prob_data
