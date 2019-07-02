!
!
!
subroutine initialize
!
      use environment
      use control
      use GMP
      use data_structure3D
      use refinements
      use frsolmod
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     open history file
      call open_history_file(trim(FILE_HISTORY))
!
!     initialize refinements arrays
      call init_refinements( trim(FILE_REFINE) )
!
!     read control file
      call read_control(     trim(FILE_CONTROL))
!
!     set Geometry Modeling Package parameters
!=======================================================================
!  NDIM   - dimension of the problem  (must be = 3)
!  MANDIM - dimension of the manifold (must be = 3)
!  MAXSU  - maximum anticipated number of surfaces
!  MAXNP  - maximum anticipated number of points
!  MAXNC  - maximum anticipated number of curves
!  MAXTR  - maximum anticipated number of triangles
!  MAXRE  - maximum anticipated number of rectangles
!  MAXBT  - maximum anticipated number of prisms
!  MAXHE  - maximum anticipated number of hexas
!  MAXTE  - maximum anticipated number of tets
!  MAXPY  - maximum anticipated number of pyramids
!=======================================================================
!
!
!///////////////////////////////////////////////////////////////////////
!                               NDIM //   MANDIM //
      call set_gmp_parameters(     3 ,         3 ,               &
!                              MAXSU //    MAXNP //    MAXNC //
                                  30 ,     10000 ,     10000 ,   &
!                              MAXTR //    MAXRE //    MAXBT //
                               10000 ,     10000 ,     10000 ,   &
!                              MAXHE //    MAXTE //    MAXPY //
                               10000 ,     10000 ,     10000    )
!///////////////////////////////////////////////////////////////////////
!
!
!     set HP3D parameters
!=======================================================================
!  NRCOMS  - number of copies of the Data Structure (DS), namely the
!            number of right-hand sides (rhs)
!  MAXEQNH - maximum number of equations posed in H1 [ Example :
!            acoustic & elasticity. (3 + 1) * number of rhs ]
!  MAXEQHE - maximum number of equations posed in H(curl)
!  MAXEQHV - maximum number of equations posed in H(div)
!  MAXEQHQ - maximum number of equations posed in L2
!  MAXNRHS - maximum number of rhs
!
!  Constraints on parameters to avoid catastrophe :
!
!  MAX...  >= 1
!  MAXNRHS >= NRCOMS
!  MAXEQNH >= NRCOMS * number of equations posed in H1
!  MAXEQNE >= NRCOMS * number of equations posed in H(curl)
!  MAXEQNV >= NRCOMS * number of equations posed in H(div)
!  MAXEQNQ >= NRCOMS * number of equations posed in L2
!=======================================================================
!
!
!  ...set GMP and hp3D parameters
      select case(INPUT_FILE)
      case(DUMPIN_LEGACY_,DUMPIN_)
        NRCOMS = 2 ; MAXNRHS = 1 ; MAXEQNH = 2 ; MAXEQNE = 8
                                   MAXEQNV = 2 ; MAXEQNQ = 24
      case default
!                                     NDIM //    MANDIM //
        call set_gmp_parameters(         3 ,          3 , &
!                                    MAXSU //     MAXNP //     MAXNC //
                                        30 ,      500 ,      600 , &
!                                    MAXTR //     MAXRE //     MAXBT //
                                    100 ,     100 ,     100 , &
!                                    MAXHE //     MAXTE //     MAXPY //
                                        50 ,      165 ,       10)
!                            NRCOMS // MAXNRHS //
        call set_parameters(      2 ,        1 ,  &
!                           MAXEQNH // MAXEQNE // MAXEQNV // MAXEQNQ //
                                  2 ,        8,         2,         24)
      end select
!
!///////////////////////////////////////////////////////////////////////
!                          NRCOMS // MAXNRHS // MAXEQNH // MAXEQNE
!      call set_parameters(      1 ,        1 ,        1 ,        1 ,  &
!                                            // MAXEQNV // MAXEQNQ
!                                                      1 ,        1 )
!///////////////////////////////////////////////////////////////////////
!
!
!     read geometry file
      call read_geometry(trim(FILE_GEOM))
!
!     print GMP parameters to screen [OPTIONAL]
      call print_GMP_parameters
!
!     initialize parameters for X11 graphics [OPTIONAL]
!
!                              MXIGTR // MXIGSTR // MXRGTRZ
      call set_x11_workspace( 200000000 ,  200000000 ,  400000000 )

!
!  ...generate initial mesh
      select case(INPUT_FILE)
      case(DUMPIN_LEGACY_)
        write(*,*)'dumping in physics...'
        call dumpin_physics_from_file('files/dumpPHYS' )
        write(*,*)'dumping in hp3d (LEGACY)...'
        !call dumpin_hp3d_LEGACY(      'files/dumpc3D')
      case(DUMPIN_)
        write(*,*)'dumping in physics...'
        call dumpin_physics_from_file('files/dumpPHYS' )
        write(*,*)'dumping in hp3d...'
        call dumpin_hp3d(             'files/dumpc3Dhp')
      case(NETGEN_,DEFAULT_)
        call hp3gen(trim(FILE_PHYS))
      endselect
!
!  ...initialize work space for the frontal solver
      call set_frsol_workspace (10000000)
!
 end subroutine initialize
