!----------------------------------------------------------------------------
!> @brief   initialize problem dependent environments, solvers, graphics and
!!          create initial mesh.
!> @date    Sep 2023
!----------------------------------------------------------------------------
subroutine initialize
!
   use environment
   use commonParam
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
!..open history file
   !call open_history_file(trim(FILE_HISTORY))
!
!..initialize refinements arrays
   call init_refinements(trim(FILE_REFINE))
!
!..read control file
   call read_control(trim(FILE_CONTROL))
!
!..set Geometry Modeling Package (GMP) parameters
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
!!                            NDIM //   MANDIM //
!   call set_gmp_parameters(     3 ,         3 ,               &
!!                           MAXSU //    MAXNP //    MAXNC //
!                               30 ,     10000 ,     10000 ,   &
!!                           MAXTR //    MAXRE //    MAXBT //
!                            10000 ,     10000 ,     10000 ,   &
!!                           MAXHE //    MAXTE //    MAXPY //
!                            10000 ,     10000 ,     10000    )
!///////////////////////////////////////////////////////////////////////
!..set GMP and hp3D parameters
!                                NDIM //    MANDIM //
   call set_gmp_parameters(         3 ,          3 , &
!                                MAXSU //     MAXNP //     MAXNC //
                                    30 ,      500 ,     600 , &
!                                MAXTR //     MAXRE //     MAXBT //
                                   100 ,      100 ,      10 , &
!                                MAXHE //     MAXTE //     MAXPY //
                                    40 ,        1 ,       1)
!
!     set HP3D parameters
!=======================================================================
!  NRCOMS  - number of solution copies
!  NRRHS   - number of right-hand sides
!=======================================================================
!
!                      NRCOMS, NRRHS
   call set_parameters(     1,     1)
!
!..read geometry file
   call read_geometry(trim(FILE_GEOM))
!
!..print GMP parameters to screen [OPTIONAL]
   !call print_GMP_parameters
!
!..initialize parameters for X11 graphics [OPTIONAL]
!                           MXIGTR  //  MXIGSTR  //  MXRGTRZ
   !call set_x11_workspace( 200000000 ,  200000000 ,  400000000 )
!
!..Overwrite MAXNODS if specified by user input via argument list
   if (MAXNODS_USER .gt. 0) MAXNODS = MAXNODS_USER
   if (.not. QUIET_MODE) then
      !write(*,*) 'User specified MAXNODS value:'
      !write(*,9999) MAXNODS
9999 format(' MAXNODS = ',i12)
   endif
!
!..generate mesh and read physics file
   call hp3gen(trim(FILE_PHYS))
!
!..initialize work space for the frontal solver
   !call set_frsol_workspace (1000000)
!
 end subroutine initialize
