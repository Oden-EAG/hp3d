!
!-------------------------------------------------------------------------------
!> @brief       Initialize workspace and parameters; read geometry; initialize mesh
!!
!> @date        July 2023
!-------------------------------------------------------------------------------
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
!..set GMP and hp3D parameters
!                             NDIM //   MANDIM //
   call set_gmp_parameters(      3,         3,             &
!                            MAXSU //    MAXNP //   MAXNC //
                                30,     10000,     10000,  &
!                            MAXTR //    MAXRE //   MAXBT //
                            10000,      10000,     10000,  &
!                            MAXHE //    MAXTE //   MAXPY //
                            10000,      10000,         1)
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
!=======================================================================
!
!                        NRCOMS // MAXNRHS //
   call set_parameters(      1 ,        1 ,  &
!                       MAXEQNH // MAXEQNE // MAXEQNV // MAXEQNQ //
                             0 ,        2,         0,         6)
!
!..read geometry file
   call read_geometry(trim(FILE_GEOM))
!
   NRELIS = NRPRISM+NRHEXAS+NRTETRA+NRPYRAM
!
!..Overwrite MAXNODS if specified by user input via argument list
   if (MAXNODS_USER .gt. 0) MAXNODS = MAXNODS_USER
!
!..generate mesh and read physics file
   call hp3gen(trim(FILE_PHYS))
!
 end subroutine initialize