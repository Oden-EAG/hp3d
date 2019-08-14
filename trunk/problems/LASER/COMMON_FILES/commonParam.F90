!----------------------------------------------------------------------
!
!   module name        - commonParam
!
!----------------------------------------------------------------------
!
!   latest revision    - Oct 2018
!
!   purpose            - problem dependent data
!
!----------------------------------------------------------------------
!
module commonParam
!
   use parameters
!
!..i
   complex*16, parameter :: ZI = (0.d0,1.d0)
!
!..pi
   real*8, parameter :: PI = 4.D0*DATAN(1.d0)
!
!..Identity matrix
   real*8, dimension(3,3), parameter :: IDENTITY = &
      reshape((/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/), (/ 3, 3 /))
!
!..Speed of light in vacuum
   real*8, parameter :: LIGHT_SPEED = 2.99792458d8
!
!..Planck's constant (hbar=h/2pi)
   real*8, parameter :: H_BAR = 1.05457266d-34
!
!..set FAST_INT=1 to activate fast integration for hexahedra
   integer :: FAST_INT = 1 
!
!..material constants
   real*8     :: MU,EPSILON
   complex*16 :: SIGMA
!
!..frequency
   real*8  :: OMEGA
!
!..impedance constant and IBCFLAG
   real*8  :: GAMMA
   integer :: IBCFLAG
!
   integer :: TIMESTEP = 0
!
!..weight for l2 term in scaled adjoint graph norm (UW Maxwell)
   real*8 :: ALPHA_NORM = 1.0d0
!
!..for PML: length of Z region and PML_REGION, PML_FRACTION
   integer :: USE_PML
   real*8  :: ZL, PML_REGION, PML_FRAC, EXP_COEFF
!
!..additional parameters including those required by the system
   integer :: ORDER_APPROX_X,ORDER_APPROX_Y,ORDER_APPROX_Z
   integer :: NPX, NPY, NPZ, ICOMP_EXACT
   integer :: ICHOOSE_DISP, ICHOOSE_COMP, ICHOOSE_SIGPUMP
   integer :: IEXACT_DISP
   integer, parameter :: MY_NR_RHS=1
!
!..Problem Number
!  ...... NO_PROBLEM = 1 Heat Solve (linear)
!  ......            = 2 Heat Solve (coupled with Maxwell)
!  ......            = 3 Maxwell Signal Solve
!  ......            = 4 Maxwell Pump Solve
!  ......            = 5
!  ......            = 6
   integer :: NO_PROBLEM
!
!..Geometry Number
!  ...... GEOM_NO = 1 Unit Cube
!  ......         = 2 Prism Core
!  ......         = 3 Prism Fiber
!  ......         = 4 Hexa Core
!  ......         = 5 Hexa Fiber
!  ......         = 6 Unit Prism
   integer :: GEOM_NO
!
!..control flag for the test inner product
   integer :: INNER_PRODUCT
!
!..choose case for exact
   integer :: ISOL
!
!..refinement type (see refine_DPG.F90)
   integer, parameter :: INOREFINEMENT = 0
   integer, parameter :: IUNIFORM      = 1
   integer, parameter :: IANISOTROPIC  = 2
   integer, parameter :: IADAPTIVE     = 3
!
!..flag for direction of anisotropic refinement
   integer :: ANISO_FLAG
   integer, parameter :: IREFINE_X  = 100;
   integer, parameter :: IREFINE_Y  = 10;
   integer, parameter :: IREFINE_Z  = 1;
   integer, parameter :: IREFINE_XY = 110;
   integer, parameter :: IREFINE_XZ = 101;
   integer, parameter :: IREFINE_YZ = 11;
!
!..I/O
   character*32 :: OUTPUT_DIR
!
!..Maximum number of iterations in nonlinear Solve
   integer, parameter :: MAX_ITER = 30
!
!..paraview parameters
   !integer                            :: IPARADAP
   !integer, dimension(:), ALLOCATABLE :: ISEL_PARAVIEW
   !integer :: IDOMAIN_SMOOTHE = 0
!
end module commonParam
