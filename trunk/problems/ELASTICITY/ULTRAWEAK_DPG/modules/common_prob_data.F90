
!------------------------------------------------------------------------------
!> @brief      Define elasticity-specific parameters and constants
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   module common_prob_data
!
      save
!
!  ...DPG
      integer :: TEST_NORM
      integer, parameter :: ADJOINT_GRAPH  = 1
      integer, parameter :: MATHEMATICIANS = 2
!
!  ...DPG test norm L2 weight
      real(8) :: ALPHA
!
!  ...initial order of approximation
      integer :: IP
!
!  ...boundary condition
      integer :: IBC_PROB
      integer, parameter :: BC_NONE      = 0
      integer, parameter :: BC_DIRICHLET = 1
      integer, parameter :: BC_NEUMANN   = 2
      integer, parameter :: BC_MIXED     = 8
!
!  ...exact solution
      integer :: IEXACT_PROB
      integer, parameter :: IEXACT_TRILINEAR   = 1
      integer, parameter :: IEXACT_POLYNOMIAL  = 2
      integer, parameter :: IEXACT_EXPONENTIAL = 3
      integer, parameter :: IEXACT_SINUSOIDAL  = 4
      integer, parameter :: IEXACT_SINGULAR    = 5
!
!  ...exponents defining exact polynomial solution
      integer :: NP1 = 3, NP2 = 2, NP3 = 2
!
!  ...pi for sinusiondal solution
      real(8), parameter :: PI = 4.d0*datan(1.d0)
!
!  ...tolerance for singular solution
      real(8), parameter :: EPS = 1.d-10
!
!  ...error
      integer :: IERROR_PROB
      integer, parameter :: IERROR_L2      = 1
      integer, parameter :: IERROR_NATURAL = 2
      integer, parameter :: IERROR_CUSTOM  = 3
!
      integer :: IERROR_ATTR
      integer, parameter :: DISPLACEMENT = 1
      integer, parameter :: STRESS       = 2
      integer, parameter :: COMBINED     = 3
      integer, parameter :: LAGRANGE     = 4
      integer, parameter :: PRESSURE     = 5
!
!  ...refinement type
      integer, parameter :: INOREFINEMENT  = 0
      integer, parameter :: IUNIFORM       = 1
      integer, parameter :: IADAPTIVE      = 2
!
!  ...display settings
      integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP
      integer :: IDOMAIN_SMOOTHE = 0
      integer, dimension(10) :: ISEL_PARAVIEW = (/0,0,0,0,0,1,1,1,1,1/)
!
!  ...tolerance for element computations
      real*8,  parameter :: SYMMETRY_TOL = 1.d-9
!
!  ...command-line defined max nodes
      integer :: MAXNODS_USER
      integer :: JOB
!
   end module common_prob_data
