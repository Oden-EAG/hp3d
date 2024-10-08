!----------------------------------------------------------------------
!
!     routine name      - common_prob_data
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2023
!
!> @brief         - module setting up the parameters for the
!                         Poisson problem
!
!
!----------------------------------------------------------------------
!
module common_prob_data
   save
!
!------------------------------------------------------------------------------
!
!..TYPE OF JOB SUBMISSION
!  0: interactive (usual main file)
!  1: stampede2 skx slurm job batch script
   integer :: JOB
!
!..User can specify parameter MAXNODS via argument list instead of input file
   integer :: MAXNODS_USER
!
!..INITIAL ORDER OF APPROXIMATION
   integer :: IP
!
!..BOUNDARY CONDITION
   integer :: IBC_PROB
   integer, parameter :: BC_NONE        = 0
   integer, parameter :: BC_DIRICHLET   = 1
!
!..EXACT SOLUTION
   integer :: ISOL
   integer, parameter :: IEXACT_POLYNOMIAL  = 0
   integer, parameter :: IEXACT_SINUSOIDAL  = 1
!
!..NUMBER OF REFINEMENTS (JOB SCRIPT)
   integer :: IMAX
!
!..order of the polynomial exact solution
   integer :: NPX, NPY, NPZ
!
!..pi for sinusoidal solution
   real(8), parameter :: PI = 4.d0*datan(1.d0)
!
!..REFINEMENT TYPE (refine_DPG.F90)
   integer, parameter :: INOREFINEMENT  = 0
   integer, parameter :: IUNIFORM       = 1
   integer, parameter :: IADAPTIVE      = 2
!
!..DISPLAY SETTINGS (soldis.F90 ...)
!..DISPLAY SETTINGS (paraview   ...)
   integer :: ICHOOSE_COMP
   integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP
   integer :: IDOMAIN_SMOOTHE = 0
!
!..DPG test norm
!  - GRAPH_NORM: scaled adjoint graph norm
!  - MATH_NORM : mathematician's norm
   integer, parameter :: GRAPH_NORM = 1
   integer, parameter ::  MATH_NORM = 2
   integer, parameter ::  TEST_NORM = GRAPH_NORM
!
!..weight for L2 term in scaled adjoint graph norm
!  ||v||_V^2 := ||A^* v||^2 + ALPHA * ||v||^2
   real(8), parameter :: ALPHA_NORM = 1.0d0
!
   character(32) :: OUTPUT_DIR
!
!
end module common_prob_data
