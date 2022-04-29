!----------------------------------------------------------------------
!
!     routine name      - common_prob_data
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2020
!
!     purpose:          - module setting up the parameters for the
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

!..DPG
   integer :: TEST_NORM
   integer, parameter :: ADJOINT_GRAPH  = 1
   integer, parameter :: MATHEMATICIANS = 2   

!..INITIAL ORDER OF APPROXIMATION
   integer :: IP
!
!..PROBLEM KIND
   integer :: PROB_KIND
   integer, parameter :: PROB_FREESPACE      = 0
   integer, parameter :: PROB_CAVITY         = 1
   integer, parameter :: PROB_SPHERE         = 2
   integer, parameter :: PROB_SCAT_CAVITY    = 3
   integer, parameter :: PROB_SCAT_SPHERE    = 4
   integer, parameter :: PROB_SCAT_CUBE_PML  = 5
   integer, parameter :: PROB_SCAT_SPHR_PML  = 6
   integer, parameter :: PROB_SCAT_PLANE_PML = 7

!..BOUNDARY CONDITION
   integer :: IBC_PROB
   integer, parameter :: BC_NONE        = 0
   integer, parameter :: BC_DIRICHLET   = 1
   integer, parameter :: BC_NEUMANN     = 2
   integer, parameter :: BC_IMPEDANCE   = 3
   integer, parameter :: BC_CAVITY      = 4
   integer, parameter :: BC_CAVITY_SCAT = 5
   integer, parameter :: BC_SPHERE_SCAT = 6
!
!..EXACT SOLUTION
   integer :: ISOL
   integer, parameter :: IEXACT_POLYNOMIAL_X  = 0
   integer, parameter :: IEXACT_SINUSOIDAL    = 1
   integer, parameter :: IEXACT_POLYNOMIAL    = 2
   integer, parameter :: IEXACT_PLANEWAVE     = 3
   integer, parameter :: IEXACT_POINTSOURCE   = 4
   integer, parameter :: IEXACT_GAUSSIANBEAM  = 5
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

!..PHYSICAL PARAMETERS
!   Angular frequency: OMEGA
   real(8) :: OMEGA   = 1.d0 * PI

!..REFINEMENT TYPE (refine_DPG.F90)
   integer, parameter :: INOREFINEMENT  = 0
   integer, parameter :: IUNIFORM       = 1
   integer, parameter :: IADAPTIVE      = 2
!
!..DISPLAY SETTINGS (soldis.F90 ...)
!..DISPLAY SETTINGS (paraview   ...)
   integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP
   integer :: IDOMAIN_SMOOTHE = 0
!
   character(32) :: OUTPUT_DIR
!
!
end module common_prob_data
