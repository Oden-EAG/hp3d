!----------------------------------------------------------------------
!
!     routine name      - common_prob_data
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2021
!
!     purpose:          - module setting up the parameters for the
!                         vector Poisson problem
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
!..MISCELLANEOUS
   integer :: NR_RHS_PROB = 1
!
!..INITIAL ORDER OF APPROXIMATION
   integer :: IP
!
!..BOUNDARY CONDITION
   integer :: IBC_PROB
   integer, parameter :: BC_NONE        = 0
   integer, parameter :: BC_DIRICHLET   = 1
   integer, parameter :: BC_IMPEDANCE   = 3
   integer, parameter :: BC_CAVITY      = 4
   integer, parameter :: BC_CAVITY_SCAT = 5
   integer, parameter :: BC_SPHERE_SCAT = 6
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
   integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP
   integer :: IDOMAIN_SMOOTHE = 0
!
   character(32) :: OUTPUT_DIR
!
end module common_prob_data
