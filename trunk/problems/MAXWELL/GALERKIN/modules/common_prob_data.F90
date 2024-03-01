!----------------------------------------------------------------------
!
!     routine name      - common_prob_data
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!> @brief         - module setting up the parameters for the
!                         Maxwell problem
!
!----------------------------------------------------------------------
!
module common_prob_data
!
   implicit none
!
   save
!
!------------------------------------------------------------------------------
!
!..TYPE OF JOB SUBMISSION
!  0: interactive (usual main file)
!  1: execute pre-defined job script
   integer :: JOB
!
!..User can specify parameter MAXNODS via argument list instead of input file
   integer :: MAXNODS_USER
!
!..INITIAL ORDER OF APPROXIMATION
   integer :: IP
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
!..imaginary unit
   complex(8), parameter :: ZI = (0.d0,1.d0)
!
!..pi for sinusoidal solution
   real(8), parameter :: PI = 4.d0*datan(1.d0)
!
!..PHYSICAL PARAMETERS
!   Angular frequency: OMEGA
!   Permittivity     : EPSILON
!   Permeability     : MU
!   Conductivity     : SIGMA
   real(8) :: OMEGA   = 1.d0 * PI
   real(8) :: EPSILON = 1.d0
   real(8) :: MU      = 1.d0
   real(8) :: SIGMA   = 0.d0
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
!
end module common_prob_data
