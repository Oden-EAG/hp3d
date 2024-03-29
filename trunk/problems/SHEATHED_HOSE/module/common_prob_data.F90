!> @brief Define problem specific parameters and constants
!> @rev - Feb.2012
!------------------------------------------------------------------------------
module common_prob_data
!
   implicit none
!
   save
!
!..TYPE OF JOB SUBMISSION
!  0: interactive (usual main file)
!  1: stampede2 skx slurm job batch script
   integer :: JOB
!
!..User can specify parameter MAXNODS via argument list instead of input file
   integer :: MAXNODS_USER

!------------------------------------------------------------------------------
! DPG
  integer :: TEST_NORM
  integer, parameter :: ADJOINT_GRAPH  = 1
  integer, parameter :: MATHEMATICIANS = 2
  integer, parameter :: H1SYMM         = 3
!------------------------------------------------------------------------------
! ORDER OF APPROXIMATION (set_initial_mesh.F90)
  integer :: IP

!------------------------------------------------------------------------------
! BOUNDARY CONDITION     (set_initial_mesh.F90)
  integer :: IBC_PROB
  integer, parameter :: BC_DIRICHLET = 1
  integer, parameter :: BC_NEUMANN   = 2
  integer, parameter :: BC_MIXED     = 3

!
!..NUMBER OF REFINEMENTS (JOB SCRIPT)
  integer :: IMAX

! pi for sinusoidal solution
  real(8),  parameter :: PI = 4.d0*datan(1.d0)
! tolerance for singular solution
  real(8),  parameter :: EPS = 1.d-10

!------------------------------------------------------------------------------
! ERROR (exact_error.F90)
  integer :: IERROR_PROB
  integer, parameter :: IERROR_L2      = 1
  integer, parameter :: IERROR_NATURAL = 2
  integer, parameter :: IERROR_CUSTOM  = 3
  ! integer, parameter :: IERROR_STRAIN = 4

  integer :: IERROR_ATTR
  integer, parameter :: DISPLACEMENT = 1
  integer, parameter :: STRESS       = 2
  integer, parameter :: COMBINED     = 3
  integer, parameter :: LAGRANGE     = 4
  integer, parameter :: PRESSURE     = 5

!------------------------------------------------------------------------------
! REFINEMENT TYPE (refine_DPG.F90)
  integer, parameter :: INOREFINEMENT  = 0
  integer, parameter :: IUNIFORM   = 1
  integer, parameter :: IADAPTIVE  = 2

!------------------------------------------------------------------------------
! DISPLAY SETTINGS (soldis.F90 ...)
  integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP
  integer :: IDOMAIN_SMOOTHE = 0
! DISPLAY SETTINGS - paraview
  integer, dimension(10) :: ISEL_PARAVIEW = (/0,0,0,0,0,1,1,1,1,1/)

! control variables/parameters of geometrical transformation on
! point coordinates
  logical :: COORD_TRANS        = .false.
  logical :: COORD_TRANS_TRAS   = .false.
  logical :: COORD_TRANS_ROT    = .true.
  logical :: COORD_TRANS_SCAL   = .false.
  double precision, dimension(3):: TRAS_VECTOR = (/1.d0,1.d0,1.d0/)
  double precision, dimension(3):: ROT_ANGLES = (/10.d0,45.d0,60.d0/)

!------------------------------------------------------------------------------
! ELEMENT CALCULATIONS (elem.F90)
  real(8),  parameter :: SYMMETRY_TOL = 1.d-9
!
end module common_prob_data
