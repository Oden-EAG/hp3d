!----------------------------------------------------------------------
!
!     routine name      - common_prob_data
!
!----------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          - module setting up the parameters for the
!                         acoustics problem
!
!
!----------------------------------------------------------------------
!
module common_prob_data
   save
!
!------------------------------------------------------------------------------
!
!..MATLAB DIR
   character(len=128) :: MATLAB_DIR
!..MISCELLANEOUS
!
   integer :: NR_RHS_PROB = 1
!
!..DPG
   integer :: TEST_NORM
   integer, parameter :: ADJOINT_GRAPH  = 1
   integer, parameter :: MATHEMATICIANS = 2
!
!..ORDER OF APPROXIMATION (set_initial_mesh.F90)
   integer :: IP
!
!..PROBLEM KIND
   integer :: PROB_KIND
   integer, parameter :: PROB_FREESPACE   = 0
   integer, parameter :: PROB_CAVITY      = 1
   integer, parameter :: PROB_SPHERE      = 2
   integer, parameter :: PROB_SCAT_CAVITY = 3
   integer, parameter :: PROB_SCAT_SPHERE = 4
   integer, parameter :: PROB_FICHERA     = 5
!
!..BOUNDARY CONDITION     (set_initial_mesh.F90)
   integer :: IBC_PROB
   integer, parameter :: BC_NONE        = 0
   integer, parameter :: BC_DIRICHLET   = 1
   integer, parameter :: BC_IMPEDANCE   = 3
   integer, parameter :: BC_CAVITY      = 4
   integer, parameter :: BC_CAVITY_SCAT = 5
   integer, parameter :: BC_SPHERE_SCAT = 6
!
!..EXACT SOLUTION (exact.F90)
   integer :: IEXACT_PROB
   integer, parameter :: IEXACT_POLYNOMIAL  = 0
   integer, parameter :: IEXACT_SINUSOIDAL  = 1
   integer, parameter :: IEXACT_POLYNOMIAL0 = 2
   integer, parameter :: IEXACT_PLANE       = 3
   integer, parameter :: IEXACT_POINT       = 4
   integer, parameter :: IEXACT_GAUSS       = 5
   integer, parameter :: IEXACT_GAUSS_SIMP  = 6
   integer, parameter :: IEXACT_POLYNOMIAL1 = 7
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
!..DISPLAY SETTINGS - paraview
   integer, dimension(10) :: ISEL_PARAVIEW = (/0,0,0,0,0,1,1,1,1,1/)
!
!..control variables/parameters of geometrical transformation on
!..point coordinates
   logical :: COORD_TRANS        = .false.
   logical :: COORD_TRANS_TRAS   = .false.
   logical :: COORD_TRANS_ROT    = .false.
   logical :: COORD_TRANS_SCAL   = .false.
   real(8), dimension(3) :: TRAS_VECTOR = (/1.d0,1.d0,1.d0/)
   real(8), dimension(3) :: ROT_ANGLES = (/10.d0,45.d0,60.d0/)
!
!..ELEMENT CALCULATIONS (elem.F90)
   real(8), parameter :: SYMMETRY_TOL = 1.d-9
!
!..number of wavelengths and angular frequency
   real(8) :: RNUM, OMEGA

!..permittivity (ε) and permeability (μ)
   real(8) :: EPS, MU

!..scaling coefficient for the adjoint test norm
   real(8) :: ALPHA
!
!
end module common_prob_data
