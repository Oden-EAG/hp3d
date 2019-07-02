!----------------------------------------------------------------------
!
!   module name        - problem
!
!----------------------------------------------------------------------
!
!   latest revision    - Jul 18
!
!   purpose            - problem dependent data
!
!----------------------------------------------------------------------
!
module CommonParam
!
   use parameters
!
!..i
   complex*16, parameter :: ZI = (0.d0,1.d0)
!
!..pi
   double precision, parameter :: PI = 4.D0*DATAN(1.D0)
!
!..material constants
   double precision :: MU,EPSILON
   complex*16       :: SIGMA
!
!..frequency
   double precision :: OMEGA
!
!..Gaussian pulse waist and intensity ratio
   double precision :: BEAM_WAIST, INTENSITY_RATIO
!
!..impedance constant and IBCFLAG
   double precision :: GAMMA
   integer          :: IBCFLAG
!
!..for PML: length of Z region and PML_REGION, PML_FRACTION
   integer          :: USE_PML
   double precision :: ZL, PML_REGION, PML_FRAC, EXP_COEFF
!
!..additional parameters including those required by the system
   integer :: ORDER_APPROX
   integer :: NPX, NPY, NPZ, ICOMP_EXACT
   integer,parameter :: MY_NR_RHS=1
   integer :: ICHOOSE_DISP,ICHOOSE_COMP,ICHOOSE_SIGPUMP
   integer :: IEXACT_DISP
!
!..Problem Number
   integer :: NO_PROBLEM
!
!..Geometry Number
   integer :: GEOM_NO
!
!..control flag for the test inner product
   integer :: INNER_PRODUCT
!
!..choose case for exact
   integer :: ISOL
!
!..refinement type (see refine_DPG.F90)
   integer, parameter :: INOREFINEMENT  = 0
   integer, parameter :: IUNIFORM       = 1
   integer, parameter :: IANISOTROPIC   = 2
   integer, parameter :: IADAPTIVE      = 3
!
!..paraview parameters
   integer                            :: IPARADAP
   integer, dimension(:), ALLOCATABLE :: IPARATTR
   integer, dimension(:), ALLOCATABLE :: ISEL_PARAVIEW
   integer :: IDOMAIN_SMOOTHE = 0
!
end module CommonParam
