!----------------------------------------------------------------------
!
!   module name        - commonParam
!
!----------------------------------------------------------------------
!
!   latest revision    - June 2021
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
   complex(8), parameter :: ZI = (0.d0,1.d0)
!
!..pi
   real(8), parameter :: PI = 4.D0*DATAN(1.d0)
!
!..Identity matrix
   real(8), dimension(3,3), parameter :: IDENTITY = &
      reshape((/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/), (/ 3, 3 /))
!
!..Speed of light in vacuum
   real(8), parameter :: LIGHT_SPEED = 2.99792458d8
!
!..Planck's constant (hbar=h/2pi)
   real(8), parameter :: H_BAR = 1.05457266d-34
!
!..material constants
   real(8) :: MU,EPSILON
!
!..frequency
   real(8) :: OMEGA
!
!..impedance BC:
!  GAMMA  : impedance constant
!..IBCFLAG: 0 (dirichlet)
!           2 (impedance via penalty method)
!           3 (impedance via elimination)
   real(8) :: GAMMA
   integer :: IBCFLAG
!
!..weight for l2 term in scaled adjoint graph norm (UW Maxwell)
!  ||v||_V^2 := ||A^* v||^2 + ALPHA * ||v||^2
   real(8) :: ALPHA_NORM
!
!..additional parameters including those required by the system
   integer :: IP
   integer :: NPX, NPY, NPZ
   integer :: ICOMP_EXACT
   integer :: ICHOOSE_DISP, ICHOOSE_COMP
   integer :: IEXACT_DISP
   integer, parameter :: MY_NR_RHS=1
!
!..choose case for exact
   integer :: ISOL
!
!..TYPE OF JOB SUBMISSION
!  0: interactive (usual main file)
!  1: stampede2 skx slurm job batch script
   integer :: JOB
!
!..USER DEFINED MAX NUMBER OF NODES
   integer :: MAXNODS_USER
!
!..NUMBER OF REFINEMENTS (JOB SCRIPT)
   integer :: IMAX
   integer :: JMAX
!
!..refinement type (see refine_DPG.F90)
   integer, parameter :: INOREFINEMENT = 0
   integer, parameter :: IUNIFORM      = 1
   integer, parameter :: IADAPTIVE     = 2
!
!..I/O
   character(32) :: OUTPUT_DIR
!
end module commonParam
