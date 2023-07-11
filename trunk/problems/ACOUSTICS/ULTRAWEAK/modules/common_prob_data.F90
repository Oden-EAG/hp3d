!
!----------------------------------------------------------------------
!> @brief      Module containing global variables for UW acoustics problem
!> @date       July 2023
!----------------------------------------------------------------------
   module common_prob_data_UW
!
      save
!
      integer :: NR_RHS_PROB = 1
!
!  ...TYPE OF JOB SUBMISSION
!     0: interactive (usual main file)
!     1: stampede2 skx slurm job batch script
      integer :: JOB
!
!  ...NUMBER OF REFINEMENTS (JOB SCRIPT)
      integer :: IMAX
!
!  ...MAXIMUM NUMBER OF NODES
      integer :: MAXNODS_USER
!
!  ...ORDER OF APPROXIMATION (set_initial_mesh.F90)
      integer :: IP
!
!  ...PROBLEM KIND
      integer :: PROB_KIND
      integer, parameter :: PROB_FREESPACE      = 0
      integer, parameter :: PROB_CAVITY         = 1

!  ...BOUNDARY CONDITION     (set_initial_mesh.F90)
      integer :: IBC_PROB
      integer, parameter :: BC_NONE        = 0
      integer, parameter :: BC_DIRICHLET   = 1
      integer, parameter :: BC_NEUMANN     = 2
      integer, parameter :: BC_IMPEDANCE   = 3
      integer, parameter :: BC_CAVITY      = 4
!   
!  ...EXACT SOLUTION (exact.F90)
      integer :: IEXACT_PROB
      integer, parameter :: IEXACT_TRIVIAL     = 0
      integer, parameter :: IEXACT_SINUSOIDAL  = 1
      integer, parameter :: IEXACT_POLYNOMIAL  = 2
      integer, parameter :: IEXACT_PLANE       = 3
      integer, parameter :: IEXACT_POINT       = 4
      integer, parameter :: IEXACT_GAUSS       = 5
!
!  ...exponents defining exact polynomial solution (IEXACT_PROB = IEXACT_POLYNOMIAL)
      integer :: NP1 = 3, NP2 = 2, NP3 = 2
!
!  ...pi for sinusoidal solution
      real*8,  parameter :: PI = 4.d0*datan(1.d0)
!
!  ...tolerance for singular solution
      real*8,  parameter :: EPS = 1.d-10
!
!  ...REFINEMENT TYPE (refine_DPG.F90)
      integer, parameter :: INOREFINEMENT  = 0
      integer, parameter :: IUNIFORM       = 1
      integer, parameter :: IADAPTIVE      = 2
!   
!  ...DISPLAY SETTINGS (soldis.F90 ...)
      integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP
!
!  ...number of wavelengths and angular frequency for Helmholtz
      real*8 :: RNUM, OMEGA
!
!  ...scaling coefficient for the adjoint test norm
      real*8 :: ALPHA
!
   end module common_prob_data_UW
