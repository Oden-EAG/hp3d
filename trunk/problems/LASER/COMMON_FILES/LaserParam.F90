!----------------------------------------------------------------------
!
!   module name        - LaserParam
!
!----------------------------------------------------------------------
!
!   latest revision    - Sept 17
!
!   purpose            - problem dependent data
!
!----------------------------------------------------------------------
!
module LaserParam
!
   use parametersDPG
   use CommonParam
!
!..number of heat steps
   integer :: NSTEPS
!
!..toggle laser mode, nonlinear problem, co/counter pumping, line search scaling
   integer :: LASER_MODE, NONLINEAR_FLAG, COPUMP
   real*8  :: LS_SCALE
!
!..for ansatz
!!  ansatz stuff not needed anymore ..
!   double precision, dimension(3) :: ANSATZVEC
!!$OMP THREADPRIVATE(ANSATZVEC)
!   double precision, dimension(3,3) :: DANSATZVEC
!!$OMP THREADPRIVATE(DANSATZVEC)
!
   integer                  :: AZINUM = 1
   double precision         :: BETA_PROPAGATION
   double precision         :: ORDER_BESSEL
!..until here (eliminate)
   double precision, parameter :: R_CORE = 0.25D0*sqrt(2.0d0)
   double precision, parameter :: R_CLAD = 2.50D0*sqrt(2.0d0)
!
! ............NON-DIMENSIONALIZATION....................
!
!..Parameters for the gain function
   double precision :: REF_INDEX_CORE
   double precision :: REF_INDEX_CLAD
   double precision :: NA,VNUM
   double precision, parameter :: REF_INDEX_AVG = 1.45075d0
!..material density for heat (not used currently)
!..now it's non-dimensionalized..
   REAL*8 RHO_0
   PARAMETER (RHO_0 = 2201.d0)
!..same for heat capacity
   REAL*8 HEAT_CAPACITY
   PARAMETER (HEAT_CAPACITY = 703.d0)
!..tau was used for active gain, not used currently
   REAL*8 TAU
   PARAMETER (TAU = 8.014D-4)
!..not currently used b/c of non-dimensionalization
   REAL*8 MU0
   PARAMETER (MU0 = 1.2566370614D-6)
   REAL*8 LIGHT_SPEED
   PARAMETER (LIGHT_SPEED = 2.99792458D8)
!..core rad not used (we have non-dim. number r_core and r_clad above)
   REAL*8 CORE_RAD
   PARAMETER (CORE_RAD = 11.232155D-6)
!..n0 not used anymore (was for active gain)
   REAL*8 N0
   PARAMETER (N0 = 8.75D31)
!..hbar not used anymore
   REAL*8 H_BAR
   PARAMETER (H_BAR = 1.05457266D-34)
   double precision, parameter :: LAMBDA_SIGNAL = 1116.0D-9
   double precision, parameter :: LAMBDA_PUMP = 1064.0D-9
   double precision, parameter :: OMEGA_SIGNAL = 2.d0*PI*LIGHT_SPEED/LAMBDA_SIGNAL
   double precision, parameter :: OMEGA_PUMP = 2.d0*PI*LIGHT_SPEED/LAMBDA_PUMP
   double precision, parameter :: OMEGA_RATIO = OMEGA_PUMP/OMEGA_SIGNAL
   double precision, parameter :: I_INC = 1000.d0/(3.14159*CORE_RAD**2)
   double precision, parameter :: E_INC = sqrt(I_INC*MU0*LIGHT_SPEED/REF_INDEX_AVG)
   double precision :: E_0
   REAL*8, parameter :: SIGMA_S_ABS    = 6.0D-27
   REAL*8, parameter :: SIGMA_S_EM     = 3.58D-25
   REAL*8, parameter :: SIGMA_P_ABS    = 1.429D-24
   REAL*8, parameter :: SIGMA_P_EM     = 1.776D-24
!
! .... RAMAN GAIN PARAMETERS
   REAL*8            :: RAMAN_GAIN    = 1.0D-3
   REAL*8, parameter :: UPSILON_RAMAN_SIGNAL     = 1.d0
   REAL*8, parameter :: UPSILON_RAMAN_PUMP    = -OMEGA_RATIO
   REAL*8, parameter :: OMEGA_RATIO_SIGNAL     = 1.d0
   REAL*8, parameter :: OMEGA_RATIO_PUMP    = OMEGA_RATIO
!
! .... HEAT EQUATION PARAMETERS
   REAL*8            :: KAPPA
   REAL*8            :: ALPHA_HEAT
   REAL*8            :: DELTAT, TMAX, HEAT_LOAD, NTHETA
!
end module LaserParam
