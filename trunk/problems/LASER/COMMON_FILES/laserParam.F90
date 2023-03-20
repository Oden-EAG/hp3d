!----------------------------------------------------------------------
!
!   module name        - laserParam
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2019
!
!   purpose            - problem dependent data
!
!----------------------------------------------------------------------
!
module laserParam
!
   use parametersDPG
   use commonParam
!
   implicit none
!
!
!..toggle heat flag, nonlinear problem, co/counter pumping, plane pump
   integer :: HEAT_FLAG, ANISO_HEAT, NONLINEAR_FLAG, COPUMP, FAKE_PUMP
!
!..HEAT FLAG = 0
!  ...NONLINEAR_FLAG = 0, Linear Maxwell
!  ...               = 1, Nonlinear Maxwell (gain)
!..HEAT_FLAG = 1
!  ...NONLINEAR_FLAG = 0, Linear Heat
!  ...               = 1, Coupled Heat/Maxwell
!
!..number of heat steps, step size, and max time
   integer :: NSTEPS
   real(8) :: DELTA_T, T_MAX
!
! ....................NON-DIMENSIONALIZATION....................
!
!..Real fiber length [m]
   real(8), parameter :: FIBER_LENGTH_0 = 10.d0
!
!..Dimensionalization constants (Maxwell)
!  ...L_0    : length scale [m]
!  ...omega_0: frequency scale [rad/s]
!  ...I_0    : irradiance scale [W/m^2]
!  ...nu_0   : Yb population concentration scale [ion/m^3]
!  ...sigma_0: Yb absorption/emission crossection scale [m^2/ion]
   real(8), parameter :: L_0     = 1.0d-5
   real(8), parameter :: OMEGA_0 = LIGHT_SPEED/L_0
   real(8), parameter :: I_0     = 1.0d10
   real(8), parameter :: NU_0    = 1.0d25
   real(8), parameter :: SIGMA_0 = 1.0d-26
!
!..Upper level radiative lifetime Yb ions [s]
   real(8), parameter :: TAU = 8.0d-4
!
!..Dimensionalization constants (Heat)
!  ...TIME_0: time scale [s]
!  ...TEMP_0: temperature scale [K]
   real(8), parameter :: TIME_0  = 1.0d-3
   real(8), parameter :: TEMP_0  = 1.0d0
!
!..Heat equation parameters
!  ...KAPPA  : Thermal conductivity of silica [W/(mK)]
!  ...RHO_0  : Mean density of silica [kg/m^3]
!  ...CP_HEAT: Specific heat of silica [Ws/(kg K)]
!  ...T_AMB_0  : Ambient temperature [K]
!  ...THERMO_OPT_COEFF_0: Thermo-optic coefficient of silica [1/K]
   real(8), parameter :: KAPPA   = 1.38d0
   real(8), parameter :: RHO_0   = 2201.d0
   real(8), parameter :: CP_HEAT = 703.d0
   real(8), parameter :: T_AMB_0 = 298.15d0
   real(8), parameter :: THERMO_OPT_COEFF_0 = 1.285d-5
!
!..all values below are non-dimensional quantities,
!  non-dimensionalized using the constants above
!
!..Fiber parameters
!  ...R_CORE: fiber core radius
!  ...R_CLAD: inner fiber cladding radius
   real(8), parameter :: R_CORE = 0.9d0*sqrt(2.0d0)
   real(8), parameter :: R_CLAD = 9.0d0*sqrt(2.0d0)
   !real(8), parameter :: R_CLAD = 14.0d0*sqrt(2.0d0)
!
!..Non-dimensional real fiber length
   real(8), parameter :: FIBER_LENGTH = FIBER_LENGTH_0/L_0
!
!..refractive index of core (n1) and cladding (n2)
   !real(8) :: REF_INDEX_CORE = 1.4515d0
   real(8) :: REF_INDEX_CORE = 1.4512d0
   real(8) :: REF_INDEX_CLAD = 1.4500d0
   real(8) :: NA, VNUM
!
!..support for anisotropic refractive index tensor
   integer :: ANISO_REF_INDEX = 0
   real(8) :: CORE_NX, CORE_NY, CORE_NZ
   real(8) :: CLAD_NX, CLAD_NY, CLAD_NZ
   real(8) :: CORE_N(3,3), CLAD_N(3,3)
!
!..induce artificial index grating (see bgpol routine)
   integer :: ART_GRATING = 0
!
!..Signal and Pump wavelengths [m]
   real(8), parameter :: LAMBDA_SIGNAL = 1064.0d-9 / L_0
   real(8), parameter :: LAMBDA_PUMP   =  976.0d-9 / L_0
!
!..Raman gain wavelengths [m] (Signal=Stokes, Pump)
   !real(8), parameter :: LAMBDA_SIGNAL = 1116.0d-9 / L_0
   !real(8), parameter :: LAMBDA_PUMP   = 1064.0d-9 / L_0
!
!..Signal and Pump angular frequency
   real(8), parameter :: OMEGA_SIGNAL = 2.d0*PI/LAMBDA_SIGNAL
   real(8), parameter :: OMEGA_PUMP   = 2.d0*PI/LAMBDA_PUMP
   real(8), parameter :: OMEGA_RATIO  = LAMBDA_SIGNAL/LAMBDA_PUMP
!
!..Active Yb gain population dynamics parameters
!  ...N_TOTAL    : Total population concentration [ion/m^3]
!  ...SIGMA_S_ABS: Signal absorption cross-section [m^2/ion]
!  ...SIGMA_S_EMS: Signal emission cross-section [m^2/ion]
!  ...SIGMA_P_ABS: Pump absorption cross-section [m^2/ion]
!  ...SIGMA_P_EMS: Pump emission cross-section [m^2/ion]
!  ...SIGMA_P_EMS: Pump emission cross-section [m^2/ion]
!  ...TAU_0      : Non-dimensional value in gain expression
!  ...ACTIVE_GAIN: Non-dimensional amplifier for active gain
   real(8), parameter :: N_TOTAL = 6.0d25 / NU_0
   real(8), parameter :: SIGMA_S_ABS = 6.000d-27 / SIGMA_0
   real(8), parameter :: SIGMA_S_EMS = 3.580d-25 / SIGMA_0
   real(8), parameter :: SIGMA_P_ABS = 1.429d-24 / SIGMA_0
   real(8), parameter :: SIGMA_P_EMS = 1.776d-24 / SIGMA_0
   real(8), parameter :: TAU_0 = H_BAR*OMEGA_0/(I_0*SIGMA_0*TAU)
   real(8)            :: ACTIVE_GAIN = 1.0d2
!
!..Raman gain parameters
!  ...RAMAN_GAIN          : Non-dimensional amplifier for Raman gain
!  ...UPSILON_RAMAN_SIGNAL:
!  ...UPSILON_RAMAN_PUMP  :
!  ...OMEGA_RATIO_SIGNAL  :
!  ...OMEGA_RATIO_PUMP    :
   real(8)            :: RAMAN_GAIN            = 1.0d-3
   real(8), parameter :: UPSILON_RAMAN_SIGNAL  = 1.d0
   real(8), parameter :: UPSILON_RAMAN_PUMP    = -OMEGA_RATIO
   real(8), parameter :: OMEGA_RATIO_SIGNAL    = 1.d0
   real(8), parameter :: OMEGA_RATIO_PUMP      = OMEGA_RATIO
!
!..Heat equation parameters
!  ...ALPHA_0         : non-dimensional diffusivity scaling
!  ...Q_0             : non-dimensional heat source scaling
!  ...T_AMB           : Ambient temperature
!  ...THERMO_OPT_COEFF: Thermo-optic coefficient
   real(8), parameter :: ALPHA_0 = KAPPA*TIME_0/(L_0*L_0*RHO_0*CP_HEAT)
   real(8), parameter :: Q_0     = NU_0*SIGMA_0*I_0*TIME_0/(RHO_0*CP_HEAT*TEMP_0)
   real(8), parameter :: T_AMB   = T_AMB_0 / TEMP_0
   real(8), parameter :: THERMO_OPT_COEFF = THERMO_OPT_COEFF_0 * TEMP_0
!
!  ALPHA_Z: non-dimensional short fiber heat diffusivity scaling in z
   real(8) :: ALPHA_Z
!
end module laserParam
