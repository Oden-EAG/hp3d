!--------------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module uweak_module

   use parametersDPG
   use parameters
!
   implicit none
!
!  ...Necessary constant since there is no dynamic allocation
      integer, parameter :: NRRHS_MOD = 1
!
!  ...MATRICES
!     stiffnes and load matrices for the enriched test space
      real(8), dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickH) :: EnrTraceDispl
!$OMP THREADPRIVATE (EnrTraceDispl)
      real(8), dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickV) :: EnrTraceStress
!$OMP THREADPRIVATE (EnrTraceStress)
      real(8), dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickQ) :: EnrFieldDispl
!$OMP THREADPRIVATE (EnrFieldDispl)
      real(8), dimension(3*MAXbrickVV+3*MAXbrickHH,6*MAXbrickQ) :: EnrFieldStress
!$OMP THREADPRIVATE (EnrFieldStress)
      real(8), dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickQ) :: EnrFieldOmega
!$OMP THREADPRIVATE (EnrFieldOmega)
      real(8), dimension(3*MAXbrickVV+3*MAXbrickHH,NRRHS_MOD) :: EnrLoad
!$OMP THREADPRIVATE (EnrLoad)
      real(8), dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickH+3*MAXbrickV+12*MAXbrickQ) :: EnrStiffness
!$OMP THREADPRIVATE (EnrStiffness)
      real(8), dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickH+3*MAXbrickV+12*MAXbrickQ+NRRHS_MOD) :: EnrEverything
!$OMP THREADPRIVATE (EnrEverything)
      real(8), dimension(3*MAXbrickH+3*MAXbrickV+12*MAXbrickQ,3*MAXbrickH+3*MAXbrickV+12*MAXbrickQ+NRRHS_MOD) :: FullDPG
!$OMP THREADPRIVATE (FullDPG)
!     Gram matrix for the local Riesz matrix in LAPACK format
      real(8), dimension((3*MAXbrickVV+3*MAXbrickHH)*(3*MAXbrickVV+3*MAXbrickHH+1)/2) :: Gram
!$OMP THREADPRIVATE (Gram)


end module uweak_module
