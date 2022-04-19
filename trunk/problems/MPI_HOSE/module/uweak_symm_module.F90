!--------------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module uweak_symm_module

  use parametersDPG
  use parameters
!
!  ...Necessary constant since there is no dynamic allocation
      integer, parameter :: MAXNRHS_MOD = 4
!
!  ...MATRICES
!     stiffnes and load matrices for the enriched test space
      real*8, dimension(9*MAXbrickHH,3*MAXbrickH) :: EnrTraceDispl!,EnrTraceDisplc
!$OMP THREADPRIVATE (EnrTraceDispl)
      real*8, dimension(9*MAXbrickHH,3*MAXbrickV) :: EnrTraceStress!,EnrTraceStressc
!$OMP THREADPRIVATE (EnrTraceStress)
      real*8, dimension(9*MAXbrickHH,3*MAXbrickQ) :: EnrFieldDispl!,EnrFieldDisplc
!$OMP THREADPRIVATE (EnrFieldDispl)
      real*8, dimension(9*MAXbrickHH,6*MAXbrickQ) :: EnrFieldStress!,EnrFieldStressc
!$OMP THREADPRIVATE (EnrFieldStress)
      real*8, dimension(9*MAXbrickHH,MAXNRHS_MOD) :: EnrLoad
!$OMP THREADPRIVATE (EnrLoad)
      real*8, dimension(9*MAXbrickHH,3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ) :: EnrStiffness
!$OMP THREADPRIVATE (EnrStiffness)
      real*8, dimension(9*MAXbrickHH,3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ+MAXNRHS_MOD) :: EnrEverything
!$OMP THREADPRIVATE (EnrEverything)
      real*8, dimension(3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ,3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ+MAXNRHS_MOD) :: FullDPG
!$OMP THREADPRIVATE (FullDPG)
!     Gram matrix for the local Riesz matrix in LAPACK format
      real*8, dimension((9*MAXbrickHH)*(9*MAXbrickHH+1)/2) :: Gram
!$OMP THREADPRIVATE (Gram)


end module uweak_symm_module
