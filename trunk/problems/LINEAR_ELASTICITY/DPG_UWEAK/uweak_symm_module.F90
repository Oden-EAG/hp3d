!--------------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module uweak_symm_module

  use parametersDPG
  use parameters
!
!  ...Necessary constant since there is no dynamic allocation
!       integer, parameter :: MAXNRHS_MOD = 4
! !
! !  ...MATRICES
! !     stiffnes and load matrices for the enriched test space
!       real*8, dimension(9*MAXbrickHH,3*MAXbrickH) :: EnrTraceDispl!,EnrTraceDisplc
!       real*8, dimension(9*MAXbrickHH,3*MAXbrickV) :: EnrTraceStress!,EnrTraceStressc
!       real*8, dimension(9*MAXbrickHH,3*MAXbrickQ) :: EnrFieldDispl!,EnrFieldDisplc
!       real*8, dimension(9*MAXbrickHH,6*MAXbrickQ) :: EnrFieldStress!,EnrFieldStressc
!       real*8, dimension(9*MAXbrickHH,MAXNRHS_MOD) :: EnrLoad
!       real*8, dimension(9*MAXbrickHH,3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ) :: EnrStiffness
!       real*8, dimension(9*MAXbrickHH,3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ+MAXNRHS_MOD) :: EnrEverything
!       real*8, dimension(3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ,3*MAXbrickH+3*MAXbrickV+9*MAXbrickQ+MAXNRHS_MOD) :: FullDPG
! !     Gram matrix for the local Riesz matrix in LAPACK format
!       real*8, dimension((9*MAXbrickHH)*(9*MAXbrickHH+1)/2) :: Gram


end module uweak_symm_module
