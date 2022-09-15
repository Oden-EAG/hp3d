!--------------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module near_incomp_module

  use parametersDPG
  use parameters
!
!  ...Necessary constant since there is no dynamic allocation
      integer, parameter :: MAXNRHS_MOD = 1
!
!  ...MATRICES
!     stiffnes matrices (and copies) for the enriched test space
      real*8, dimension(3*MAXbrickHH+MAXbrickQQ,3*MAXbrickH) :: EnrFieldDispl
!$OMP THREADPRIVATE (EnrFieldDispl)
      real*8, dimension(3*MAXbrickHH+MAXbrickQQ,MAXbrickQ)   :: EnrFieldPress
!$OMP THREADPRIVATE (EnrFieldPress)
      real*8, dimension(3*MAXbrickHH+MAXbrickQQ,3*MAXbrickV) :: EnrTraceStress
!$OMP THREADPRIVATE (EnrTraceStress)
      real*8, dimension(3*MAXbrickHH+MAXbrickQQ,MAXNRHS_MOD) :: EnrLoad
!$OMP THREADPRIVATE (EnrLoad)
      real*8, dimension(3*MAXbrickHH+MAXbrickQQ,3*MAXbrickH+3*MAXbrickV+MAXbrickQ) :: EnrStiffness
!$OMP THREADPRIVATE (EnrStiffness)
      real*8, dimension(3*MAXbrickHH+MAXbrickQQ,3*MAXbrickH+3*MAXbrickV+MAXbrickQ+MAXNRHS_MOD) :: EnrEverything
!$OMP THREADPRIVATE (EnrEverything)
!     Gram matrix for the local Riesz matrix in LAPACK format
     real*8, dimension((3*MAXbrickHH+MAXbrickQQ)*(3*MAXbrickHH+MAXbrickQQ+1)/2) :: Gram
!$OMP THREADPRIVATE (Gram)
      real*8, dimension(3*MAXbrickHH,6) :: GramRBM
!$OMP THREADPRIVATE (GramRBM)

end module near_incomp_module
