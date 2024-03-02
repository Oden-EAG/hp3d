!--------------------------------------------------------------------------------
!> @brief define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module primal_module
!
   use parametersDPG
   use parameters
!
   implicit none
!
!
!  ...Necessary constant since there is no dynamic allocation
      integer, parameter :: NRRHS_MOD = 1
!      
!  ...MATRICES
!     stiffnes matrices (and copies) for the enriched test space
      real(8), dimension(3*MAXbrickHH,3*MAXbrickH) :: EnrField
!$OMP THREADPRIVATE (EnrField)
      real(8), dimension(3*MAXbrickHH,3*MAXbrickV) :: EnrTrace
!$OMP THREADPRIVATE (EnrTrace)
      real(8), dimension(3*MAXbrickHH,NRRHS_MOD) :: EnrLoad
!$OMP THREADPRIVATE (EnrLoad)
      real(8), dimension(3*MAXbrickHH,3*MAXbrickH+3*MAXbrickV) :: EnrStiffness
!$OMP THREADPRIVATE (EnrStiffness)
      real(8), dimension(3*MAXbrickHH,3*MAXbrickH+3*MAXbrickV+NRRHS_MOD) :: EnrEverything
!$OMP THREADPRIVATE (EnrEverything)
      real(8), dimension(3*MAXbrickH+3*MAXbrickV,3*MAXbrickH+3*MAXbrickV+NRRHS_MOD) :: FullDPG
!$OMP THREADPRIVATE (FullDPG)
!     Gram matrix for the local Riesz matrix in LAPACK format
      real(8), dimension(3*MAXbrickHH*(3*MAXbrickHH+1)/2) :: Gram
!$OMP THREADPRIVATE (Gram)

end module primal_module
