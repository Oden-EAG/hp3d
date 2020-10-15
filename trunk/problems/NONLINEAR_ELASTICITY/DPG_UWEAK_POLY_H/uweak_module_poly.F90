!--------------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!--------------------------------------------------------------------------------
!
module uweak_module_poly
!
  use parameters
  use parametersDPG  

  ! implicit none
!
!  ...Necessary constant since there is no dynamic allocation
  integer, parameter :: MAXNRHS_MOD = 1
!
!  ...MATRICES
!     stiffnes and load matrices for the enriched test space
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),(3*MAXtriaQ*MAX_NRFC)) :: EnrTraceDispl!,EnrTraceDisplc
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtriaQ*MAX_NRFC) :: EnrTraceStress!,EnrTraceStressc
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtetraQ) :: EnrFieldDispl!,EnrFieldDisplc
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),6*MAXtetraQ) :: EnrFieldStress!,EnrFieldStressc
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),3*MAXtetraQ) :: EnrFieldOmega!,EnrFieldOmegac
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),MAXNRHS_MOD) :: EnrLoad
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),6*MAXtriaQ*MAX_NRFC+12*MAXtetraQ) :: EnrStiffness
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH),6*MAXtriaQ*MAX_NRFC+12*MAXtetraQ+MAXNRHS_MOD) :: EnrEverything
!   real*8, dimension((6*MAXtriaQ*MAX_NRFC+12*MAXtetraQ),6*MAXtriaQ*MAX_NRFC+12*MAXtetraQ+MAXNRHS_MOD) :: FullDPG
! !     Gram matrix for the local Riesz matrix in LAPACK format
!   real*8, dimension((3*MAXtetraVV+3*MAXtetraHH)*(3*MAXtetraVV+3*MAXtetraHH+1)/2) :: Gram

end module uweak_module_poly
