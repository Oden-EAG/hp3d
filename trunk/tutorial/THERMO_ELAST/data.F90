!-----------------------------------------------------------------------
!> Purpose:          - routine returns body force and heat source
!!                     data
!! @param[in]   Mdle - middile node
!! @param[in]   X    - physical coordinates of a point
!! @param[out]  Zf   - body force
!! @param[out]  Zq   - heat source
!-----------------------------------------------------------------------
#include "implicit_none.h"
subroutine getf(Mdle,X, Zf,Zq)
  !
  use parameters
  use thermo_elast
  implicit none
  integer,              intent(in)  :: Mdle
  real*8, dimension(3), intent(in)  :: X
  VTYPE,                intent(out) :: Zf(3,MAXNRHS),Zq(MAXNRHS)
  !
  Zf=ZERO; Zq=ZERO
  !
end subroutine getf

!-----------------------------------------------------------------------
!> Purpose:              - routine returns material data
!! @param[in]   Mdle     - middle node
!! @param[in]   X        - physical coordinates of a point
!! @param[out]  Zmu,Zlam - elasticity constants
!! @param[out]  Zk       - thermal conductivity
!! @param[out]  Zalpha   - thermal expansion coeeficient
!-----------------------------------------------------------------------
#include "implicit_none.h"
subroutine getmat(Mdle,X, Zmu,Zlam,Zk,Zalpha)
  !
  use thermo_elast
  implicit none
  integer,              intent(in)  :: Mdle
  real*8, dimension(3), intent(in)  :: X
  VTYPE,                intent(out) :: Zmu, Zlam, Zk, Zalpha

  !     mu  =  E/2/(1+niu)
  !     lam =  E*niu/(1+niu)/(1-2*niu)
  !     E=200000; niu=0.3 --> mu=76923.076923077; lam=115384.615384615 
  Zmu  =  76923.076923077d0
  Zlam = 115384.615384615d0
  !
  Zk     = 50.d0
  Zalpha =  1.d-5
  !
end subroutine getmat

!-----------------------------------------------------------------------
!> Purpose:          - routine returns robin boundary
!! @param[in]  Mdle  - middle node
!! @param[in]  Xi    - master element coordinates
!! @param[in]  X     - physical coordinates of a point
!! @param[in]  Rn    - unit outward normal vector
!! @param[out] Zgval - traction
!! @param[out] Zhflx - heat flux
!
!-----------------------------------------------------------------------
#include "implicit_none.h"
subroutine getg(Mdle,ibc,Xi,X,Rn, Zgval,Zhflx)

  use parameters
  use thermo_elast
  implicit none
  integer,              intent(in)  :: Mdle, ibc
  real*8, dimension(3), intent(in)  :: X, Xi, Rn
  VTYPE,                intent(out) :: Zgval(3,MAXNRHS),Zhflx(MAXNRHS)
  !
  real*8 :: rho
  !
  Zgval=ZERO; Zhflx=ZERO
  !
  rho = sqrt(X(2)**2+X(3)**2)
  if (0.08d0 .le. rho  .AND. rho .le. 0.13d0) then
     Zgval(1,1) = -1.d3*Rn(1)
     Zhflx(1)   = 2.d2 
  endif
  !
end subroutine getg
