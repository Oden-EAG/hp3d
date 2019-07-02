!-------------------------------------------------------------------------------
!> Purpose : set the 
!!
!! @param[in]  Nxigtr
!! @param[in]  Nxigstr
!! @param[in]  Nxrgtrz
!!
!-------------------------------------------------------------------------------
!
subroutine set_x11_workspace(Nxigtr, Nxigstr, Nxrgtrz)
!     
  use graphmod , only : MXIGTR,MXIGSTR,MXRGTRZ,INITIALIZED

  implicit none
  ! ** Arguments
  !------------------------------------------------------------
  integer, intent(in) :: Nxigtr, Nxigstr, Nxrgtrz

  !  ...set parameters for graphics package
  MXIGTR  = Nxigtr
  MXIGSTR = Nxigstr
  MXRGTRZ = Nxrgtrz
 
  INITIALIZED = .TRUE.

end subroutine set_x11_workspace
