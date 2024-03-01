!-------------------------------------------------------------------------------
!> Purpose : set the
!!
!! @param[in]  Nxigtr
!! @param[in]  Nxigstr
!! @param[in]  Nxrgtrz
!!
!-------------------------------------------------------------------------------
!
#if HP3D_USE_X11

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

  INITIALIZED = .true.

end subroutine set_x11_workspace

#else

subroutine set_x11_workspace(Nxigtr, Nxigstr, Nxrgtrz)
  implicit none
  integer, intent(in) :: Nxigtr, Nxigstr, Nxrgtrz
end subroutine set_x11_workspace

#endif
