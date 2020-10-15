!> Purpose : calculate modified stiffness matrix and load vector
!! @param[in]  Mdle   - middle node of an element
!! @param[in]  Idec   - 1 nodes, 2 matrix
!!
!! @param[out] Nrdofs - # of element local dof for each physics
!! @param[out] Nrdofm - # of modified element dof in the expanded
!! @param[out] Nrdofc - # of modified element dof after compression
!!
!! @param[out] Nodm   - actual nodes in the order
!!
!! @param[out] NdofmH - the number of H1 dof
!! @param[out] NdofmE - the number of H(curl) dof
!! @param[out] NdofmV - the number of H(div) dof
!! @param[out] NdofmQ - the number of L2 dof
!!
!! @param[out] Nrnodm - number of the modified element nodes
!!
!! @param[out] Bload  - 1D array containing the modified load vector
!! @param[out] Astif  - 1D array containing the modified stiffness matrix
!!
subroutine celem( &
     Mdle,Idec, &
     Nrdofs,Nrdofm,Nrdofc, &
     Nodm, &
     NdofmH,NdofmE,NdofmV,NdofmQ, &
     Nrnodm, &
     Bload,Astif)
  use physics
  use parameters
  !
  implicit none
  !
  integer,                         intent(in)  :: Mdle,Idec
  integer,    dimension(NR_PHYSA), intent(out) :: Nrdofs
  integer,                         intent(out) :: Nrdofm,Nrdofc
  integer,    dimension(MAXNODM),  intent(out) :: Nodm
  integer,    dimension(MAXNODM),  intent(out) :: NdofmH,NdofmE,NdofmV,NdofmQ
  integer,                         intent(out) :: Nrnodm
  real*8,                          intent(out) :: Bload(*),Astif(*)
  !
  ! temporary regular system celem for single physics
  call celem_system(Mdle,Idec, &
       Nrdofs,Nrdofm,Nrdofc, &
       Nodm, &
       NdofmH,NdofmE,NdofmV,NdofmQ, &
       Nrnodm, &
       Bload,Astif)
  !
end subroutine celem