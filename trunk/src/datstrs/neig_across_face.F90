!------------------------------------------------------------------------------
!> Purpose : routine finds YOUNGEST neighbor Neig across face Iface of element
!!           middle node Mdle
!!
!! @param[in]  Mdle  - middle node
!! @param[in]  Iface - face number
!! @param[out] Neig  - neighbor (0 if no neighbor exists, i.e. element is on
!!                     the boundary)
!!
!! @revision : May 12
!------------------------------------------------------------------------------
subroutine neig_across_face(Mdle,Iface, Neig)
  implicit none

  !  ...Arguments
  integer, intent(in)  :: Mdle
  integer, intent(in)  :: Iface
  integer, intent(out) :: Neig

  !  ...Locals
  integer,dimension(2) ::  nlist,nvoid
  integer :: mdlf,nrneig,ipos,ifound
!------------------------------------------------------------------------------

  !  ...initialze to NO neighbor 
  Neig=0

  !  ...determine face node
  call elem_face(Mdle,Iface, mdlf)
 
  !  ...determine neighboring element(s)
  call neig_face(mdlf, nrneig,nlist,nvoid,nvoid)

  !  ...if only one neighbor, return 
  if (nrneig.ne.2) return
 
  !  ...locate element other than Mdle on the list
  call locate(Mdle,nlist,2, ifound)
  select case(ifound)
  case(0)
    write(*,1)Mdle,Iface,mdlf
1   format(' neig_across_face: inconsistency! Mdle,Iface,mdlf = ',3(i12,2x))
    write(*,*)'nlist = ',nlist(1:2)
    call result
    stop
  case(1) ; ipos=2
  case(2) ; ipos=1
  endselect
  Neig=nlist(ipos)
  
end subroutine neig_across_face
