!------------------------------------------------------------------------------
!> Purpose : for an ACTIVE element Mdle, routine finds a big neighbor Neig
!!           across face Iface
!!
!! @param[in]  Mdle  - middle node
!! @param[in]  Iface - face number
!! @param[out] Neig  - equal or bigger neighbor across face
!!                     (0 if no such neighbor exists)
!!
!! @revision : May 12
!------------------------------------------------------------------------------
subroutine big_neig_across_face(Mdle,Iface, Neig)
  use data_structure3D
  implicit none
  !  ...Arguments
  integer, intent(in)  :: Mdle
  integer, intent(in)  :: Iface
  integer, intent(out) :: Neig

  !  ...Locals
  character(len=4) :: etype
  integer,dimension(2) :: neig_ls,iface_ls,orient_ls
  integer,dimension(27) :: nodesl,void
  integer :: mdlf,nrneig,ipos,i,inode,ifound
!------------------------------------------------------------------------------

  !  ...initialze to NO neighbor
  Neig=0

  !  ...determine face node
  call elem_face(Mdle,Iface, mdlf)

  !  ...if the face is refined, NO big neighbor on the other side can exist
  if (NODES(mdlf)%ref_kind.ne.0) return

  !  ...determine neighboring element(s)
  call neig_face(mdlf, nrneig,neig_ls,iface_ls,orient_ls)

  !  ...if only one neighbor, return
  if (Nrneig.ne.2) return

  !  ...locate element other than Mdle on the list
  call locate(Mdle,neig_ls,2, ifound)
  select case(ifound)
  case(0)
    write(*,1)Mdle,Iface,mdlf
1   format(' big_neig_across_face: inconsistency! Mdle,Iface,mdlf = ',3(i12,2x))
    write(*,*)'neig_ls = ',neig_ls(1:2)
    call result
    stop
  case(1) ; ipos=2
  case(2) ; ipos=1
  endselect
  Neig=neig_ls(ipos)

  !  ...check that YOUNGEST neighbor across face is unrefined
  !     [refer to subroutine "neig_face" for details]
  if (NODES(Neig)%ref_kind.ne.0) Neig=0


end subroutine big_neig_across_face
