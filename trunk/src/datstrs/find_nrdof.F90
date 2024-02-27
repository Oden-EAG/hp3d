#if HP3D_DEBUG

!------------------------------------------------------------------------------
!> Purpose : total number of H1,H(curl),H(div),L2 dof for a SINGLE
!!           component
!!
!! @param[out] NrdofH - total number of H1      dof for a single component
!! @param[out] NrdofE - total number of H(curl) dof for a single component
!! @param[out] NrdofV - total number of H(div)  dof for a single component
!! @param[out] NrdofQ - total number of L2      dof for a single component
!!
!> @revision Nov 12
!------------------------------------------------------------------------------
!  REMARK : the number of geometry dof is equal to the total number of H1 dof
!           for a single component
!------------------------------------------------------------------------------
subroutine find_nrdof(NrdofH, NrdofE, NrdofV, NrdofQ)
  !
  use data_structure3D
  !
  implicit none
  integer, intent(out) :: NrdofH, NrdofE, NrdofV, NrdofQ
  !
  integer :: iel, ino
  integer :: mdle, nod, nrnodm, ndofH,ndofE,ndofV,ndofQ
  integer :: nodesl(27), norientl(27), nodm(MAXNODM)
  !------------------------------------------------------------------------------
  !
  !  initialize
  NrdofH=0 ; NrdofE=0 ; NrdofV=0 ; NrdofQ=0
  !
  !  loop over active elements
  do iel=1,NRELES
     mdle = ELEM_ORDER(iel)
     !
     !  nodes of modified element
     call elem_nodes( mdle, nodesl,norientl)
     call logic_nodes(mdle, nodesl,nodm,nrnodm)
     !
     !  loop over nodes of modified element
     do ino=1,nrnodm
        nod=nodm(ino)
        !
        if ( Is_active(nod).and.(NODES(nod)%visit.eq.0) ) then
           !
           !  number of H1,H(curl),H(div),L2 dof for a SINGLE component
           call find_ndof(nod, ndofH,ndofE,ndofV,ndofQ)
           NrdofH = NrdofH + ndofH
           NrdofE = NrdofE + ndofE
           NrdofV = NrdofV + ndofV
           NrdofQ = NrdofQ + ndofQ
           !
           NODES(nod)%visit = 1
        endif
        !
        !  end loop over nodes of modified element
     enddo
     !  end loop over active elements
  enddo
  !
  !  lower visitation flag
  call reset_visit
  !
  !
end subroutine find_nrdof

#endif
