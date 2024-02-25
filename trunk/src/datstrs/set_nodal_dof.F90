#if HP3D_DEBUG

!----------------------------------------------------------------------
!> @brief   routine sets interactively dof for nodes
!> @date    Feb 2023
!----------------------------------------------------------------------
! REMARK : routine should eventually be moved to TEST_PROJ, since it is
!          merely a debugging routine
!----------------------------------------------------------------------
subroutine set_nodal_dof
!
      use data_structure3D , only : NODES,NRNODS,Is_inactive,find_ndof
      use parameters       , only : N_COMS,ZERO,ZONE
      use node_types
!
      implicit none
      integer :: inod,ndofH,ndofE,ndofV,ndofQ,ispace,nod,idof
!
!----------------------------------------------------------------------
!
!     1. Print header
!
      write(*,*) '-- Active nodes --'

!     loop over active nodes and print info
      do inod=1,NRNODS
!
!       skip inactive nodes
        if (Is_inactive(inod))  cycle
!
!       find number of dofs associated to node
        call find_ndof(inod, ndofH,ndofE,ndofV,ndofQ)
!
!       print
        write(*,1000) inod,S_Type(NODES(inod)%ntype),ndofH,ndofE,ndofV,ndofQ
 1000   format(' node = ',i3,' ; type = ',a4,' ; dof H,E,V,Q = ',4(i3,' ; '))
!
      enddo
      write(*,*) ''
!
!     2. Select node and dof
      write(*,*) 'Select space : 1 - H1 ; 2 - H(curl) ; 3 - H(div) ; 4 - L2'
      read( *,*) ispace
!
      write(*,*) 'Select node and dof'
      read( *,*) nod,idof
!
!     3. Reset dofs
!
!     loop over active nodes
      do inod=1,NRNODS
!
!       skip inactive nodes
        if (Is_inactive(inod)) cycle
        if (.not. associated(NODES(inod)%dof)) cycle
!
!       reset all dofs to zero
        if (associated(NODES(inod)%dof%zdofH)) NODES(inod)%dof%zdofH=ZERO
        if (associated(NODES(inod)%dof%zdofE)) NODES(inod)%dof%zdofE=ZERO
        if (associated(NODES(inod)%dof%zdofV)) NODES(inod)%dof%zdofV=ZERO
        if (associated(NODES(inod)%dof%zdofQ)) NODES(inod)%dof%zdofQ=ZERO
!
!       node of interest
        if (inod == nod) then
!
!         select space
          select case(ispace)
!         H1
          case(1) ; NODES(inod)%dof%zdofH(:,idof,N_COMS)=ZONE
!         H(curl)
          case(2) ; NODES(inod)%dof%zdofE(:,idof,N_COMS)=ZONE
!         H(div)
          case(3) ; NODES(inod)%dof%zdofV(:,idof,N_COMS)=ZONE
!         L2
          case(4) ; NODES(inod)%dof%zdofQ(:,idof,N_COMS)=ZONE
          endselect
        endif
!
      enddo
!
end subroutine set_nodal_dof

#endif
