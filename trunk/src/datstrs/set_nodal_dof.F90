!----------------------------------------------------------------------
!> Purpose - routine sets interactively dof for nodes
!!      
!> @data Oct 14
!----------------------------------------------------------------------
! REMARK : routine should eventually be moved to TEST_PROJ, since it is
!          merely a debugging routine
!----------------------------------------------------------------------
!
subroutine set_nodal_dof
!
      use data_structure3D , only : NODES,NRNODS,Is_inactive,find_ndof
      use parameters       , only : ZERO,ZONE
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
        write(*,1000) inod,NODES(inod)%type,ndofH,ndofE,ndofV,ndofQ
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
        if (Is_inactive(inod))  cycle
!               
!       reset all dofs to zero
        if (associated(NODES(inod)%zdofH))  NODES(inod)%zdofH=ZERO
        if (associated(NODES(inod)%zdofE))  NODES(inod)%zdofE=ZERO
        if (associated(NODES(inod)%zdofV))  NODES(inod)%zdofV=ZERO
        if (associated(NODES(inod)%zdofQ))  NODES(inod)%zdofQ=ZERO
!
!       node of interest
        if (inod == nod) then
!
!         select space 
          select case(ispace)
!         H1
          case(1) ; NODES(inod)%zdofH(:,idof)=ZONE
!         H(curl)
          case(2) ; NODES(inod)%zdofE(:,idof)=ZONE
!         H(div)
          case(3) ; NODES(inod)%zdofV(:,idof)=ZONE
!         L2                   
          case(4) ; NODES(inod)%zdofQ(:,idof)=ZONE
          endselect
        endif
!        
      enddo
!
!
endsubroutine set_nodal_dof
