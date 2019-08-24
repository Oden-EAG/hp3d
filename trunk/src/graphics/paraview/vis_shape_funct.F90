subroutine vis_shape_funct
!      
      use environment , only : QUIET_MODE
!      
      implicit none
      logical :: mode_save
      integer :: icountH,icountE,icountV,icountQ,idec
!      
!--------------------------------------------------------------------------------
!
!     set quiet mode to TRUE
      mode_save = QUIET_MODE
      QUIET_MODE = .TRUE.
!
      write(*,*)'Select: 0 - All spaces ; 1 - H1 ; 2 - H(curl) ; 3 - H(div) ; 4 - L2'
      read( *,*) idec
!
      select case(idec)
      case(0)
        call vis_shape_h1(   icountH)
        call vis_shape_hcurl(icountE)
        call vis_shape_hdiv( icountV)
        call vis_shape_l2(   icountQ)
!      
        write(*,*)'-- Paraview Enumeration --'
        write(*,8000) 0,                       icountH                        -1
        write(*,8001) icountH,                 icountH+icountE                -1
        write(*,8002) icountH+icountE,         icountH+icountE+icountV        -1
        write(*,8003) icountH+icountE+icountV, icountH+icountE+icountV+icountQ-1
 8000   format(' H1      : ',i3,' - ',i3)
 8001   format(' H(curl) : ',i3,' - ',i3)
 8002   format(' H(div)  : ',i3,' - ',i3)
 8003   format(' L2      : ',i3,' - ',i3)
        write(*,*) ''
!        
      case(1) ; call vis_shape_h1(   icountH)
      case(2) ; call vis_shape_hcurl(icountE)
      case(3) ; call vis_shape_hdiv( icountV)
      case(4) ; call vis_shape_l2(   icountQ)
      endselect
!      
!     reset quiet mode
      QUIET_MODE = mode_save
!
!      
endsubroutine vis_shape_funct
!
!
!
!----------------------------------------------------------------------
!> Purpose - routine dumps H1 shape functions to Paraview
!!      
!> @data Oct 14
!----------------------------------------------------------------------
!
subroutine vis_shape_h1(Icount)
!
      use data_structure3D , only : NODES,NRNODS,Is_inactive,find_ndof
      use parameters       , only : ZERO,ZONE
!
      implicit none
      integer, intent(out) :: Icount
      integer :: inod,ndofH,ndofE,ndofV,ndofQ,idof
!      
!----------------------------------------------------------------------
!
!     Step 1 : reset all dofs to zero
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
      enddo
!
!     Step 2 : activate a single dof a time and export to Paraview
!
      icount=0
!
!     loop over active nodes and print info
      do inod=1,NRNODS
!
!       skip inactive nodes      
        if (Is_inactive(inod))  cycle
!
!       find number of dofs associated to node
        call find_ndof(inod, ndofH,ndofE,ndofV,ndofQ)        
!
!       loop over nodes
        do idof=1,ndofH
!
          icount=icount+1
!          
!         activate
          NODES(inod)%zdofH(:,idof)=ZONE
!
!          write(*,*) 'i,inod,type,idof = ',icount-1,inod,NODES(inod)%type,idof

!         export
          call paraview_driver
!
!         deactivate
          NODES(inod)%zdofH(:,idof)=ZERO
! 
        enddo
      enddo
!
!
endsubroutine vis_shape_h1
!
!
!
!----------------------------------------------------------------------
!> Purpose - routine dumps H(curl) shape functions to Paraview
!!      
!> @data Oct 14
!----------------------------------------------------------------------
!
subroutine vis_shape_hcurl(Icount)
!
      use data_structure3D , only : NODES,NRNODS,Is_inactive,find_ndof
      use parameters       , only : ZERO,ZONE
!
      implicit none
      integer,intent(out) :: Icount
      integer :: inod,ndofH,ndofE,ndofV,ndofQ,idof
!
!----------------------------------------------------------------------
!
!     Step 1 : reset all dofs to zero
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
      enddo
!
!     Step 2 : activate a single dof a time and export to Paraview
!
      icount=0
!      
!     loop over active nodes and print info
      do inod=1,NRNODS
!
!       skip inactive nodes      
        if (Is_inactive(inod))  cycle
!
!       find number of dofs associated to node
        call find_ndof(inod, ndofH,ndofE,ndofV,ndofQ)        
!
!       loop over nodes
        do idof=1,ndofE
!
          icount=icount+1
!
!         activate
          NODES(inod)%zdofE(:,idof)=ZONE
!
!         export
          call paraview_driver
!
!         deactivate
          NODES(inod)%zdofE(:,idof)=ZERO
! 
        enddo
      enddo
!
!
endsubroutine vis_shape_hcurl
!
!
!
!----------------------------------------------------------------------
!> Purpose - routine dumps H(div) shape functions to Paraview
!!      
!> @data Oct 14
!----------------------------------------------------------------------
!
subroutine vis_shape_hdiv(Icount)
!
      use data_structure3D , only : NODES,NRNODS,Is_inactive,find_ndof
      use parameters       , only : ZERO,ZONE
!
      implicit none
      integer,intent(out) :: Icount
      integer :: inod,ndofH,ndofE,ndofV,ndofQ,idof
!
!----------------------------------------------------------------------
!
!     Step 1 : reset all dofs to zero
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
      enddo
!
!     Step 2 : activate a single dof a time and export to Paraview
!
      icount=0
!      
!     loop over active nodes and print info
      do inod=1,NRNODS
!
!       skip inactive nodes      
        if (Is_inactive(inod))  cycle
!
!       find number of dofs associated to node
        call find_ndof(inod, ndofH,ndofE,ndofV,ndofQ)        
!
!       loop over nodes
        do idof=1,ndofV
!
          icount=icount+1
!          
!         activate
          NODES(inod)%zdofV(:,idof)=ZONE
!
!         export
          call paraview_driver
!
!         deactivate
          NODES(inod)%zdofV(:,idof)=ZERO
! 
        enddo
      enddo
!
!
endsubroutine vis_shape_hdiv
!
!
!
!----------------------------------------------------------------------
!> Purpose - routine dumps L2 shape functions to Paraview
!!      
!> @data Oct 14
!----------------------------------------------------------------------
!
subroutine vis_shape_l2(Icount)
!
      use data_structure3D , only : NODES,NRNODS,Is_inactive,find_ndof
      use parameters       , only : ZERO , ZONE
!
      implicit none
      integer, intent(out) :: Icount
      integer :: inod,ndofH,ndofE,ndofV,ndofQ,idof
!
!----------------------------------------------------------------------
!
!     Step 1 : reset all dofs to zero
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
      enddo
!
!     Step 2 : activate a single dof a time and export to Paraview
!
      icount=0
!      
!     loop over active nodes and print info
      do inod=1,NRNODS
!
!       skip inactive nodes      
        if (Is_inactive(inod))  cycle
!
!       find number of dofs associated to node
        call find_ndof(inod, ndofH,ndofE,ndofV,ndofQ)        
!
!       loop over nodes
        do idof=1,ndofQ
!
          icount=icount+1
!          
!         activate
          NODES(inod)%zdofQ(:,idof)=ZONE
!
!         export
          call paraview_driver
!
!         deactivate
          NODES(inod)%zdofQ(:,idof)=ZERO
! 
        enddo
      enddo
!
!
endsubroutine vis_shape_l2
