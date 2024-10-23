#include "typedefs.h"
!-----------------------------------------------------------------------
!
!    routine name       - update_deg
!
!-----------------------------------------------------------------------
!
!    latest revision    - May 2023
!
!    purpose            - routine updates values of solution degrees of freedom
!                         for active vertex, edge and face nodes of an 
!                         element, using the GMP parametrizations and
!                         Dirichlet data
!
!    arguments:
!      in:      Mdle    - (midlle node of an) element
!
!------------------------------------------------------------------------
!
subroutine update_deg(Mdle, visit_flag_comm)
!
   use data_structure3D
   use element_data
   use environment, only: QUIET_MODE
!
   implicit none
!
   integer,intent(in)     :: Mdle
   integer,intent(out)    :: visit_flag_comm(MAXNODS)  
!
   integer                :: ntype
   character(len=6)       :: mdltype
!
!..orientation of element nodes
   integer, dimension(12) :: nedge_orient
   integer, dimension(6)  :: nface_orient
!
!..element nodes
   integer, dimension(27) :: nodesl, norientl
!
!..order of approximation for an element
   integer, dimension(19) :: norder
!
!..reference coordinates for an element
   real(8), dimension(3,8) :: xsub
!
!..solution dof for an element
   VTYPE, dimension(MAXEQNH,MAXbrickH,N_COMS) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE,N_COMS) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV,N_COMS) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ,N_COMS) :: zdofQ
!
!..auxiliary variables
   integer :: iv,ie,ifc,ind,iflag,i,kH,kV,kE,k
   integer :: no,nod,nod1,noda(8),idec
!
!..offsets for geometry and solution dof of higher order nodes
   integer :: nH(18),nE(18),nV(6),ndofH(18),ndofE(18),ndofV(18),ndofQ

#if HP3D_DEBUG
integer :: iprint
   iprint=0
#endif
!-----------------------------------------------------------------------
!
!..skip the element if it has already been updated
   if (NODES(Mdle)%visit.eq.1) then
      return
   endif
!
   call find_elem_type(Mdle, mdltype)
   ntype = NODES(Mdle)%ntype
   call refel(Mdle, iflag,no,xsub)
   call elem_nodes(Mdle, nodesl,norientl)
   zdofH = 0.d0
   zdofE = 0.d0
   zdofV = 0.d0
   zdofQ = 0.d0     
!
!..loop through element vertex nodes and update vertex dofs for dirichlet nodes
   do iv = 1,nvert(ntype)
      nod = nodesl(iv)
      if (.not.associated(NODES(nod)%dof))       cycle
!
      if (.not.associated(NODES(nod)%dof%zdofH)) cycle
!
      if (NODES(nod)%visit .eq. 0) then
            if(is_Dirichlet_attr(nod,CONTIN)) then
               call dhpvert(Mdle,iflag,no,xsub(1:3,iv),NODES(nod)%case, &
                           NODES(nod)%bcond, NODES(nod)%dof%zdofH(:,:,N_COMS))
!
            endif
            NODES(nod)%visit=1
            visit_flag_comm(nod) = 1
      endif
      zdofH(1:MAXEQNH,iv,1:N_COMS) = NODES(nod)%dof%zdofH(1:MAXEQNH,1,1:N_COMS)
!
!..end of loop through element vertices
   enddo
!
!-----------------------------------------------------------------------
!
   call find_orient_from_list(ntype,norientl, nedge_orient, nface_orient)
   call find_order_from_list(ntype,nodesl, norder)
!
!..determine offsets for solution dofs
   kH=nvert(ntype)
   kV=0
   kE=0
   do i=1,nedge(ntype)+nface(ntype)
      nH(i)=kH
      nE(i)=kE
      nod = nodesl(nvert(ntype)+i)
      call ndof_nod(NODES(nod)%ntype,norder(i), ndofH(i),ndofE(i),ndofV(i),ndofQ)
      kH=kH+ndofH(i)
      kE=kE+ndofE(i)
   enddo
!  
   do i=1,nface(ntype)
      nV(i)=kV
      kV=kV+ndofV(nedge(ntype)+i)
   enddo
!
!..loop through element edge nodes
   do ie=1,nedge(ntype)
      ind = nvert(ntype)+ie
      nod = nodesl(ind)
      if (.not.associated(NODES(nod)%dof)) then
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
         cycle
      endif
!
      if (.not.associated(NODES(nod)%dof%coord)) then
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
         cycle
      endif
!
!....skip if no solution dofs
      if((.not.associated(NODES(nod)%dof%zdofH)) .and. (.not.associated(NODES(nod)%dof%zdofE))) then 
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
         cycle
      endif

      if (Is_inactive(nod)) then
         write(*,*) 'update_deg: INCONSISTENCY 1'
         stop 1
      endif
!
!....check consistency, the corresponding vertex nodes should be updated
      call edge_to_vert(ntype,ie, noda(1),noda(2))
      do i=1,2
         nod1 = nodesl(noda(i))
         if (NODES(nod1)%visit.eq.0) then
               write(*,*) 'update_deg: INCONSISTENCY 2'
         endif
      enddo
!
      if (NODES(nod)%visit .eq. 0) then
!
         if (is_Dirichlet_attr(nod,CONTIN)) then
!
!        ...update H1 Dirichlet dof
            if (associated(NODES(nod)%dof%zdofH)) then 
               call dhpedgeH(Mdle,iflag,no,xsub,                     &
                           ntype,NODES(nod)%case,NODES(nod)%bcond,  &
                           nedge_orient,nface_orient,norder,ie,     &
                           zdofH, NODES(nod)%dof%zdofH(:,:,N_COMS))
!
               do k=1,ndofH(ie)
                     zdofH(1:MAXEQNH,nH(ie)+k,1:N_COMS)= NODES(nod)%dof%zdofH(1:MAXEQNH,k,1:N_COMS)
               enddo
!
            endif              
         endif
!
         if (is_Dirichlet_attr(nod,TANGEN)) then
!
!        ...update H(Curl) Dirichlet dof	
            if (associated(NODES(nod)%dof%zdofE)) then
               call dhpedgeE(mdle,iflag,no,xsub,                     &
                              ntype,NODES(nod)%case,NODES(nod)%bcond, &
                              nedge_orient,nface_orient,norder,ie,    &
                              NODES(nod)%dof%zdofE(:,:,N_COMS))
!
               do k=1,ndofE(ie)
                  zdofE(1:MAXEQNE,nE(ie)+k,1:N_COMS)= NODES(nod)%dof%zdofE(1:MAXEQNE,k,1:N_COMS)
               enddo
!
            endif
         endif
!
         NODES(nod)%visit = 1
         visit_flag_comm(nod) = 1
!
      endif
!
!..end of loop through element edges
   enddo
!
!-----------------------------------------------------------------------
!
!..loop through element face nodes
   
   do ifc=1,nface(ntype)
!
      ind = nvert(ntype) + nedge(ntype) + ifc
      nod = nodesl(ind)
!
!...skip if no solution dofs
      if (.not.associated(NODES(nod)%dof)) then
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
         cycle
      endif
      if (.not.associated(NODES(nod)%dof%coord)) then
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
         cycle
      endif
!
      if((.not.associated(NODES(nod)%dof%zdofH)) .and. (.not.associated(NODES(nod)%dof%zdofE)) & 
         .and. (.not.associated(NODES(nod)%dof%zdofV))) then
!
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
         cycle
      endif
      if (Is_inactive(nod)) then
         write(*,*) 'update_deg: INCONSISTENCY 3'
         stop 1
      endif
!
!  ...the face can be updated only if the corresponding vertex and 
!     edge nodes have already been updated
      idec=1
      call face_to_vert(ntype,ifc, noda(1),noda(2),noda(3),noda(4))
      call face_to_edge(ntype,ifc, noda(5),noda(6),noda(7),noda(8))
      do i=1,8
         nod1 = nodesl(noda(i))
         if (NODES(nod1)%visit.eq.0) idec=0
      enddo
      if (idec.eq.0) then
         write(*,*) 'update_deg: INCONSISTENCY 4'
         call result
      endif
!
      if (NODES(nod)%visit .eq. 0) then
!
         if (is_Dirichlet_attr(nod,CONTIN)) then
!        ...update H1 Dirichlet dofs
            if (associated(NODES(nod)%dof%zdofH)) then
!
               call dhpfaceH(Mdle,iflag,no,xsub,                     &
                           ntype,NODES(nod)%case,NODES(nod)%bcond, &
                           nedge_orient,nface_orient,norder,ifc,   &
                           zdofH, NODES(nod)%dof%zdofH(:,:,N_COMS))
!
               do k = 1,ndofH(nedge(ntype)+ifc)
                  zdofH(1:MAXEQNH,nH(nedge(ntype)+ifc)+k,1:N_COMS) = NODES(nod)%dof%zdofH(1:MAXEQNH,k,1:N_COMS)
               enddo
            endif
         endif
   !
         if (is_Dirichlet_attr(nod,TANGEN)) then
!           ...update H(Curl) Dirichlet dof	
               if (associated(NODES(nod)%dof%zdofE)) then
                  call dhpfaceE(Mdle,iflag,no,xsub,                     &
                              ntype,NODES(nod)%case,NODES(nod)%bcond, &
                              nedge_orient,nface_orient,norder,ifc,   &
                              zdofE, NODES(nod)%dof%zdofE(:,:,N_COMS))
!
                  do k = 1,ndofE(nedge(ntype)+ifc)
                     zdofE(1:MAXEQNE,nE(nedge(ntype)+ifc)+k,1:N_COMS) = NODES(nod)%dof%zdofE(1:MAXEQNE,k,1:N_COMS)
                  enddo
!
               endif
         endif
!
         if (is_Dirichlet_attr(nod,NORMAL)) then
!           ...update H(div) Dirichlet dof		
               if(associated(NODES(nod)%dof%zdofV)) then
!
                  call dhpfaceV(mdle,iflag,no,xsub,                     &
                              ntype,NODES(nod)%case,NODES(nod)%bcond, &
                              nedge_orient,nface_orient,norder,ifc,   &
                              NODES(nod)%dof%zdofV(:,:,N_COMS))
!
                  do k = 1,ndofV(nedge(ntype)+ifc)
                     zdofV(1:MAXEQNV,nV(ifc)+k,1:N_COMS) = NODES(nod)%dof%zdofV(1:MAXEQNV,k,1:N_COMS)
                  enddo
!	
               endif
!
         endif
!	
         NODES(nod)%visit = 1
         visit_flag_comm(nod) = 1
      endif
!
!..end of loop through element faces
   enddo
!
!..raise the flag indicating that the element has been updated
   NODES(Mdle)%visit=1
   visit_flag_comm(Mdle) = 1
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
   do i=1,nvert(ntype)+nedge(ntype)+nface(ntype)
      nod = nodesl(i)
      write(*,7100) i,nod,NODES(nod)%visit
7100      format('update_deg: i,nod,NODES(nod)%visit = ', i2,i6,2x,i1)
   enddo
   endif
#endif
!
end subroutine update_deg

    
