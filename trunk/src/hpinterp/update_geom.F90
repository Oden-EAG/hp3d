!-----------------------------------------------------------------------
!
!    routine name       - update_geom
!
!-----------------------------------------------------------------------
!
!    latest revision    - May 2023
!
!    purpose            - routine updates values of geometry degrees
!                         for active vertex, edge and face nodes of an 
!                         element, using the GMP parametrizations
!
!    arguments:
!      in:      Mdle    - (midlle node of an) element
!
!------------------------------------------------------------------------
!
subroutine update_geom(Mdle, visit_flag_comm)
!
   use data_structure3D
   use element_data
   use environment, only: QUIET_MODE
!
   implicit none
!
   integer,intent(in) :: Mdle
   integer,intent(out):: visit_flag_comm(MAXNODS)      
!
   integer            :: ntype
   character(len=6)   :: mdltype
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
!..geometry dofs for an element
   real(8), dimension(3,MAXbrickH) :: xnod
!
!..auxiliary variables
   integer :: iv,ie,ifc,ind,iflag,i,k
   integer :: no,nod,nod1,noda(8),idec,iprint
!
!..offsets for geometry dof of higher order nodes
   integer na(18),ndofH(18),ndofE,ndofV,ndofQ
!
#if HP3D_DEBUG
integer :: iprint
   iprint=0
#endif
!-----------------------------------------------------------------------
!
!
!..skip the element if it has alreday been updated
   if (NODES(Mdle)%visit.eq.1) return
!
   call find_elem_type(Mdle, mdltype)
   ntype = NODES(Mdle)%ntype
   call refel(Mdle, iflag,no,xsub)
   call elem_nodes(Mdle, nodesl,norientl)
   xnod = 0.d0
!
!..loop through the element vertex nodes
   do iv=1,nvert(ntype)
      nod = nodesl(iv)
      if (.not.associated(NODES(nod)%dof)) cycle
      if (.not.associated(NODES(nod)%dof%coord)) cycle
      if (NODES(nod)%visit.eq.0) then
         call hpvert(iflag,no,xsub(1:3,iv), NODES(nod)%dof%coord)
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
      endif
      xnod(1:3,iv) = NODES(nod)%dof%coord(1:3,1)
!
!..end of loop through element vertices
   enddo
!
!-----------------------------------------------------------------------
!
   call find_orient_from_list(ntype,norientl, nedge_orient, nface_orient)
   call find_order_from_list(ntype,nodesl, norder)
!
!..determine offsets for geometry dof
   k=nvert(ntype)
   do i=1,nedge(ntype)+nface(ntype)
      na(i)=k
      nod = nodesl(nvert(ntype)+i)
      call ndof_nod(NODES(nod)%ntype,norder(i), ndofH(i),ndofE,ndofV,ndofQ)
      k=k+ndofH(i)
   enddo
!
!..loop through element edge nodes
   do ie=1,nedge(ntype)
      ind = nvert(ntype)+ie
      nod = nodesl(ind)
!
!  ...skip if no gdof
      if (.not.associated(NODES(nod)%dof)) then
         NODES(nod)%visit=1;
         visit_flag_comm(nod) = 1 
         cycle
      endif
      if (.not.associated(NODES(nod)%dof%coord)) then
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
         cycle
      endif
      if (Is_inactive(nod)) then
         write(*,*) 'update_geom: INCONSISTENCY 1'
         stop 1
      endif
!
!  ...check consistency, the corresponding vertex nodes should be updated
      call edge_to_vert(ntype,ie, noda(1),noda(2))
      do i=1,2
         nod1 = nodesl(noda(i))
         if (NODES(nod1)%visit.eq.0) then
         write(*,*) 'update_geom: INCONSISTENCY 2'
         endif
      enddo
!
      if (NODES(nod)%visit.eq.0) then
         if (mdltype.ne.'Linear') then
            call hpedge(mdle,iflag,no,xsub,ntype,             &
                        nedge_orient,nface_orient,norder,ie,  &
                        xnod,NODES(nod)%dof%coord)
         else
            NODES(nod)%dof%coord = 0.d0
         endif
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
      endif
      do k=1,ndofH(ie)
         xnod(1:3,na(ie)+k) = NODES(nod)%dof%coord(1:3,k)
      enddo
!
!..end of loop through element edges
   enddo
!
!-----------------------------------------------------------------------
!
!  ...loop through element face nodes
   do ifc=1,nface(ntype)
      ind = nvert(ntype)+nedge(ntype)+ifc
      nod = nodesl(ind)
!
!  ...skip if no gdof
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
      if (Is_inactive(nod)) then
         write(*,*) 'update_geom: INCONSISTENCY 3'
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
         write(*,*) 'update_geom: INCONSISTENCY 4'
         call result
      endif
!
      if (NODES(nod)%visit.eq.0) then
         if (mdltype.ne.'Linear') then
            call hpface(mdle,iflag,no,xsub,ntype,             &
                        nedge_orient,nface_orient,norder,ifc, &
                        xnod,NODES(nod)%dof%coord)
         else
            NODES(nod)%dof%coord = 0.d0
         endif
         NODES(nod)%visit=1
         visit_flag_comm(nod) = 1
      endif
      do k=1,ndofH(nedge(ntype)+ifc)
         xnod(1:3,na(nedge(ntype)+ifc)+k) = NODES(nod)%dof%coord(1:3,k)
      enddo
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
7100      format('update_geom: i,nod,NODES(nod)%visit = ', i2,i6,2x,i1)
      enddo
      call pause
   endif
#endif
!
end subroutine update_geom