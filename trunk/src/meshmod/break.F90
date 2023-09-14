!-------------------------------------------------------------------------
!> @brief     routine breaks an element middle node according to a given
!!            refinement flag. Compatibility is assumed across the faces
!!            and the edges.
!!
!> @param[in] Mdle - middle node
!> @param[in] Kref - refinement flag
!!
!> @date      Feb 2023
!-------------------------------------------------------------------------
subroutine break(Mdle,Kref)
!
   use data_structure3D
   use element_data
   use refinements
!
   implicit none
!
   integer, intent(in) :: Mdle, Kref
!
   integer :: ntype
   integer :: nodesl(27),norientl(27)
   integer :: kreff(6)
   integer :: krefe(12)
   integer :: i, is, iprint, iface, ipass
   integer :: kref_face, kref_edge, nod, nrsons
   logical :: iact
!
!-------------------------------------------------------------------------
!
   iprint = 0
!
!..nodal connectivities
   call elem_nodes(Mdle, nodesl,norientl)
!
!..determine refinements for the element faces and edges
   ntype = NODES(Mdle)%ntype
   call find_face_ref_flags(ntype,Kref, kreff)
   call find_edge_ref_flags(ntype,Kref, krefe)
!
#if DEBUG_MODE
   if (iprint.eq.2) then
      call elem_show(Mdle,ntype,nodesl,norientl)
 7000 format('break: Mdle = ',i7)
      write(*,7011) Kref
 7011 format('       Kref  = ',i3)
      write(*,7012) kreff(1:NFACE(ntype))
 7012 format('       kreff = ',6i3)
      write(*,7013) krefe(1:NEDGE(ntype))
 7013 format('       krefe = ',12i2)
      call pause
   endif
#endif
!
!=========================================================================
! E D G E S                                                              |
!=========================================================================
!
!..loop over edges and generate INACTIVE son nodes
   do i=1,NEDGE(ntype)
      nod=nodesl(NVERT(ntype)+i)
      if (iprint .eq. 2) write(*,*) 'edge nod = ', nod
!  ...check if edge is unrefined and needs to be refined
      if ((krefe(i).ne.0) .and. (is_leaf(nod))) then
!     ...break edge
         kref_edge = 1; iact = .false.
         call nodbreak(nod,kref_edge,iact)
      endif
   enddo
!
!=========================================================================
! F A C E S                                                              |
!=========================================================================
!
!..loop over faces and generate ACTIVE/INACTIVE son nodes
   do i=1,NFACE(ntype)
      iface=NVERT(ntype)+NEDGE(ntype)+i
      nod=nodesl(iface)
      if (iprint .eq. 2) write(*,*) 'face nod = ', nod
!
!  ...face should be refined
      if (kreff(i).ne.0) then
!
!     ...modify refinement according to orientation
         call change_ref_flag('l2g',NODES(nod)%ntype,kreff(i), &
                              norientl(iface), kref_face)
!
!     ...check
         call check_ref(NODES(nod)%ntype,NODES(nod)%ref_kind, &
                        kref_face, ipass)
         if (ipass.eq.0) then
            write(*,7002) mdle,i,NODES(nod)%ref_kind,kref_face
7002        format('break: mdle,i,NODES(nod)%ref_kind,kref_face = ',i6,2x,3(2x,i2))
            call result
            return
         endif
!
!     ...compare desired refinement "kref_face" to existing one "%ref_kind"
         kref_face=kref_face-NODES(nod)%ref_kind
!
!     ...by default, generate INACTIVE son nodes
         iact=.false.
!
!     ...existing refinement coincides with desired refinement
         if (kref_face.eq.0) then
!        ...activate original constrained faces (sons) that are now unconstrained
!        ...face activation needed for 'refine' routine
            call activate_sons(nod)
!
!     ...existing refinement does NOT coincide with desired refinement
         else
            select case(NODES(nod)%ntype)
!
!        ...triangular face
            case(MDLT)
               call nodbreak(nod,kref_face,iact)
!
!        ...quadrilateral face
            case(MDLQ)
               select case(NODES(nod)%ref_kind)
!           ...UNREFINED face : just break
               case(0)
                  call nodbreak(nod,kref_face,iact)
!           ...ANISOTROPICALLY REFINED face : break son quads and son edge
               case(10,01)
!              ...activate original constrained faces (sons) that are now unconstrained
!              ...face activation needed for 'refine' routine
                  call activate_sons(nod)
!              ...loop over son quads and generate INACTIVE son nodes
                  do is=1,2
                     call nodbreak(Son(nod,is),kref_face,iact)
                  enddo
!              ...break edge node
                  kref_edge = 1; call nodbreak(Son(nod,3),kref_edge,iact)
               end select
            end select
         endif
      endif
!..end loop over faces
   enddo
!
!=========================================================================
! M I D D L E                                                            |
!=========================================================================
!
!..break middle node
   iact = .false.
   call nodbreak(Mdle,Kref,iact)
!
   call activate_sons(Mdle)
!
!..update the number of active elements
   call nr_mdle_sons(ntype,Kref, nrsons)
   NRELES=NRELES+nrsons-1
!
#if DEBUG_MODE
   if (iprint.ge.1) then
      write(*,7100) Mdle, S_Type(NODES(Mdle)%ntype), Kref
 7100 format('break: Mdle ',i7,' ',a4, ' HAS BEEN BROKEN WITH Kref = ', i3)
   endif
#endif
!
end subroutine break
!
!
!------------------------------------------------------------------------
!> @brief     activate son nodes
!!            remark: in distributed mesh, this routines relies upon that
!!                    the Nod's subdomain was set correctly previously
!!
!> @param[in] Nod - node number
!> @date      Feb 2023
!------------------------------------------------------------------------
subroutine activate_sons(Nod)
   use data_structure3D
   implicit none
   integer, intent(in) :: Nod
   integer             :: nrsons,subd,i
!
   call find_nsons(Nod, nrsons)
   call get_subd(Nod, subd)
!..loop over sons
   do i=1,nrsons
!  ...set subdomain for sons
      call set_subd(Son(Nod,i),subd)
!  ...activate sons
      if (Is_inactive(Son(Nod,i))) then
         call activate(Son(Nod,i), NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ)
      endif
   enddo
!..deactivate nod
   call deactivate(Nod, NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ)
!
end subroutine activate_sons
