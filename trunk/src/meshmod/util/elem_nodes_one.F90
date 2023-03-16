!----------------------------------------------------------------------------
!> @brief Reconstruct nodal connectivity with given parent connectivity
!!
!! @param[in]  Nfath         - father node
!! @param[in]  Nodesl_fath   - nodal connectivity of father
!! @param[in]  Norientl_fath - orientation of father
!! @param[in]  Ison          - son number of this node
!! @param[out] Nod           - node number
!! @param[out] Nodesl        - nodal connectivity of node
!! @param[out] Norientl      - orientation of node
!
!> @date Feb 2023
!----------------------------------------------------------------------------
subroutine elem_nodes_one(Nfath,Nodesl_fath,Norientl_fath,Ison, &
                          Nod,Nodesl,Norientl)
   use element_data
   use data_structure3D
   use refinements
   implicit none
!
   integer, intent(in)  :: Nfath,Ison
   integer, intent(in)  :: Nodesl_fath(27),Norientl_fath(27)
   integer, intent(out) :: Nod
   integer, intent(out) :: Nodesl(27),Norientl(27)
!
   integer :: npar_refs(27),nson_refs(27),nort_refs(27)
!
   integer :: ntype_fath,ntype_cur,nref_fath
   integer :: nvert_fath,nedge_fath
!
   integer :: kref_face(6)
   integer :: iref, ireff, iref1, iref2, iref3, nort
   integer :: j, jp, nodp, nodpp, is, is1, n_nodes
!
#if DEBUG_MODE
   integer :: iprint
   iprint=0
#endif
!
!----------------------------------------------------------------------------
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,*) '------------------------------------------------'
      write(*,*) 'elem_nodes_one: Begin'
      write(*,*) 'elem_nodes_one: Nfath, Ison = ', Nfath, Ison
      write(*,7001) 'nodes ', Nodesl_fath
      write(*,7001) 'orient', Norientl_fath
 7001 format('elem_nodes_one: Nfath ', a5, ' = ',27(i6))
   endif
#endif
!
!-----------------------------------------------------------
!..initialize output
   Nodesl(1:27) = 0; Norientl(1:27) = 0
!
   ntype_fath = NODES(Nfath)%ntype
   nvert_fath = nvert(ntype_fath)
   nedge_fath = nedge(ntype_fath)
   nref_fath  = NODES(Nfath)%ref_kind
   Nod        = NODES(Nfath)%first_son+Ison-1
   ntype_cur  = NODES(Nod)%ntype
!
!-----------------------------------------------------------
!..one step down the tree reconstructing connectivities
   call find_face_ref_flags(ntype_fath,nref_fath, kref_face)
   call decode_ref(ntype_fath,nref_fath, iref1,iref2,iref3)
!
!..pre-compute info for loop
   call npar_ref_all(ntype_fath,Ison,iref1,iref2,iref3, npar_refs)
   call nson_ref_all(ntype_fath,Ison,iref1,iref2,iref3, nson_refs)
   call nort_ref_all(ntype_fath,Ison,iref1,iref2,iref3, nort_refs)
!
   n_nodes = nvert(ntype_cur)+nedge(ntype_cur)+nface(ntype_cur)+1
   do j=1,n_nodes
!
      jp   = npar_refs(j)
      is   = nson_refs(j)
      nort = nort_refs(j)
!
      if (is.eq.0) then
!
!     ...inheritance rule
         Nodesl(j)   = Nodesl_fath(jp)
         Norientl(j) = Norientl_fath(jp)
!
      else
         nodp = Nodesl_fath(jp)
         select case (TYPE_NOD(jp,ntype_fath))
         case (MEDG)
            call rotate_edge(Norientl_fath(jp),is,nort)
            Nodesl(j) = Son(nodp,is)
         case (MDLT)
            ! local and global
            iref  = kref_face(jp-nvert_fath-nedge_fath)
            ireff = NODES(nodp)%ref_kind
            call rotate_trian(iref,ireff,Norientl_fath(jp),is,nort)
            Nodesl(j) = Son(nodp,is)
         case (MDLQ)
            ! local and global
            iref  = kref_face(jp-nvert_fath-nedge_fath)
            ireff = NODES(nodp)%ref_kind
            call rotate_quad(iref,ireff,Norientl_fath(jp), is,is1,nort)
            Nodesl(j) = Son(nodp,is)
            if (is1.ne.0) then
               nodpp     = Nodesl(j)
               Nodesl(j) = Son(nodpp,is1)
            endif
         case (VERT)
            write(*,*) 'elem_nodes_one; VERTEX CANNOT BE PARENT'
            stop
         case default
            Nodesl(j) = Son(nodp,is)
         end select
         Norientl(j) = nort
      endif
   enddo
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7031) Nod
 7031 format('elem_nodes_one: NODES AND ORIENT FOR mdle = ',i6)
      select case(ntype_cur)
      case(MDLN)
         write(*,7103) Nodesl(1:15)
         write(*,7103) Norientl(1:15)
      case(MDLP)
         write(*,7104) Nodesl(1:21)
         write(*,7104) Norientl(1:21)
      case(MDLD)
         write(*,7105) Nodesl(1:19)
         write(*,7105) Norientl(1:19)
      case(MDLB)
         write(*,7106) Nodesl(1:27)
         write(*,7106) Norientl(1:27)
      end select
   endif
 7103 format(4i6,',', 6i6,',',4i6,',', i6)
 7104 format(6i6,',', 9i6,',',2i6,',',3i6,',',i6)
 7105 format(5i6,',', 8i6,',', i6,',',4i6,',',i6)
 7106 format(8i6,',',12i6,',',6i6,',', i6)
#endif
!
end subroutine elem_nodes_one
