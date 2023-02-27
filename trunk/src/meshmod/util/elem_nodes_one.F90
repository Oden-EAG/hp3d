!----------------------------------------------------------------------------
!> @brief reconstruct nodal connectivity with given parent connectivity
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
subroutine elem_nodes_one( &
     Nfath,Nodesl_fath,Norientl_fath,Ison, &
     Nod,Nodesl,Norientl)
   use element_data
   use data_structure3D
   use refinements
   implicit none
!
   integer,                intent(in)  :: Nfath,Ison
   integer, dimension(27), intent(in)  :: Nodesl_fath,Norientl_fath
   integer,                intent(out) :: Nod
   integer, dimension(27), intent(out) :: Nodesl,Norientl
!
   type(node) :: fath, cur
   integer, dimension(6) :: kref_face
   integer :: iref, ireff, iref1, iref2, iref3, nort
   integer :: j, jp, nodp, nodpp, is, is1, n_nodes
!
#if DEBUG_MODE
   integer :: iprint = 0
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
   Nodesl   = 0
   Norientl = 0
   fath     = NODES(Nfath) ! TODO: MAKES COPY OF NODE (use pointer?)
   Nod      = fath%first_son+Ison-1
   cur      = NODES(Nod) ! TODO: MAKES COPY OF NODE (use pointer?)
!
!-----------------------------------------------------------
!..one step down the tree reconstructing connectivities
   call find_face_ref_flags(fath%ntype,fath%ref_kind, kref_face)
   call decode_ref(fath%ntype,fath%ref_kind, iref1,iref2,iref3)
!
   n_nodes = nvert(cur%ntype)+nedge(cur%ntype)+nface(cur%ntype)+1
   do j=1, n_nodes
!
      jp   = npar_ref(fath%ntype,j,Ison,iref1,iref2,iref3)
      is   = nson_ref(fath%ntype,j,Ison,iref1,iref2,iref3)
      nort = nort_ref(fath%ntype,j,Ison,iref1,iref2,iref3)
!
      if (is.eq.0) then
!
!     ...inheritance rule
         Nodesl(j)   = Nodesl_fath(jp)
         Norientl(j) = Norientl_fath(jp)
!
      else
         nodp = Nodesl_fath(jp)
         select case (Type_nod(fath%ntype, jp))
         case (MEDG)
            call rotate_edge(Norientl_fath(jp),is,nort)
            Nodesl(j) = Son(nodp,is)
         case (MDLT)
            ! local and global
            iref  = kref_face(jp-nvert(fath%ntype)-nedge(fath%ntype))
            ireff = NODES(nodp)%ref_kind
            call rotate_trian(iref,ireff,Norientl_fath(jp),is,nort)
            Nodesl(j) = Son(nodp,is)
         case (MDLQ)
            ! local and global
            iref  = kref_face(jp-nvert(fath%ntype)-nedge(fath%ntype))
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
      select case(cur%ntype)
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

