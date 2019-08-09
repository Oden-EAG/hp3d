!> Purpose : reconstruct nodal connectivity with given parent connectivity
!! @param[in]  Nfath         - father node
!! @param[in]  Nodesl_fath   - nodal connectivity of father
!! @param[in]  Norientl_fath - orientation of father
!! @param[in]  Ison          - son number of this node
!! @param[out] Nod           - node number
!! @param[out] Nodesl        - nodal connectivity of node
!! @param[out] Norinetl      - orientation of node
!-----------------------------------------------------------
subroutine elem_nodes_one( &
     Nfath, Nodesl_fath, Norientl_fath, Ison, &
     Nod, Nodesl, Norientl)
  use element_data
  use data_structure3D
  use refinements
  implicit none
  ! ** Arguments
  integer,                intent(in)  :: Nfath,Ison
  integer, dimension(27), intent(in)  :: Nodesl_fath,Norientl_fath
  integer,                intent(out) :: Nod
  integer, dimension(27), intent(out) :: Nodesl,Norientl
  ! ** Locals
  type(node) :: fath, cur
  integer, dimension(6) :: kref_face
  integer :: iprint, iref, ireff, iref1, iref2, iref3, nort
  integer :: j, jp, nodp, nodpp, is, is1, n_nodes
  !-----------------------------------------------------------
  iprint = 0
  if (iprint.eq.1) then
     write(*,*) '------------------------------------------------'
     write(*,*) 'elem_nodes_one: Begin'
     write(*,*) 'elem_nodes_one: Nfath, Ison = ', Nfath, Ison
     write(*,7001) 'nodes ', Nodesl_fath
     write(*,7001) 'orient', Norientl_fath
7001 format('elem_nodes_one: Nfath ', a5, ' = ',27(i6))
  endif

  !-----------------------------------------------------------
  ! initialize output
  Nodesl   = 0
  Norientl = 0
  fath     = NODES(Nfath)
  Nod      = fath%first_son+Ison-1
  cur      = NODES(Nod)

  !-----------------------------------------------------------
  ! one step down the tree reconstructing connectivities
  call find_face_ref_flags(fath%type, fath%ref_kind, kref_face)
  call decode_ref(fath%type, fath%ref_kind, iref1, iref2, iref3)
  !
  n_nodes = nvert(cur%type)+nedge(cur%type)+nface(cur%type)+1
  do j=1, n_nodes

     jp   = npar_ref(fath%type, j, Ison, iref1,iref2,iref3)
     is   = nson_ref(fath%type, j, Ison, iref1,iref2,iref3)
     nort = nort_ref(fath%type, j, Ison, iref1,iref2,iref3)

     if (is.eq.0) then
        !
        ! inheritance rule
        Nodesl(j)   = Nodesl_fath(jp)
        Norientl(j) = Norientl_fath(jp)
        !
     else
        nodp = Nodesl_fath(jp)
        select case (Type_nod(fath%type, jp))
        case ('vert')
           write(*,*) 'elem_nodes_one; VERTEX CANNOT BE PARENT'
           stop 1
        case ('medg')
           call rotate_edge(Norientl_fath(jp),is,nort)
           !Nodesl(j) = NODES( nodp )%sons(is)
           Nodesl(j) = Son(nodp,is)
        case ('mdlt')
           ! local and global
           iref  = kref_face(jp-nvert(fath%type)-nedge(fath%type))
           ireff = NODES(nodp)%ref_kind
           call rotate_trian(iref,ireff,Norientl_fath(jp),is,nort)
           !Nodesl(j) = NODES( nodp )%sons(is)
           Nodesl(j) = Son(nodp,is)
        case('mdlq')
           ! local and global
           iref  = kref_face(jp-nvert(fath%type)-nedge(fath%type))
           ireff = NODES(nodp)%ref_kind
           call rotate_quad(iref,ireff,Norientl_fath(jp), is,is1,nort)
           !Nodesl(j) = NODES(nodp)%sons(is)
           Nodesl(j) = Son(nodp,is)
           if (is1.ne.0) then
              nodpp     = Nodesl(j)
              !Nodesl(j) = NODES(nodpp)%sons(is1)
              Nodesl(j) = Son(nodpp,is1)
           endif
        case default
           !Nodesl(j) = NODES(nodp)%sons(is)
           Nodesl(j) = Son(nodp,is)
        end select
        Norientl(j) = nort
     endif
  enddo
  !
  if (iprint.eq.1) then
     write(*,7031) Nod
7031 format('elem_nodes: NODES AND ORIENT FOR mdle = ',i6)
     select case(cur%type)
     case('mdln')
        write(*,7103) Nodesl(1:15)
        write(*,7103) Norientl(1:15)
     case('mdlp')
        write(*,7104) Nodesl(1:21)
        write(*,7104) Norientl(1:21)
     case('mdld')
        write(*,7105) Nodesl(1:19)
        write(*,7105) Norientl(1:19)
     case('mdlb')
        write(*,7106) Nodesl(1:27)
        write(*,7106) Norientl(1:27)
     end select
  endif
7103 format(4i6,2x,6i6,2x,4i6,2x,i6)
7104 format(6i6,2x,9i6,2x,2i6,2x,3i6,2x,i6)
7105 format(5i6,2x,8i6,2x,i6,2x,4i6,2x,i6)
7106 format(8i6,2x,12i6,2x,6i6,2x,i6)
end subroutine elem_nodes_one

