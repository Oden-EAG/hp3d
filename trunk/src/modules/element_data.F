c----------------------------------------------------------------------
c> Purpose : element data for elements of different type
c
c> @revision Nov 12
c----------------------------------------------------------------------
c
      module element_data
c
!       save
c
c----------------------------------------------------------------------
c  MASTER ELEMENTS VERTEX COORDINATES (vertex enumeration)            |
c----------------------------------------------------------------------
      real*8, parameter, dimension(2,3) :: TRIAN_COORD =
     .reshape(
     .(/0.d0,0.d0, 1.d0,0.d0, 0.d0,1.d0/)
     .,(/2,3/))
c
      real*8, parameter, dimension(2,4) :: QUADR_COORD =
     .reshape(
     .(/0.d0,0.d0, 1.d0,0.d0, 1.d0,1.d0, 0.d0,1.d0/)
     .,(/2,4/))
c
      real*8, parameter, dimension(3,6) :: PRISM_COORD =
     .reshape(
     .(/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0,
     .  0.d0,0.d0,1.d0, 1.d0,0.d0,1.d0, 0.d0,1.d0,1.d0/)
     .,(/3,6/))
c
      real*8, parameter, dimension(3,8) :: BRICK_COORD =
     .reshape(
     .(/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 1.d0,1.d0,0.d0, 0.d0,1.d0,0.d0,
     .  0.d0,0.d0,1.d0, 1.d0,0.d0,1.d0, 1.d0,1.d0,1.d0, 0.d0,1.d0,1.d0/)
     .,(/3,8/))
c
      real*8, parameter, dimension(3,4) :: TETRA_COORD =
     .reshape(
     .(/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0,
     .  0.d0,0.d0,1.d0/)
     .,(/3,4/))
c
      real*8, parameter, dimension(3,5) :: PYRAM_COORD =
     .reshape(
     .(/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 1.d0,1.d0,0.d0, 0.d0,1.d0,0.d0,
     .  0.d0,0.d0,1.d0/)
     .,(/3,5/))
c
c----------------------------------------------------------------------
c  EDGE_2_VERT CONNECTIVITIES (edge enumeration & orientation)        |
c----------------------------------------------------------------------
      integer, parameter, dimension(2,3) :: TRIAN_EDGE_TO_VERT =
     .reshape(
     .(/1,2, 2,3, 1,3/)
     .,(/2,3/))
c
      integer, parameter, dimension(2,4) :: QUADR_EDGE_TO_VERT =
     .reshape(
     .(/1,2, 2,3, 4,3, 1,4/)
     .,(/2,4/))
c
      integer, parameter, dimension(2,9) :: PRISM_EDGE_TO_VERT =
     .reshape(
     .(/1,2, 2,3, 1,3, 4,5, 5,6, 4,6, 1,4, 2,5, 3,6/)
     .,(/2,9/))
c
      integer, parameter, dimension(2,12) :: BRICK_EDGE_TO_VERT =
     .reshape(
     .(/1,2, 2,3, 4,3, 1,4, 5,6, 6,7, 8,7, 5,8, 1,5, 2,6, 3,7, 4,8/)
     .,(/2,12/))
c
      integer, parameter, dimension(2,6) :: TETRA_EDGE_TO_VERT =
     .reshape(
     .(/1,2, 2,3, 1,3, 1,4, 2,4, 3,4/)
     .,(/2,6/))
c
      integer, parameter, dimension(2,8) :: PYRAM_EDGE_TO_VERT =
     .reshape(
     .(/1,2, 2,3, 4,3, 1,4, 1,5, 2,5, 3,5, 4,5/)
     .,(/2,8/))
c
c----------------------------------------------------------------------
c  FACE_2_VERT CONNECTIVITIES (face enumeration & orientation)        |
c----------------------------------------------------------------------
      integer, parameter, dimension(4,5) :: PRISM_FACE_TO_VERT =
     .reshape(
     .(/1,2,3,1, 4,5,6,4, 1,2,5,4, 2,3,6,5, 1,3,6,4/)
     .,(/4,5/))
c
      integer, parameter, dimension(4,6) :: BRICK_FACE_TO_VERT =
     .reshape(
     .(/1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,4,8,5/)
     .,(/4,6/))
c
      integer, parameter, dimension(3,4) :: TETRA_FACE_TO_VERT =
     .reshape(
     .(/1,2,3, 1,2,4, 2,3,4, 1,3,4/)
     .,(/3,4/))
c
      integer, parameter, dimension(4,5) :: PYRAM_FACE_TO_VERT =
     .reshape(
     .(/1,2,3,4, 1,2,5,1, 2,3,5,2, 4,3,5,4, 1,4,5,1/)
     .,(/4,5/))
c
c----------------------------------------------------------------------
c  FACE_2_EDGE CONNECTIVITIES (redundant!)                            |
c----------------------------------------------------------------------
      integer, parameter, dimension(4,5) :: PRISM_FACE_TO_EDGE =
     .reshape(
     .(/1,2,3,1, 4,5,6,4, 1,8,4,7, 2,9,5,8, 3,9,6,7/)
     .,(/4,5/))
c
      integer, parameter, dimension(4,6) :: BRICK_FACE_TO_EDGE =
     .reshape(
     .(/1,2,3,4, 5,6,7,8, 1,10,5,9, 2,11,6,10, 3,11,7,12, 4,12,8,9/)
     .,(/4,6/))
c
      integer, parameter, dimension(3,4) :: TETRA_FACE_TO_EDGE =
     .reshape(
     .(/1,2,3, 1,5,4, 2,6,5, 3,6,4/)
     .,(/3,4/))
c
      integer, parameter, dimension(4,5) :: PYRAM_FACE_TO_EDGE =
     .reshape(
     .(/1,2,3,4, 1,6,5,1, 2,7,6,2, 3,7,8,3, 4,8,5,4/)
     .,(/4,5/))
c
c----------------------------------------------------------------------
c  VERT_2_EDGE CONNECTIVITIES (redundant!)                            |
c----------------------------------------------------------------------
      integer, parameter, dimension(3,4) :: TETRA_VERT_TO_EDGE =
     .reshape(
     .(/1,3,4, 1,2,5, 2,3,6, 4,5,6/)
     .,(/3,4/))
c
      integer, parameter, dimension(3,8) :: BRICK_VERT_TO_EDGE =
     .reshape(
     .(/4,1,9, 1,2,10, 2,3,11, 3,4,12, 8,5,9, 5,6,10, 6,7,11, 7,8,12/)
     .,(/3,8/))
c
      integer, parameter, dimension(3,6) :: PRISM_VERT_TO_EDGE =
     .reshape(
     .(/1,3,7, 1,2,8, 2,3,9, 4,6,7, 4,5,8, 5,6,9/)
     .,(/3,6/))
c
c----------------------------------------------------------------------
c  VERT_2_FACE CONNECTIVITIES (redundant!)                            |
c----------------------------------------------------------------------
      integer, parameter, dimension(3,6) :: PRISM_VERT_TO_FACE =
     .reshape(
     .(/1,3,5, 1,3,4, 1,4,5, 2,3,5, 2,3,4, 2,4,5/)
     .,(/3,6/))
c
      integer, parameter, dimension(3,8) :: BRICK_VERT_TO_FACE =
     .reshape(
     .(/1,6,3, 1,3,4, 1,4,5, 1,5,6, 2,6,3, 2,3,4, 2,4,5, 2,5,6/)
     .,(/3,8/))
c
      integer, parameter, dimension(3,4) :: TETRA_VERT_TO_FACE =
     .reshape(
     .(/1,2,4, 1,2,3, 1,3,4, 2,3,4/)
     .,(/3,4/))
c
      integer, parameter, dimension(4,5) :: PYRAM_VERT_TO_FACE =
     .reshape(
     .(/1,2,5,0, 1,2,3,0, 1,3,4,0, 1,4,5,0, 2,3,4,5/)
     .,(/4,5/))
c
c----------------------------------------------------------------------
c  EDGE_2_FACE CONNECTIVITIES (redundant!)                            |
c----------------------------------------------------------------------
      integer, parameter, dimension(2,9) :: PRISM_EDGE_TO_FACE =
     .reshape(
     .(/1,3, 1,4, 1,5, 2,3, 2,4, 2,5, 3,5, 3,4, 4,5/)
     .,(/2,9/))
c
      integer, parameter, dimension(2,12) :: BRICK_EDGE_TO_FACE =
     .reshape(
     .(/1,3, 1,4, 1,5, 1,6, 2,3, 2,4, 2,5, 2,6, 3,6, 3,4, 4,5, 5,6/)
     .,(/2,12/))
c
      integer, parameter, dimension(2,6) :: TETRA_EDGE_TO_FACE =
     .reshape(
     .(/1,2, 1,3, 1,4, 2,4, 2,3, 3,4/)
     .,(/2,6/))
c
      integer, parameter, dimension(2,8) :: PYRAM_EDGE_TO_FACE =
     .reshape(
     .(/1,2, 1,3, 1,4, 1,5, 2,5, 2,3, 3,4, 4,5/)
     .,(/2,8/))
c
c----------------------------------------------------------------------
c  LOCAL_2_GLOBAL ORIENTATIONS                                        |
c----------------------------------------------------------------------
      integer, parameter, dimension(2,0:1) :: EDGE_L2G =
     .reshape(
     .(/1,2, 2,1/)
     .,(/2,2/))
c
      integer, parameter, dimension(3,0:5) :: TRIAN_L2G =
     .reshape(
     .(/1,2,3, 2,3,1, 3,1,2, 1,3,2, 2,1,3, 3,2,1/)
     .,(/3,6/))
c
      integer, parameter, dimension(3,0:5) :: TRIAN_ORIENT_L2G =
     .reshape(
     .(/0,0,0, 0,1,1, 1,0,1, 0,1,0, 1,0,0, 1,1,1/)
     .,(/3,6/))

      integer, parameter, dimension(3,0:5) :: TRIAN_L2G_EDGE =
     .reshape(
     .(/1,2,3, 2,3,1, 3,1,2, 3,2,1, 1,3,2, 2,1,3/)
     .,(/3,6/))
c
      integer, parameter, dimension(4,0:7) :: QUADR_L2G_EDGE =
     .reshape(
     .(/1,2,3,4, 2,3,4,1, 3,4,1,2, 4,1,2,3,
     .  4,3,2,1, 1,4,3,2, 2,1,4,3, 3,2,1,4/)
     .,(/4,8/))


      integer, parameter, dimension(4,0:7) :: QUADR_L2G =
     .reshape(
     .(/1,2,3,4, 2,3,4,1, 3,4,1,2, 4,1,2,3,
     .  1,4,3,2, 2,1,4,3, 3,2,1,4, 4,3,2,1/)
     .,(/4,8/))
c
      integer, parameter, dimension(4,0:7) :: QUADR_ORIENT_L2G =
     .reshape(
     .(/0,0,0,0, 0,1,0,1, 1,1,1,1, 1,0,1,0,
     .  0,0,0,0, 1,0,1,0, 1,1,1,1, 0,1,0,1/)
     .,(/4,8/))
c
c----------------------------------------------------------------------
c  QUAD_2_AXES ORIENTATIONS                                           |
c                                                                     |
c  NFAXES(3,j)=0,[1] : x and y axes have not been [HAVE BEEN] swapped |
c                      for j-th orientation                           |
c                                                                     |
c  NFAXES(1,j)=0,[1] : consistent [INCONSISTENT] orientation for      |
c                      horizontal axis, for j-th orientation          |
c                                                                     |
c  NFAXES(2,j)=0,[1] : consistent [INCONSISTENT] orientation for      |
c                      vertical axis, for j-th orientation            |
c----------------------------------------------------------------------
      integer, parameter, dimension(3,0:7) :: NFAXES =
     .reshape(
     .(/0,0,0, 1,0,1, 1,1,0, 0,1,1, 0,0,1, 1,0,0, 1,1,1, 0,1,0/)
     .,(/3,8/))
c
c----------------------------------------------------------------------
c  ...defines vertex shape functions; for vertex iv,
c     IJKV(ixi,iv) = number of 1D linear shape function in Xi(ixi)
c
c  ...defines edge shape functions; for edge ie
c     IXIEDGE(ie) = edge coordinate
c     IBLENDE(2,ie) = blending directions
c     NBLENDE(2,ie) = numbers of blending linear shape function
c
c  ...defines face shape functions; for face if
c     IXIFACE(2,if) = face coordinates
c     IBLENDF(if) = blending directions
c     NBLENDF(if) = numbers of blending linear shape function
c
      integer, parameter, dimension(3,8) :: IJKV =
     .reshape(
     .(/1,1,1, 2,1,1, 2,2,1, 1,2,1,
     .  1,1,2, 2,1,2, 2,2,2, 1,2,2/)
     .,(/3,8/))
c
      integer, parameter, dimension(12) :: IXIEDGE =
     .(/1,2,1,2, 1,2,1,2, 3,3,3,3/)
c
      integer, parameter, dimension(2,12) :: IBLENDE =
     .reshape(
     .(/2,3, 1,3, 2,3, 1,3, 2,3, 1,3,
     .  2,3, 1,3, 1,2, 1,2, 1,2, 1,2/)
     .,(/2,12/))
c
      integer, parameter, dimension(2,12) :: NBLENDE =
     .reshape(
     .(/1,1, 2,1, 2,1, 1,1, 1,2, 2,2,
     .  2,2, 1,2, 1,1, 2,1, 2,2, 1,2/)
     .,(/2,12/))
c
      integer, parameter, dimension(2,6) :: IXIFACE =
     .reshape(
     .(/1,2, 1,2, 1,3, 2,3, 1,3, 2,3/)
     .,(/2,6/))
c
      integer, parameter, dimension(6) :: IBLENDF =
     .(/3,3,2,1,2,1/)
c
      integer, parameter, dimension(6) :: NBLENDF =
     . (/1,2,1,2,2,1/)
c
c
c
      contains
c
      function Ivtk_type(Type)
      character(len=4) :: Type
c
      select case(Type)
      case ('vert'); Ivtk_type = 1
      case ('medg'); Ivtk_type = 2
      case ('mdlt'); Ivtk_type = 5
      case ('mdlq'); Ivtk_type = 9
      case ('mdln'); Ivtk_type = 10
      case ('mdlb'); Ivtk_type = 12
      case ('mdlp'); Ivtk_type = 13
      case ('mdld'); Ivtk_type = 14
      end select
c
      end function
c
      function Ixdmf_type(Type)
      character(len=4) :: Type
c
      select case(Type)
      case ('vert'); Ixdmf_type = 0
      case ('medg'); Ixdmf_type = 0
      case ('mdlt'); Ixdmf_type = 4
      case ('mdlq'); Ixdmf_type = 5
      case ('mdln'); Ixdmf_type = 6
      case ('mdlb'); Ixdmf_type = 9
      case ('mdlp'); Ixdmf_type = 8
      case ('mdld'); Ixdmf_type = 7
      end select
c
      end function
c
      function Type_nod(Type, J)
      integer :: J
      character(len=4) :: Type_nod, Type;
      select case(Type)
      case ('mdln')
        select case(J)
        case (1,2,3,4);       Type_nod = 'vert'
        case (5,6,7, 8,9,10); Type_nod = 'medg'
        case (11,12,13,14);   Type_nod = 'mdlt'
        case (15);            Type_nod = 'mdln'
        end select
      case ('mdlp')
        select case(J)
        case (1,2,3, 4,5,6);               Type_nod = 'vert'
        case (7,8,9, 10,11,12, 13,14,15);  Type_nod = 'medg'
        case (16,17);                      Type_nod = 'mdlt'
        case (18,19,20);                   Type_nod = 'mdlq'
        case (21);                         Type_nod = 'mdlp'
        end select
      case ('mdlb')
        select case(J)
        case (1,2,3,4, 5,6,7,8);           Type_nod = 'vert'
        case (9,10,11,12, 13,14,15,16, 17,18,19,20)
          Type_nod = 'medg'
        case (21,22,23,24,25,26);          Type_nod = 'mdlq'
        case (27);                         Type_nod = 'mdlb'
        end select
      case ('mdld')
        select case(J)
        case (1,2,3,4,5);            Type_nod = 'vert'
        case (6,7,8,9, 10,11,12,13); Type_nod = 'medg'
        case (14);                   Type_nod = 'mdlq'
        case (15,16,17,18);          Type_nod = 'mdlt'
        case (19);                   Type_nod = 'mdld'
        end select
      case default
        write(*,*) 'Type_nod: wrong input type [', Type,']'
        stop 1
      end select
c
      end function Type_nod
c
c  ...determine number of vertices for a middle node or element type
      function Nvert(Type)
c
      character(len=4) :: Type
      integer Nvert
c
      select case(Type)
      case('medg'); Nvert=2
      case('mdlt','tria'); Nvert=3
      case('mdlq','quad'); Nvert=4
      case('mdlp','pris'); Nvert=6
      case('mdlb','bric'); Nvert=8
      case('mdln','tetr'); Nvert=4
      case('mdld','pyra'); Nvert=5
      end select
c
      end function Nvert
c
c
c  ...determine number of edges for a middle node or element  type
      function Nedge(Type)
c
      character(len=4) :: Type
      integer Nedge
c
      select case(Type)
      case('medg'); Nedge=0
      case('mdlt','tria'); Nedge=3
      case('mdlq','quad'); Nedge=4
      case('mdlp','pris'); Nedge=9
      case('mdlb','bric'); Nedge=12
      case('mdln','tetr'); Nedge=6
      case('mdld','pyra'); Nedge=8
      end select
c
      end function Nedge
c
c
c  ...determine number of faces for a middle node or element  type
      function Nface(Type)
c
      character(len=4) :: Type
      integer Nface
c
      select case(Type)
      case('mdlp','pris'); Nface=5
      case('mdlb','bric'); Nface=6
      case('mdln','tetr'); Nface=4
      case('mdld','pyra'); Nface=5
      end select
c
      end function Nface
c
c  ...return vertex numbers endpoints of an edge......
      subroutine edge_to_vert(Type,Ie, Nv1,Nv2)
      implicit none
      character(len=4), intent(in) :: Type
      integer, intent(in)  :: Ie
      integer, intent(out) :: Nv1,Nv2
c
      select case(Type)
      case('mdlt','tria')
        Nv1 = TRIAN_EDGE_TO_VERT(1,Ie); Nv2 = TRIAN_EDGE_TO_VERT(2,Ie)
      case('mdlq','rect')
        Nv1 = QUADR_EDGE_TO_VERT(1,Ie); Nv2 = QUADR_EDGE_TO_VERT(2,Ie)
      case('mdlp','pris')
        Nv1 = PRISM_EDGE_TO_VERT(1,Ie); Nv2 = PRISM_EDGE_TO_VERT(2,Ie)
      case('mdlb','bric')
        Nv1 = BRICK_EDGE_TO_VERT(1,Ie); Nv2 = BRICK_EDGE_TO_VERT(2,Ie)
      case('mdln','tetr')
        Nv1 = TETRA_EDGE_TO_VERT(1,Ie); Nv2 = TETRA_EDGE_TO_VERT(2,Ie)
      case('mdld','pyra')
        Nv1 = PYRAM_EDGE_TO_VERT(1,Ie); Nv2 = PYRAM_EDGE_TO_VERT(2,Ie)
      end select
      end subroutine edge_to_vert
c
c  ...return vertex numbers for a face
      subroutine face_to_vert(Type,Ifc, Nv1,Nv2,Nv3,Nv4)
      implicit none
      character(len=4), intent(in) :: Type
      integer, intent(in)  :: Ifc
      integer, intent(out) :: Nv1,Nv2,Nv3,Nv4
c
      select case(Type)
      case('mdlp','pris')
        Nv1 = PRISM_FACE_TO_VERT(1,Ifc);Nv2 = PRISM_FACE_TO_VERT(2,Ifc);
        Nv3 = PRISM_FACE_TO_VERT(3,Ifc);Nv4 = PRISM_FACE_TO_VERT(4,Ifc);
      case('mdlb','bric')
        Nv1 = BRICK_FACE_TO_VERT(1,Ifc);Nv2 = BRICK_FACE_TO_VERT(2,Ifc);
        Nv3 = BRICK_FACE_TO_VERT(3,Ifc);Nv4 = BRICK_FACE_TO_VERT(4,Ifc);
      case('mdln','tetr')
        Nv1 = TETRA_FACE_TO_VERT(1,Ifc);Nv2 = TETRA_FACE_TO_VERT(2,Ifc);
        Nv3 = TETRA_FACE_TO_VERT(3,Ifc);Nv4 = TETRA_FACE_TO_VERT(1,Ifc);
      case('mdld','pyra')
        Nv1 = PYRAM_FACE_TO_VERT(1,Ifc);Nv2 = PYRAM_FACE_TO_VERT(2,Ifc);
        Nv3 = PYRAM_FACE_TO_VERT(3,Ifc);Nv4 = PYRAM_FACE_TO_VERT(4,Ifc);
      end select
      end subroutine face_to_vert
c
c  ...return edge numbers for a face
      subroutine face_to_edge(Type,Ifc, Ne1,Ne2,Ne3,Ne4)
      implicit none
      character(len=4), intent(in) :: Type
      integer, intent(in)  :: Ifc
      integer, intent(out) :: Ne1,Ne2,Ne3,Ne4
c
      select case(Type)
      case('mdlp','pris')
        Ne1 = PRISM_FACE_TO_EDGE(1,Ifc);Ne2 = PRISM_FACE_TO_EDGE(2,Ifc);
        Ne3 = PRISM_FACE_TO_EDGE(3,Ifc);Ne4 = PRISM_FACE_TO_EDGE(4,Ifc);
      case('mdlb','bric')
        Ne1 = BRICK_FACE_TO_EDGE(1,Ifc);Ne2 = BRICK_FACE_TO_EDGE(2,Ifc);
        Ne3 = BRICK_FACE_TO_EDGE(3,Ifc);Ne4 = BRICK_FACE_TO_EDGE(4,Ifc);
      case('mdln','tetr')
        Ne1 = TETRA_FACE_TO_EDGE(1,Ifc);Ne2 = TETRA_FACE_TO_EDGE(2,Ifc);
        Ne3 = TETRA_FACE_TO_EDGE(3,Ifc);Ne4 = TETRA_FACE_TO_EDGE(1,Ifc);
      case('mdld','pyra')
        Ne1 = PYRAM_FACE_TO_EDGE(1,Ifc);Ne2 = PYRAM_FACE_TO_EDGE(2,Ifc);
        Ne3 = PYRAM_FACE_TO_EDGE(3,Ifc);Ne4 = PYRAM_FACE_TO_EDGE(4,Ifc);
      end select
      end subroutine face_to_edge
c
c----------------------------------------------------------------------
c
c> @param[in ] Ftype    - 2D element type = 3D element face type
c> @param[in ] Nface_or - the face orientation as seen by the element
c> @param[out] Nfver    - face edge numbers from the point of view
c                         of the face
c> @revision Feb 18
c----------------------------------------------------------------------
      subroutine face_to_vert_nos(Ftype,Nface_or, Nfver)
      implicit none
      character(len=4)     ,intent(in ) :: Ftype
      integer              ,intent(in ) :: Nface_or
      integer,dimension(4) ,intent(out) :: Nfver
c
      integer :: iv
c
      select case(Ftype)
      case('mdlt')
        do iv=1,3
          Nfver(iv) = TRIAN_L2G(iv,Nface_or)
        enddo
      case('mdlq')
        do iv=1,4
          Nfver(iv) = QUADR_L2G(iv,Nface_or)
        enddo
      end select
c
      end subroutine face_to_vert_nos

      subroutine face_to_edge_nos(Ftype,Nface_or, Nfedg)
      implicit none
      character(len=4)     ,intent(in ) :: Ftype
      integer              ,intent(in ) :: Nface_or
      integer,dimension(4) ,intent(out) :: Nfedg
c
      integer :: ie
c
      select case(Ftype)
      case('mdlt')
        do ie=1,3
          Nfedg(ie) = TRIAN_L2G_EDGE(ie,Nface_or)
        enddo
      case('mdlq')
        do ie=1,4
          Nfedg(ie) = QUADR_L2G_EDGE(ie,Nface_or)
        enddo
      end select
c
      end subroutine face_to_edge_nos




c
c----------------------------------------------------------------------
c
c> @param[in ] Ftype    - 2D element type = 3D element face type
c> @param[in ] Nface_or - the face orientation as seen by the element
c> @param[in ] Nedge_or - orientations for the face edges from
c                         the element point of view
c> @param[out] Nfedg_or - orientations for the face edges from
c                         the face point of view
c> @revision Feb 18
c----------------------------------------------------------------------
      subroutine face_to_edge_orient(Ftype,Nface_or,Nedge_or, Nfedg_or)
      implicit none
      character(len=4)     ,intent(in ) :: Ftype
      integer              ,intent(in ) :: Nface_or
      integer,dimension(4) ,intent(in ) :: Nedge_or
      integer,dimension(4) ,intent(out) :: Nfedg_or
c
      integer :: ie
c
      select case(Ftype)
      case('mdlt')
        do ie=1,3
          Nfedg_or(ie) = 
     .            mod(Nedge_or(ie)+TRIAN_ORIENT_L2G(ie,Nface_or),2)
        enddo
      case('mdlq')
        do ie=1,4
          Nfedg_or(ie) = 
     .            mod(Nedge_or(ie)+QUADR_ORIENT_L2G(ie,Nface_or),2)
        enddo
      end select
c
      end subroutine face_to_edge_orient
c
c----------------------------------------------------------------------
c
c  ...return local edge parametrizations......
      subroutine edge_param(Type,Ie,T, Xi,Dxidt)
c
      character(len=4) :: Type
      double precision T,Xi(*),Dxidt(*)
      double precision xi1(3),xi2(3)
c
      select case(Type)
      case('mdlt','tria')
        n1 = TRIAN_EDGE_TO_VERT(1,Ie);  n2 = TRIAN_EDGE_TO_VERT(2,Ie)
        xi1(1:2) = TRIAN_COORD(1:2,n1); xi2(1:2) = TRIAN_COORD(1:2,n2)
        nvar=2
      case('mdlq','quad')
        n1 = QUADR_EDGE_TO_VERT(1,Ie);  n2 = QUADR_EDGE_TO_VERT(2,Ie)
        xi1(1:2) = QUADR_COORD(1:2,n1); xi2(1:2) = QUADR_COORD(1:2,n2)
        nvar=2
      case('mdlp','pris')
        n1 = PRISM_EDGE_TO_VERT(1,Ie);  n2 = PRISM_EDGE_TO_VERT(2,Ie)
        xi1(1:3) = PRISM_COORD(1:3,n1); xi2(1:3) = PRISM_COORD(1:3,n2)
        nvar=3
      case('mdlb','bric')
        n1 = BRICK_EDGE_TO_VERT(1,Ie);  n2 = BRICK_EDGE_TO_VERT(2,Ie)
        xi1(1:3) = BRICK_COORD(1:3,n1); xi2(1:3) = BRICK_COORD(1:3,n2)
        nvar=3
      case('mdln','tetr')
        n1 = TETRA_EDGE_TO_VERT(1,Ie);  n2 = TETRA_EDGE_TO_VERT(2,Ie)
        xi1(1:3) = TETRA_COORD(1:3,n1); xi2(1:3) = TETRA_COORD(1:3,n2)
        nvar=3
      case('mdld','pyra')
        n1 = PYRAM_EDGE_TO_VERT(1,Ie);  n2 = PYRAM_EDGE_TO_VERT(2,Ie)
        xi1(1:3) = PYRAM_COORD(1:3,n1); xi2(1:3) = PYRAM_COORD(1:3,n2)
        nvar=3
      end select
c
      Dxidt(1:nvar) = xi2(1:nvar) - xi1(1:nvar)
      Xi(1:nvar) = xi1(1:nvar) + T*Dxidt(1:nvar)
      end subroutine edge_param
c
c  ...return local face parametrizations......
      subroutine face_param(Type,If,T, Xi,Dxidt)
c
      use control , only : GEOM_TOL
c
      character(len=4) :: Type
      double precision T(2),Xi(3),Dxidt(3,2)
      double precision xi1(3),xi2(3),xi3(3)
c
c  ...check that point is inside master face
      if (face_type(Type,If).eq.'tria') then
        if ((T(1)     .lt.-GEOM_TOL    ).or.
     .      (T(2)     .lt.-GEOM_TOL    ).or.
     .      (T(1)+T(2).gt.1.d0+GEOM_TOL)    ) then
          write(*,*)'face_param: point outside of tria face'
          write(*,1000)T(1:2)
 1000     format(' T = ',2(e12.5,x))
          call pause
        endif
          k=3
      endif
c
      if (face_type(Type,If).eq.'rect') then
        if ((T(1).lt.-GEOM_TOL).or.(T(1).gt.1.d0+GEOM_TOL).or.
     .      (T(2).lt.-GEOM_TOL).or.(T(2).gt.1.d0+GEOM_TOL)) then
          write(*,*)'face_param: point outside of rect face'
          write(*,1000)T(1:2)
          call pause
        endif
          k=4
      endif
c
      select case(Type)
      case('mdlp','pris')
        n = PRISM_FACE_TO_VERT(1,If);  xi1(1:3) = PRISM_COORD(1:3,n)
        n = PRISM_FACE_TO_VERT(2,If);  xi2(1:3) = PRISM_COORD(1:3,n)
        n = PRISM_FACE_TO_VERT(k,If);  xi3(1:3) = PRISM_COORD(1:3,n)
      case('mdlb','bric')
        n = BRICK_FACE_TO_VERT(1,If);  xi1(1:3) = BRICK_COORD(1:3,n)
        n = BRICK_FACE_TO_VERT(2,If);  xi2(1:3) = BRICK_COORD(1:3,n)
        n = BRICK_FACE_TO_VERT(k,If);  xi3(1:3) = BRICK_COORD(1:3,n)
      case('mdln','tetr')
        n = TETRA_FACE_TO_VERT(1,If);  xi1(1:3) = TETRA_COORD(1:3,n)
        n = TETRA_FACE_TO_VERT(2,If);  xi2(1:3) = TETRA_COORD(1:3,n)
        n = TETRA_FACE_TO_VERT(k,If);  xi3(1:3) = TETRA_COORD(1:3,n)
      case('mdld','pyra')
        n = PYRAM_FACE_TO_VERT(1,If);  xi1(1:3) = PYRAM_COORD(1:3,n)
        n = PYRAM_FACE_TO_VERT(2,If);  xi2(1:3) = PYRAM_COORD(1:3,n)
        n = PYRAM_FACE_TO_VERT(k,If);  xi3(1:3) = PYRAM_COORD(1:3,n)
      end select
c
      Dxidt(1:3,1) = xi2(1:3)-xi1(1:3)
      Dxidt(1:3,2) = xi3(1:3)-xi1(1:3)
      Xi(1:3) = xi1(1:3) + T(1)*Dxidt(1:3,1) + T(2)*Dxidt(1:3,2)
c
      end subroutine face_param
c
c  ...return normal signs corresponding to face parametrizations......
      integer function Nsign_param(Type,If)
c
      character(len=4) :: Type
      integer :: If
c
      select case(Type)
      case('mdlp','pris')
        select case(If)
        case(1,5)
          Nsign_param=-1
        case(2,3,4)
          Nsign_param= 1
        end select
      case('mdlb','bric')
        select case(If)
        case(1,5,6)
          Nsign_param=-1
        case(2,3,4)
          Nsign_param= 1
        end select
      case('mdln','tetr')
        select case(If)
        case(1,4)
          Nsign_param=-1
        case(2,3)
          Nsign_param= 1
        end select
      case('mdld','pyra')
        select case(If)
        case(1,4,5)
          Nsign_param=-1
        case(2,3)
          Nsign_param= 1
        end select
      end select
c
      end function Nsign_param
c
c  ...check if element face is triangle or rectangle
      function Face_type(Type,If)
c
      character(len=4) :: Face_type,Type
c
      select case(Type)
      case('mdlp','pris')
        select case(If)
          case(1,2);   Face_type = 'tria'
          case(3,4,5); Face_type = 'rect'
        end select
      case('mdlb','bric')
        Face_type = 'rect'
      case('mdln','tetr')
        Face_type = 'tria'
      case('mdld','pyra')
        select case(If)
          case(1); Face_type = 'rect'
          case(2,3,4,5); Face_type = 'tria'
        end select
      end select
      end function Face_type
c
c  ...return local face node numbers
      subroutine face_nodes(Type,If, Nface_nodes,Nrfn)
c
      character(len=4) :: Type
      integer If, Nface_nodes(9),Nrfn
c
      select case(Type)
      case('mdlp','pris')
        select case(If)
        case(1,2)
          Nface_nodes(1:3) = PRISM_FACE_TO_VERT(1:3,If)
          Nface_nodes(4:6) = 6+PRISM_FACE_TO_EDGE(1:3,If)
          Nface_nodes(7) = 15+If
          Nrfn=7
        case(3,4,5)
          Nface_nodes(1:4) = PRISM_FACE_TO_VERT(1:4,If)
          Nface_nodes(5:8) = 6+PRISM_FACE_TO_EDGE(1:4,If)
          Nface_nodes(9) = 15+If
          Nrfn=9
        end select
      case('mdlb','bric')
        Nface_nodes(1:4) = BRICK_FACE_TO_VERT(1:4,If)
        Nface_nodes(5:8) = 8+BRICK_FACE_TO_EDGE(1:4,If)
        Nface_nodes(9) = 20+If
        Nrfn=9
      case('mdln','tetr')
        Nface_nodes(1:3) = TETRA_FACE_TO_VERT(1:3,If)
        Nface_nodes(4:6) = 4+TETRA_FACE_TO_EDGE(1:3,If)
        Nface_nodes(7) = 10+If
        Nrfn=7
      case('mdld','pyra')
        select case(If)
        case(1)
          Nface_nodes(1:4) = PYRAM_FACE_TO_VERT(1:4,If)
          Nface_nodes(5:8) = 5+PYRAM_FACE_TO_EDGE(1:4,If)
          Nface_nodes(9) = 13+If
          Nrfn=9
        case(2,3,4,5)
          Nface_nodes(1:3) = PYRAM_FACE_TO_VERT(1:3,If)
          Nface_nodes(4:6) = 5+PYRAM_FACE_TO_EDGE(1:3,If)
          Nface_nodes(7) = 13+If
          Nrfn=7
        end select
      end select
      end subroutine face_nodes
c
c  ...return local face order of approximation
      subroutine face_order(Type,If,Norder, Nordf)
c
      character(len=4) :: Type
      integer If, Nordf(5),Norder(19)
c
      select case(Type)
      case('mdlp','pris')
        select case(If)
        case(1,2)
          do i=1,3
            k = PRISM_FACE_TO_EDGE(i,If)
            Nordf(i) = Norder(k)
          enddo
          Nordf(4) = Norder(9+If)
        case(3,4,5)
          do i=1,4
            k = PRISM_FACE_TO_EDGE(i,If)
            Nordf(i) = Norder(k)
          enddo
          Nordf(5) = Norder(9+If)
        end select
      case('mdlb','bric')
        do i=1,4
          k = BRICK_FACE_TO_EDGE(i,If)
          Nordf(i) = Norder(k)
        enddo
        Nordf(5) = Norder(12+If)
      case('mdln','tetr')
        do i=1,3
          k = TETRA_FACE_TO_EDGE(i,If)
          Nordf(i) = Norder(k)
        enddo
        Nordf(4) = Norder(6+If)
      case('mdld','pyra')
        select case(If)
        case(1)
          do i=1,4
            k = PYRAM_FACE_TO_EDGE(i,If)
            Nordf(i) = Norder(k)
          enddo
          Nordf(5) = Norder(8+If)
        case(2,3,4,5)
          do i=1,3
            k = PYRAM_FACE_TO_EDGE(i,If)
            Nordf(i) = Norder(k)
          enddo
          Nordf(4) = Norder(8+If)
        end select
      end select
      end subroutine face_order
c
c  ...return number of dof for a single component
      subroutine ndof_nod(Type,Nord, NdofH,NdofE,NdofV,NdofQ)
c
      character(len=4) :: Type
      integer  Nord,NdofH,NdofE,NdofV,NdofQ
      integer nordx,nordy,nordz,naux
c
      select case(Type)
c  ...VERTEX
      case('vert')
        NdofH=1; NdofE=0; NdofV=0; NdofQ=0
c  ...EDGE
      case('medg')
        NdofH=Nord-1; NdofE=Nord; NdofV=0; NdofQ=0
c  ...TRIANGLE
      case('tria','mdlt')
        NdofH=(Nord-2)*(Nord-1)/2
        NdofE=(Nord-1)*Nord
        NdofV=Nord*(Nord+1)/2
        NdofQ=0
c  ...QUAD
      case('rect','mdlq')
        call decode(Nord, nordx,nordy)
        NdofH=(nordx-1)*(nordy-1)
        NdofE=nordx*(nordy-1)+nordy*(nordx-1)
        NdofV=nordx*nordy
        NdofQ=0
c  ...PRISM
      case('pris','mdlp')
        call decode(Nord, nordx,nordz)
        NdofH=(nordx-2)*(nordx-1)/2*(nordz-1)
        NdofE= (nordx-1)*nordx*(nordz-1) + (nordx-2)*(nordx-1)/2*nordz
        NdofV=(nordx-1)*nordx*nordz + nordx*(nordx+1)/2*(nordz-1)
        NdofQ=(nordx+1)*nordx/2*nordz
c  ...HEXA
      case('bric','mdlb')
        call decode(Nord, naux,nordz)
        call decode(naux, nordx,nordy)
        NdofH=(nordx-1)*(nordy-1)*(nordz-1)
        NdofE= nordx*(nordy-1)*(nordz-1) + (nordx-1)*nordy*(nordz-1)
     .           + (nordx-1)*(nordy-1)*nordz
        NdofV=(nordx-1)*nordy*nordz + nordx*(nordy-1)*nordz
     .           + nordx*nordy*(nordz-1)
        NdofQ=nordx*nordy*nordz
c  ...TETRA
      case('tetr','mdln')
        NdofH=(Nord-3)*(Nord-2)*(Nord-1)/6
        NdofE=(Nord-2)*(Nord-1)*Nord/2
        NdofV=(Nord-1)*Nord*(Nord+1)/2
        NdofQ=Nord*(Nord+1)*(Nord+2)/6
c  ...PYRAMID
      case('pyra','mdld')
        n = Nord-1
c  LD: value for the original pyramid of Paolo
ccc        NdofH=(2*n**3+3*n**2+n)/6
        NdofH=(Nord-1)**3
        NdofE=3*(Nord-1)**2*Nord
c        NdofV=2*(Nord-1)*Nord**2+ Nord**3
        NdofV=3*(Nord-1)*Nord**2
c        NdofQ=Nord**2*(Nord+1)
        NdofQ=Nord**3

      case default
        write(*,*) 'Error!!! ndof_nod: Type = ',Type
        NdofH=0
        NdofE= 0
        NdofV=0
        NdofQ=0
        call pause
        return
      end select
ccc      write(*,*) 'ndof_nod: Type,Nord, NdofH,NdofE,NdofV,NdofQ = ',
ccc     .                      Type,Nord, NdofH,NdofE,NdofV,NdofQ
ccc      call pause
c
      end subroutine ndof_nod
c
c
      end module element_data
c
c
c
c======================================================================
c  REMARK : the following subroutines are somehow related to module   |
c           "element_data". They are (meant to be) robust wrt the     |
c           info contained in "element_data".                         |
c======================================================================
c
c
c
c----------------------------------------------------------------------
c> Purpose : given a point in a 2D master element having coordinates
c            T(1:2) wrt a LOCAL system of coordinates, routine
c            returns coordinates Eta(1:2) wrt a GLOBAL system having
c            orientation Norient wrt the LOCAL system. Orientations
c            are defined by TRIAN,QUADR_L2G.
c
c> @param[in ] Type    - 2D element type
c> @param[in ] T       - coordinates in the LOCAL system
c> @param[in ] Norient - orientation of GLOBAL system wrt LOCAL system
c> @param[out] Eta     - coordinates in the GLOBAL system
c> @param[out] Detadt  - derivatives
c
c> @revision Dec 12
c----------------------------------------------------------------------
      subroutine local2global(Type,T,Norient, Eta,Detadt)
c
      use element_data , only : TRIAN_L2G,QUADR_L2G,QUADR_EDGE_TO_VERT
c
      implicit none
      character(len=4)     ,intent(in ) :: Type
      real*8,dimension(2  ),intent(in ) :: T
      integer              ,intent(in ) :: Norient
      real*8,dimension(2  ),intent(out) :: Eta
      real*8,dimension(2,2),intent(out) :: Detadt
c
      real*8,dimension(  4) :: rlam
      real*8,dimension(2,4) :: drlam
      integer :: iv,iv1,iv2,iv3,iv4,iedge,i
c
c----------------------------------------------------------------------
c
c  ...initialize
      Eta(1:2)=0.d0 ; Detadt(1:2,1:2)=0.d0
c
c  ...vertex shape functions
      call vshape2(Type,T, rlam,drlam)
c
      select case(Type)
c
c  ...TRIANGLE.........................................................
      case('tria')
        iv2=TRIAN_L2G(2,Norient) ; iv3=TRIAN_L2G(3,Norient)
        Eta(1)=rlam(iv2) ; Detadt(1,1:2)=drlam(1:2,iv2)
        Eta(2)=rlam(iv3) ; Detadt(2,1:2)=drlam(1:2,iv3)
c
c  ...QUAD.............................................................
      case('quad')
        iv1=QUADR_L2G(1,Norient)
        iv2=QUADR_L2G(2,Norient)
        iv4=QUADR_L2G(4,Norient)
c
c  .....projection on (iv1,iv2) edge
        call quad_aux(iv1,iv2, iedge)
        do i=1,2
          iv=QUADR_EDGE_TO_VERT(i,iedge)
          Eta   (1)    =Eta   (1)    +rlam     (iv)
          Detadt(1,1:2)=Detadt(1,1:2)+drlam(1:2,iv)
        enddo
c  .....projection on (iv1,iv4) edge
        call quad_aux(iv1,iv4, iedge)
        do i=1,2
          iv=QUADR_EDGE_TO_VERT(i,iedge)
          Eta(   2    )=Eta(   2    )+rlam(     iv)
          Detadt(2,1:2)=Detadt(2,1:2)+drlam(1:2,iv)
        enddo
c
      case default
        write(*,7001) Type
 7001   format(' local2global: invalid Type = ',a10)
        stop
      endselect
c
c
      endsubroutine local2global
c
c
c
c----------------------------------------------------------------------
c> Purpose : auxiliary routine needed by "hp3gen"
c
c> @param[in ] Type    - element type
c> @param[in ] Ivert   - vertex number
c> @param[out] Nrfaces - number of adjacent faces
c> @param[out] Nofaces - adjacent face numbers
c
c> @revision Nov 12
c----------------------------------------------------------------------
c
      subroutine vert_to_faces(Type,Ivert, Nrfaces,Nofaces)
c
      use element_data , only : PRISM_VERT_TO_FACE,
     .                          BRICK_VERT_TO_FACE,
     .                          TETRA_VERT_TO_FACE,
     .                          PYRAM_VERT_TO_FACE
c
      implicit none
      character(len=4)    ,intent(in ) :: Type
      integer             ,intent(in ) :: Ivert
      integer             ,intent(out) :: Nrfaces
      integer,dimension(4),intent(out) :: Nofaces
c----------------------------------------------------------------------
c
c  ...number of connected faces
      Nrfaces=3
      select case(Type)
      case('bric','mdlb')
        Nofaces(1:3)=BRICK_VERT_TO_FACE(1:3,Ivert)
      case('tetr','mdln')
        Nofaces(1:3)=TETRA_VERT_TO_FACE(1:3,Ivert)
      case('pris','mdlp')
        Nofaces(1:3)=PRISM_VERT_TO_FACE(1:3,Ivert)
      case('pyra','mdld')
        select case(Ivert)
c  .....lateral triangular faces
        case(1,2,3,4)
          Nofaces(1:3)=PYRAM_VERT_TO_FACE(1:3,Ivert)
c  .....rectangular base
        case(5)
          Nrfaces=4
          Nofaces(1:4)=PYRAM_VERT_TO_FACE(1:4,Ivert)
        endselect
      endselect
c
      end subroutine vert_to_faces
c
c
c----------------------------------------------------------------------
c> Purpose : auxiliary routine needed by "hp3gen"
c
c> @param[in ] Type    - element type
c> @param[in ] Iedge   - edge number
c> @param[out] Nofaces - adjacent face numbers
c
c> @revision Nov 12
c----------------------------------------------------------------------
c
      subroutine edge_to_faces(Type,Iedge, Nofaces)
c
      use element_data , only : PRISM_EDGE_TO_FACE,
     .                          BRICK_EDGE_TO_FACE,
     .                          TETRA_EDGE_TO_FACE,
     .                          PYRAM_EDGE_TO_FACE
c
      implicit none
      character(len=4)    ,intent(in ) :: Type
      integer             ,intent(in ) :: Iedge
      integer,dimension(2),intent(out) :: Nofaces
c----------------------------------------------------------------------
c
      select case(Type)
      case('bric','mdlb')
        Nofaces(1:2)=BRICK_EDGE_TO_FACE(1:2,Iedge)
      case('tetr','mdln')
        Nofaces(1:2)=TETRA_EDGE_TO_FACE(1:2,Iedge)
      case('pris','mdlp')
        Nofaces(1:2)=PRISM_EDGE_TO_FACE(1:2,Iedge)
      case('pyra','mdld')
        Nofaces(1:2)=PYRAM_EDGE_TO_FACE(1:2,Iedge)
      end select
c
c
      end subroutine edge_to_faces
c
c
c----------------------------------------------------------------------
      subroutine quad_aux(Iv1,Iv2, Iedge)
c
      use element_data , only : QUADR_EDGE_TO_VERT
c
      implicit none
      integer,intent(in ) :: Iv1,Iv2
      integer,intent(out) :: Iedge
c
      integer,save,dimension(4,4) :: AUX
      logical,save                :: INITIALIZED=.FALSE.
      integer,       dimension(2) :: iedge1,iedge2,list
      integer                     :: i,j,i1,i2,ifound
c----------------------------------------------------------------------
c
      if (.NOT.INITIALIZED) then
c
c  ...initialize
      AUX(1:4,1:4)=0
c
c  ...double loop over edges
      do i=1,4
c
c  .....1st edge endpoints
        i1=QUADR_EDGE_TO_VERT(1,i)
        i2=QUADR_EDGE_TO_VERT(2,i)
        iedge1(1)=i1 ; iedge1(2)=i2
c
        do j=1,4
c
c  .......2nd edge endpoints
          iedge2(1)=QUADR_EDGE_TO_VERT(1,j)
          iedge2(2)=QUADR_EDGE_TO_VERT(2,j)
c
c  .......compare and record
          call list_minus(iedge1(1:2),iedge2(1:2),2, list,ifound)
          if (ifound.ne.1) cycle
c
          if (list(1).eq.i1)  AUX(i1,i2)=j
          if (list(1).eq.i2)  AUX(i2,i1)=j
        enddo
      enddo
c
ccc  ...printing
cc      do i=1,4 ; do j=1,4
ccc  .....print only relevant entries
cc        if (AUX(i,j).eq.0) cycle
cc        write(*,1000) i,j,AUX(i,j)
cc 1000   format(' iv1,2, edge = '2(i1,2x),2x,i1)
cc      enddo    ; enddo
c
c  ...update initialization flag
      INITIALIZED=.TRUE.
c
      endif
c
c  ...retrive edge number
      Iedge=AUX(Iv1,Iv2)
c
c
      endsubroutine quad_aux
c
c
c----------------------------------------------------------------------
c> Pupose : {List1} - {List2} = {List3}
c
c> @param[in ] List1,2 - lists of integers
c> @param[in ] N       - length of List1,2
c> @param[out] List3   - difference
c> @param[out] J       - length of List3
c
c> @revision Nov 12
c----------------------------------------------------------------------
      subroutine list_minus(List1,List2,N, List3,J)
c
      implicit none
      integer,dimension(N),intent(in ) :: List1,List2
      integer             ,intent(in ) :: N
      integer,dimension(N),intent(out) :: List3
      integer             ,intent(out) :: J
      integer :: i,ifound
c----------------------------------------------------------------------
c
c  ...initialize
      J=0 ; List3(1:N)=0
c
c  ...loop over entries of List1
      do i=1,N
c  .....locate entry of List1 on List2
        call locate(List1(i),List2(1:N),N, ifound)
c  .....if entry was not found, add it to List3
        if (ifound.eq.0) then
          J=J+1 ; List3(J)=List1(i)
        endif
      enddo
c
c
      endsubroutine list_minus
c
c
c----------------------------------------------------------------------
      subroutine hexa_aux(Iv1,Iv2, Iface)
c
      use element_data , only : BRICK_EDGE_TO_VERT,
     .                          BRICK_VERT_TO_FACE
c
      implicit none
      integer,intent(in ) :: Iv1,Iv2
      integer,intent(out) :: Iface
c
      integer,save,dimension(8,8) :: AUX
      logical,save                :: INITIALIZED=.FALSE.
      integer,       dimension(3) :: iface1,iface2,list
      integer                     :: i,i1,i2,ifound
c----------------------------------------------------------------------
c
      if (.NOT.INITIALIZED) then
c
c  ...initialize
      AUX(1:8,1:8)=0
c
c  ...loop over edges
      do i=1,12
        i1=BRICK_EDGE_TO_VERT(1,i)
        i2=BRICK_EDGE_TO_VERT(2,i)
c
c  .....faces connected to 1st and 2nd vertex
        iface1(1:3)=BRICK_VERT_TO_FACE(1:3,i1)
        iface2(1:3)=BRICK_VERT_TO_FACE(1:3,i2)
c
c  .....determine non-common face, and store
        call list_minus(iface2(1:3),iface1(1:3),3, list,ifound)
        if (ifound.eq.0) then
          write(*,*)'hexa_aux: no face found!'
          stop
        endif
        AUX(i1,i2)=list(1)
c
        call list_minus(iface1(1:3),iface2(1:3),3, list,ifound)
        if (ifound.eq.0) then
          write(*,*)'hexa_aux: no face found!'
          stop
        endif
        AUX(i2,i1)=list(1)
      enddo
c
cccc  ...printing
ccc      do i=1,8 ; do j=1,8
cccc  .....print only relevant entries
ccc        if (AUX(i,j).eq.0) cycle
ccc        write(*,1000) i,j,AUX(i,j)
ccc 1000   format(' iv1,2, face = '2(i1,2x),2x,i1)
ccc      enddo    ; enddo
c
c  ...update initialization flag
      INITIALIZED=.TRUE.
c
      endif
c
c  ...retrive face number
      Iface=AUX(Iv1,Iv2)
c
c
      end subroutine hexa_aux
