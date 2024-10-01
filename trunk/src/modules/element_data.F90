!----------------------------------------------------------------------
!> @brief   element data for elements of different type
!> @date    Feb 2024
!----------------------------------------------------------------------
!
   module element_data
!
      use node_types
!
      implicit none
!
!----------------------------------------------------------------------
!  MASTER ELEMENTS VERTEX COORDINATES (vertex enumeration)            |
!----------------------------------------------------------------------
      real(8), parameter, dimension(2,3) :: TRIAN_COORD = &
      reshape( &
      (/0.d0,0.d0, 1.d0,0.d0, 0.d0,1.d0/) &
      ,(/2,3/))
!
      real(8), parameter, dimension(2,4) :: QUADR_COORD = &
      reshape( &
      (/0.d0,0.d0, 1.d0,0.d0, 1.d0,1.d0, 0.d0,1.d0/) &
      ,(/2,4/))
!
      real(8), parameter, dimension(3,6) :: PRISM_COORD = &
      reshape( &
      (/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, &
        0.d0,0.d0,1.d0, 1.d0,0.d0,1.d0, 0.d0,1.d0,1.d0/) &
      ,(/3,6/))
!
      real(8), parameter, dimension(3,8) :: BRICK_COORD = &
      reshape( &
      (/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 1.d0,1.d0,0.d0, 0.d0,1.d0,0.d0, &
        0.d0,0.d0,1.d0, 1.d0,0.d0,1.d0, 1.d0,1.d0,1.d0, 0.d0,1.d0,1.d0/) &
      ,(/3,8/))
!
      real(8), parameter, dimension(3,4) :: TETRA_COORD = &
      reshape( &
      (/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, &
        0.d0,0.d0,1.d0/) &
      ,(/3,4/))
!
      real(8), parameter, dimension(3,5) :: PYRAM_COORD = &
      reshape( &
      (/0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 1.d0,1.d0,0.d0, 0.d0,1.d0,0.d0, &
        0.d0,0.d0,1.d0/) &
      ,(/3,5/))
!
!----------------------------------------------------------------------
!  EDGE_2_VERT CONNECTIVITIES (edge enumeration & orientation)        |
!----------------------------------------------------------------------
      integer, parameter, dimension(2,3) :: TRIAN_EDGE_TO_VERT = &
      reshape( &
      (/1,2, 2,3, 1,3/) &
      ,(/2,3/))
!
      integer, parameter, dimension(2,4) :: QUADR_EDGE_TO_VERT = &
      reshape( &
      (/1,2, 2,3, 4,3, 1,4/) &
      ,(/2,4/))
!
      integer, parameter, dimension(2,9) :: PRISM_EDGE_TO_VERT = &
      reshape( &
      (/1,2, 2,3, 1,3, 4,5, 5,6, 4,6, 1,4, 2,5, 3,6/) &
      ,(/2,9/))
!
      integer, parameter, dimension(2,12) :: BRICK_EDGE_TO_VERT = &
      reshape( &
      (/1,2, 2,3, 4,3, 1,4, 5,6, 6,7, 8,7, 5,8, 1,5, 2,6, 3,7, 4,8/) &
      ,(/2,12/))
!
      integer, parameter, dimension(2,6) :: TETRA_EDGE_TO_VERT = &
      reshape( &
      (/1,2, 2,3, 1,3, 1,4, 2,4, 3,4/) &
      ,(/2,6/))
!
      integer, parameter, dimension(2,8) :: PYRAM_EDGE_TO_VERT = &
      reshape( &
      (/1,2, 2,3, 4,3, 1,4, 1,5, 2,5, 3,5, 4,5/) &
      ,(/2,8/))
!
!----------------------------------------------------------------------
!  FACE_2_VERT CONNECTIVITIES (face enumeration & orientation)        |
!----------------------------------------------------------------------
      integer, parameter, dimension(4,5) :: PRISM_FACE_TO_VERT = &
      reshape( &
      (/1,2,3,1, 4,5,6,4, 1,2,5,4, 2,3,6,5, 1,3,6,4/) &
      ,(/4,5/))
!
      integer, parameter, dimension(4,6) :: BRICK_FACE_TO_VERT = &
      reshape( &
      (/1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,4,8,5/) &
      ,(/4,6/))
!
      integer, parameter, dimension(3,4) :: TETRA_FACE_TO_VERT = &
      reshape( &
      (/1,2,3, 1,2,4, 2,3,4, 1,3,4/) &
      ,(/3,4/))
!
      integer, parameter, dimension(4,5) :: PYRAM_FACE_TO_VERT = &
      reshape( &
      (/1,2,3,4, 1,2,5,1, 2,3,5,2, 4,3,5,4, 1,4,5,1/) &
      ,(/4,5/))
!
!----------------------------------------------------------------------
!  FACE_2_EDGE CONNECTIVITIES (redundant!)                            |
!----------------------------------------------------------------------
      integer, parameter, dimension(4,5) :: PRISM_FACE_TO_EDGE = &
      reshape( &
      (/1,2,3,1, 4,5,6,4, 1,8,4,7, 2,9,5,8, 3,9,6,7/) &
      ,(/4,5/))
!
      integer, parameter, dimension(4,6) :: BRICK_FACE_TO_EDGE = &
      reshape( &
      (/1,2,3,4, 5,6,7,8, 1,10,5,9, 2,11,6,10, 3,11,7,12, 4,12,8,9/) &
      ,(/4,6/))
!
      integer, parameter, dimension(3,4) :: TETRA_FACE_TO_EDGE = &
      reshape( &
      (/1,2,3, 1,5,4, 2,6,5, 3,6,4/) &
      ,(/3,4/))
!
      integer, parameter, dimension(4,5) :: PYRAM_FACE_TO_EDGE = &
      reshape( &
      (/1,2,3,4, 1,6,5,1, 2,7,6,2, 3,7,8,3, 4,8,5,4/) &
      ,(/4,5/))
!
!----------------------------------------------------------------------
!  VERT_2_EDGE CONNECTIVITIES (redundant!)                            |
!----------------------------------------------------------------------
      integer, parameter, dimension(3,4) :: TETRA_VERT_TO_EDGE = &
      reshape( &
      (/1,3,4, 1,2,5, 2,3,6, 4,5,6/) &
      ,(/3,4/))
!
      integer, parameter, dimension(3,8) :: BRICK_VERT_TO_EDGE = &
      reshape( &
      (/4,1,9, 1,2,10, 2,3,11, 3,4,12, 8,5,9, 5,6,10, 6,7,11, 7,8,12/) &
      ,(/3,8/))
!
      integer, parameter, dimension(3,6) :: PRISM_VERT_TO_EDGE = &
      reshape( &
      (/1,3,7, 1,2,8, 2,3,9, 4,6,7, 4,5,8, 5,6,9/) &
      ,(/3,6/))
!
!----------------------------------------------------------------------
!  VERT_2_FACE CONNECTIVITIES (redundant!)                            |
!----------------------------------------------------------------------
      integer, parameter, dimension(3,6) :: PRISM_VERT_TO_FACE = &
      reshape( &
      (/1,3,5, 1,3,4, 1,4,5, 2,3,5, 2,3,4, 2,4,5/) &
      ,(/3,6/))
!
      integer, parameter, dimension(3,8) :: BRICK_VERT_TO_FACE = &
      reshape( &
      (/1,6,3, 1,3,4, 1,4,5, 1,5,6, 2,6,3, 2,3,4, 2,4,5, 2,5,6/) &
      ,(/3,8/))
!
      integer, parameter, dimension(3,4) :: TETRA_VERT_TO_FACE = &
      reshape( &
      (/1,2,4, 1,2,3, 1,3,4, 2,3,4/) &
      ,(/3,4/))
!
      integer, parameter, dimension(4,5) :: PYRAM_VERT_TO_FACE = &
      reshape( &
      (/1,2,5,0, 1,2,3,0, 1,3,4,0, 1,4,5,0, 2,3,4,5/) &
      ,(/4,5/))
!
!----------------------------------------------------------------------
!  EDGE_2_FACE CONNECTIVITIES (redundant!)                            |
!----------------------------------------------------------------------
      integer, parameter, dimension(2,9) :: PRISM_EDGE_TO_FACE = &
      reshape( &
      (/1,3, 1,4, 1,5, 2,3, 2,4, 2,5, 3,5, 3,4, 4,5/) &
      ,(/2,9/))
!
      integer, parameter, dimension(2,12) :: BRICK_EDGE_TO_FACE = &
      reshape( &
      (/1,3, 1,4, 1,5, 1,6, 2,3, 2,4, 2,5, 2,6, 3,6, 3,4, 4,5, 5,6/) &
      ,(/2,12/))
!
      integer, parameter, dimension(2,6) :: TETRA_EDGE_TO_FACE = &
      reshape( &
      (/1,2, 1,3, 1,4, 2,4, 2,3, 3,4/) &
      ,(/2,6/))
!
      integer, parameter, dimension(2,8) :: PYRAM_EDGE_TO_FACE = &
      reshape( &
      (/1,2, 1,3, 1,4, 1,5, 2,5, 2,3, 3,4, 4,5/) &
      ,(/2,8/))
!
!----------------------------------------------------------------------
!  LOCAL_2_GLOBAL ORIENTATIONS                                        |
!----------------------------------------------------------------------
      integer, parameter, dimension(2,0:1) :: EDGE_L2G = &
      reshape( &
      (/1,2, 2,1/) &
      ,(/2,2/))
!
      integer, parameter, dimension(3,0:5) :: TRIAN_L2G = &
      reshape( &
      (/1,2,3, 2,3,1, 3,1,2, 1,3,2, 2,1,3, 3,2,1/) &
      ,(/3,6/))
!
      integer, parameter, dimension(3,0:5) :: TRIAN_ORIENT_L2G = &
      reshape( &
      (/0,0,0, 0,1,1, 1,0,1, 0,1,0, 1,0,0, 1,1,1/) &
      ,(/3,6/))

      integer, parameter, dimension(3,0:5) :: TRIAN_L2G_EDGE = &
      reshape( &
      (/1,2,3, 2,3,1, 3,1,2, 3,2,1, 1,3,2, 2,1,3/) &
      ,(/3,6/))
!
      integer, parameter, dimension(4,0:7) :: QUADR_L2G_EDGE = &
      reshape( &
      (/1,2,3,4, 2,3,4,1, 3,4,1,2, 4,1,2,3, &
        4,3,2,1, 1,4,3,2, 2,1,4,3, 3,2,1,4/) &
      ,(/4,8/))


      integer, parameter, dimension(4,0:7) :: QUADR_L2G = &
      reshape( &
      (/1,2,3,4, 2,3,4,1, 3,4,1,2, 4,1,2,3, &
        1,4,3,2, 2,1,4,3, 3,2,1,4, 4,3,2,1/) &
      ,(/4,8/))
!
      integer, parameter, dimension(4,0:7) :: QUADR_ORIENT_L2G = &
      reshape( &
      (/0,0,0,0, 0,1,0,1, 1,1,1,1, 1,0,1,0, &
        0,0,0,0, 1,0,1,0, 1,1,1,1, 0,1,0,1/) &
      ,(/4,8/))
!
!----------------------------------------------------------------------
!  QUAD_2_AXES ORIENTATIONS                                           |
!                                                                     |
!  NFAXES(3,j)=0,[1] : x and y axes have not been [HAVE BEEN] swapped |
!                      for j-th orientation                           |
!                                                                     |
!  NFAXES(1,j)=0,[1] : consistent [INCONSISTENT] orientation for      |
!                      horizontal axis, for j-th orientation          |
!                                                                     |
!  NFAXES(2,j)=0,[1] : consistent [INCONSISTENT] orientation for      |
!                      vertical axis, for j-th orientation            |
!----------------------------------------------------------------------
      integer, parameter, dimension(3,0:7) :: NFAXES = &
      reshape( &
      (/0,0,0, 1,0,1, 1,1,0, 0,1,1, 0,0,1, 1,0,0, 1,1,1, 0,1,0/) &
      ,(/3,8/))
!
!----------------------------------------------------------------------
!  ...defines vertex shape functions; for vertex iv,
!     IJKV(ixi,iv) = number of 1D linear shape function in Xi(ixi)
!
!  ...defines edge shape functions; for edge ie
!     IXIEDGE(ie) = edge coordinate
!     IBLENDE(2,ie) = blending directions
!     NBLENDE(2,ie) = numbers of blending linear shape function
!
!  ...defines face shape functions; for face if
!     IXIFACE(2,if) = face coordinates
!     IBLENDF(if) = blending directions
!     NBLENDF(if) = numbers of blending linear shape function
!
      integer, parameter, dimension(3,8) :: IJKV = &
      reshape( &
      (/1,1,1, 2,1,1, 2,2,1, 1,2,1, &
        1,1,2, 2,1,2, 2,2,2, 1,2,2/) &
      ,(/3,8/))
!
      integer, parameter, dimension(12) :: IXIEDGE = &
      (/1,2,1,2, 1,2,1,2, 3,3,3,3/)
!
      integer, parameter, dimension(2,12) :: IBLENDE = &
      reshape( &
      (/2,3, 1,3, 2,3, 1,3, 2,3, 1,3, &
        2,3, 1,3, 1,2, 1,2, 1,2, 1,2/) &
      ,(/2,12/))
!
      integer, parameter, dimension(2,12) :: NBLENDE = &
      reshape( &
      (/1,1, 2,1, 2,1, 1,1, 1,2, 2,2, &
        2,2, 1,2, 1,1, 2,1, 2,2, 1,2/) &
      ,(/2,12/))
!
      integer, parameter, dimension(2,6) :: IXIFACE = &
      reshape( &
      (/1,2, 1,2, 1,3, 2,3, 1,3, 2,3/) &
      ,(/2,6/))
!
      integer, parameter, dimension(6) :: IBLENDF = &
      (/3,3,2,1,2,1/)
!
      integer, parameter, dimension(6) :: NBLENDF = &
       (/1,2,1,2,2,1/)
!
!----------------------------------------------------------------------
!
      contains
!
!----------------------------------------------------------------------
!
      function Ivtk_type(Ntype)
!
      integer :: Ivtk_type
      integer :: Ntype
!
      select case(Ntype)
        case (VERT); Ivtk_type = 1
        case (MEDG); Ivtk_type = 2
        case (MDLT); Ivtk_type = 5
        case (MDLQ); Ivtk_type = 9
        case (MDLN); Ivtk_type = 10
        case (MDLB); Ivtk_type = 12
        case (MDLP); Ivtk_type = 13
        case (MDLD); Ivtk_type = 14
        case default
          write(*,*) 'Ivtk_type'; stop
      end select
!
      end function Ivtk_type
!
!----------------------------------------------------------------------
!
      function Ixdmf_type(Ntype)
!
      integer :: Ixdmf_type
      integer :: Ntype
!
      select case(Ntype)
        case (VERT); Ixdmf_type = 0
        case (MEDG); Ixdmf_type = 0
        case (MDLT); Ixdmf_type = 4
        case (MDLQ); Ixdmf_type = 5
        case (MDLN); Ixdmf_type = 6
        case (MDLB); Ixdmf_type = 9
        case (MDLP); Ixdmf_type = 8
        case (MDLD); Ixdmf_type = 7
        case default
          write(*,*) 'Ixdmf_type'; stop
      end select
!
      end function Ixdmf_type
!
!----------------------------------------------------------------------
!
!  ...return vertex numbers endpoints of an edge......
   subroutine edge_to_vert(Ntype,Ie, Nv1,Nv2)
!
      integer, intent(in)  :: Ntype,Ie
      integer, intent(out) :: Nv1,Nv2
!
      select case(Ntype)
      case(MDLT,TRIA)
        Nv1 = TRIAN_EDGE_TO_VERT(1,Ie); Nv2 = TRIAN_EDGE_TO_VERT(2,Ie)
      case(MDLQ,RECT)
        Nv1 = QUADR_EDGE_TO_VERT(1,Ie); Nv2 = QUADR_EDGE_TO_VERT(2,Ie)
      case(MDLP,PRIS)
        Nv1 = PRISM_EDGE_TO_VERT(1,Ie); Nv2 = PRISM_EDGE_TO_VERT(2,Ie)
      case(MDLB,BRIC)
        Nv1 = BRICK_EDGE_TO_VERT(1,Ie); Nv2 = BRICK_EDGE_TO_VERT(2,Ie)
      case(MDLN,TETR)
        Nv1 = TETRA_EDGE_TO_VERT(1,Ie); Nv2 = TETRA_EDGE_TO_VERT(2,Ie)
      case(MDLD,PYRA)
        Nv1 = PYRAM_EDGE_TO_VERT(1,Ie); Nv2 = PYRAM_EDGE_TO_VERT(2,Ie)
      end select
!
   end subroutine edge_to_vert
!
!----------------------------------------------------------------------
!
!  ...return vertex numbers for a face
   subroutine face_to_vert(Ntype,Ifc, Nv1,Nv2,Nv3,Nv4)
!
      integer, intent(in)  :: Ntype,Ifc
      integer, intent(out) :: Nv1,Nv2,Nv3,Nv4
!
      select case(Ntype)
      case(MDLP,PRIS)
        Nv1 = PRISM_FACE_TO_VERT(1,Ifc);Nv2 = PRISM_FACE_TO_VERT(2,Ifc);
        Nv3 = PRISM_FACE_TO_VERT(3,Ifc);Nv4 = PRISM_FACE_TO_VERT(4,Ifc);
      case(MDLB,BRIC)
        Nv1 = BRICK_FACE_TO_VERT(1,Ifc);Nv2 = BRICK_FACE_TO_VERT(2,Ifc);
        Nv3 = BRICK_FACE_TO_VERT(3,Ifc);Nv4 = BRICK_FACE_TO_VERT(4,Ifc);
      case(MDLN,TETR)
        Nv1 = TETRA_FACE_TO_VERT(1,Ifc);Nv2 = TETRA_FACE_TO_VERT(2,Ifc);
        Nv3 = TETRA_FACE_TO_VERT(3,Ifc);Nv4 = TETRA_FACE_TO_VERT(1,Ifc);
      case(MDLD,PYRA)
        Nv1 = PYRAM_FACE_TO_VERT(1,Ifc);Nv2 = PYRAM_FACE_TO_VERT(2,Ifc);
        Nv3 = PYRAM_FACE_TO_VERT(3,Ifc);Nv4 = PYRAM_FACE_TO_VERT(4,Ifc);
      end select
!
   end subroutine face_to_vert
!
!----------------------------------------------------------------------
!
!  ...return edge numbers for a face
   subroutine face_to_edge(Ntype,Ifc, Ne1,Ne2,Ne3,Ne4)
!
      integer, intent(in)  :: Ntype,Ifc
      integer, intent(out) :: Ne1,Ne2,Ne3,Ne4
!
      select case(Ntype)
      case(MDLP,PRIS)
        Ne1 = PRISM_FACE_TO_EDGE(1,Ifc);Ne2 = PRISM_FACE_TO_EDGE(2,Ifc);
        Ne3 = PRISM_FACE_TO_EDGE(3,Ifc);Ne4 = PRISM_FACE_TO_EDGE(4,Ifc);
      case(MDLB,BRIC)
        Ne1 = BRICK_FACE_TO_EDGE(1,Ifc);Ne2 = BRICK_FACE_TO_EDGE(2,Ifc);
        Ne3 = BRICK_FACE_TO_EDGE(3,Ifc);Ne4 = BRICK_FACE_TO_EDGE(4,Ifc);
      case(MDLN,TETR)
        Ne1 = TETRA_FACE_TO_EDGE(1,Ifc);Ne2 = TETRA_FACE_TO_EDGE(2,Ifc);
        Ne3 = TETRA_FACE_TO_EDGE(3,Ifc);Ne4 = TETRA_FACE_TO_EDGE(1,Ifc);
      case(MDLD,PYRA)
        Ne1 = PYRAM_FACE_TO_EDGE(1,Ifc);Ne2 = PYRAM_FACE_TO_EDGE(2,Ifc);
        Ne3 = PYRAM_FACE_TO_EDGE(3,Ifc);Ne4 = PYRAM_FACE_TO_EDGE(4,Ifc);
      end select
!
   end subroutine face_to_edge
!
!----------------------------------------------------------------------
!
!> @param[in ] Ftype    - 2D element type = 3D element face type
!> @param[in ] Nface_or - the face orientation as seen by the element
!> @param[out] Nfver    - face edge numbers from the point of view
!                         of the face
!> @date       Feb 2023
!----------------------------------------------------------------------
   subroutine face_to_vert_nos(Ftype,Nface_or, Nfver)
!
      integer, intent(in ) :: Ftype,Nface_or
      integer, intent(out) :: Nfver(4)
!
      integer :: iv
!
      select case(Ftype)
      case(MDLT)
        do iv=1,3
          Nfver(iv) = TRIAN_L2G(iv,Nface_or)
        enddo
      case(MDLQ)
        do iv=1,4
          Nfver(iv) = QUADR_L2G(iv,Nface_or)
        enddo
      end select
!
   end subroutine face_to_vert_nos
!
!----------------------------------------------------------------------
!
   subroutine face_to_edge_nos(Ftype,Nface_or, Nfedg)
!
      integer, intent(in ) :: Ftype,Nface_or
      integer, intent(out) :: Nfedg(4)
!
      integer :: ie
!
      select case(Ftype)
      case(MDLT)
        do ie=1,3
          Nfedg(ie) = TRIAN_L2G_EDGE(ie,Nface_or)
        enddo
      case(MDLQ)
        do ie=1,4
          Nfedg(ie) = QUADR_L2G_EDGE(ie,Nface_or)
        enddo
      end select
!
   end subroutine face_to_edge_nos
!
!----------------------------------------------------------------------
!
!> @param[in ] Ftype    - 2D element type = 3D element face type
!> @param[in ] Nface_or - the face orientation as seen by the element
!> @param[in ] Nedge_or - orientations for the face edges from
!                         the element point of view
!> @param[out] Nfedg_or - orientations for the face edges from
!                         the face point of view
!> @date       Feb 2023
!----------------------------------------------------------------------
   subroutine face_to_edge_orient(Ftype,Nface_or,Nedge_or, Nfedg_or)
!
      integer, intent(in ) :: Ftype
      integer, intent(in ) :: Nface_or
      integer, intent(in ) :: Nedge_or(4)
      integer, intent(out) :: Nfedg_or(4)
!
      integer :: ie
!
      Nfedg_or(1:4) = 0
      select case(Ftype)
      case(MDLT)
        do ie=1,3
          Nfedg_or(ie) = &
                  mod(Nedge_or(ie)+TRIAN_ORIENT_L2G(ie,Nface_or),2)
        enddo
      case(MDLQ)
        do ie=1,4
          Nfedg_or(ie) = &
                  mod(Nedge_or(ie)+QUADR_ORIENT_L2G(ie,Nface_or),2)
        enddo
      end select
!
   end subroutine face_to_edge_orient
!
!----------------------------------------------------------------------
!
!  ...return local edge parametrizations......
   subroutine edge_param(Ntype,Ie,T, Xi,Dxidt)
!
      integer, intent(in)  :: Ntype,Ie
      real(8), intent(in)  :: T
      real(8), intent(out) :: Xi(*),Dxidt(*)
!
      real(8) :: xi1(3),xi2(3)
      integer :: n1,n2,nvar
!
      select case(Ntype)
      case(MDLT,TRIA)
        n1 = TRIAN_EDGE_TO_VERT(1,Ie);  n2 = TRIAN_EDGE_TO_VERT(2,Ie)
        xi1(1:2) = TRIAN_COORD(1:2,n1); xi2(1:2) = TRIAN_COORD(1:2,n2)
        nvar=2
      case(MDLQ,QUAD)
        n1 = QUADR_EDGE_TO_VERT(1,Ie);  n2 = QUADR_EDGE_TO_VERT(2,Ie)
        xi1(1:2) = QUADR_COORD(1:2,n1); xi2(1:2) = QUADR_COORD(1:2,n2)
        nvar=2
      case(MDLP,PRIS)
        n1 = PRISM_EDGE_TO_VERT(1,Ie);  n2 = PRISM_EDGE_TO_VERT(2,Ie)
        xi1(1:3) = PRISM_COORD(1:3,n1); xi2(1:3) = PRISM_COORD(1:3,n2)
        nvar=3
      case(MDLB,BRIC)
        n1 = BRICK_EDGE_TO_VERT(1,Ie);  n2 = BRICK_EDGE_TO_VERT(2,Ie)
        xi1(1:3) = BRICK_COORD(1:3,n1); xi2(1:3) = BRICK_COORD(1:3,n2)
        nvar=3
      case(MDLN,TETR)
        n1 = TETRA_EDGE_TO_VERT(1,Ie);  n2 = TETRA_EDGE_TO_VERT(2,Ie)
        xi1(1:3) = TETRA_COORD(1:3,n1); xi2(1:3) = TETRA_COORD(1:3,n2)
        nvar=3
      case(MDLD,PYRA)
        n1 = PYRAM_EDGE_TO_VERT(1,Ie);  n2 = PYRAM_EDGE_TO_VERT(2,Ie)
        xi1(1:3) = PYRAM_COORD(1:3,n1); xi2(1:3) = PYRAM_COORD(1:3,n2)
        nvar=3
      case default
        write(*,*) 'edge_param'; stop
      end select
!
      Dxidt(1:nvar) = xi2(1:nvar) - xi1(1:nvar)
      Xi(1:nvar) = xi1(1:nvar) + T*Dxidt(1:nvar)
!
   end subroutine edge_param
!
!----------------------------------------------------------------------
!
!  ...return local face parametrizations......
   subroutine face_param(Ntype,Iface,T, Xi,Dxidt)
!
      use control , only : GEOM_TOL
!
      integer, intent(in)  :: Ntype,Iface
      real(8), intent(in)  :: T(2)
      real(8), intent(out) :: Xi(3),Dxidt(3,2)
!
      real(8) :: xi1(3),xi2(3),xi3(3)
      integer :: n,k
!
!  ...check that point is inside master face
      if (face_type(Ntype,Iface).eq.TRIA) then
        if ((T(1)     .lt.-GEOM_TOL    ).or. &
            (T(2)     .lt.-GEOM_TOL    ).or. &
            (T(1)+T(2).gt.1.d0+GEOM_TOL)    ) then
          write(*,*)'face_param: point outside of tria face'
          write(*,1000) T(1:2)
 1000     format(' T = ',2(e12.5,1x))
          call pause
        endif
          k=3
      endif
!
      if (face_type(Ntype,Iface).eq.RECT) then
        if ((T(1).lt.-GEOM_TOL).or.(T(1).gt.1.d0+GEOM_TOL).or. &
            (T(2).lt.-GEOM_TOL).or.(T(2).gt.1.d0+GEOM_TOL)) then
          write(*,*)'face_param: point outside of rect face'
          write(*,1000) T(1:2)
          call pause
        endif
          k=4
      endif
!
      select case(Ntype)
      case(MDLP,PRIS)
        n = PRISM_FACE_TO_VERT(1,Iface); xi1(1:3) = PRISM_COORD(1:3,n)
        n = PRISM_FACE_TO_VERT(2,Iface); xi2(1:3) = PRISM_COORD(1:3,n)
        n = PRISM_FACE_TO_VERT(k,Iface); xi3(1:3) = PRISM_COORD(1:3,n)
      case(MDLB,BRIC)
        n = BRICK_FACE_TO_VERT(1,Iface); xi1(1:3) = BRICK_COORD(1:3,n)
        n = BRICK_FACE_TO_VERT(2,Iface); xi2(1:3) = BRICK_COORD(1:3,n)
        n = BRICK_FACE_TO_VERT(k,Iface); xi3(1:3) = BRICK_COORD(1:3,n)
      case(MDLN,TETR)
        n = TETRA_FACE_TO_VERT(1,Iface); xi1(1:3) = TETRA_COORD(1:3,n)
        n = TETRA_FACE_TO_VERT(2,Iface); xi2(1:3) = TETRA_COORD(1:3,n)
        n = TETRA_FACE_TO_VERT(k,Iface); xi3(1:3) = TETRA_COORD(1:3,n)
      case(MDLD,PYRA)
        n = PYRAM_FACE_TO_VERT(1,Iface); xi1(1:3) = PYRAM_COORD(1:3,n)
        n = PYRAM_FACE_TO_VERT(2,Iface); xi2(1:3) = PYRAM_COORD(1:3,n)
        n = PYRAM_FACE_TO_VERT(k,Iface); xi3(1:3) = PYRAM_COORD(1:3,n)
      end select
!
      Dxidt(1:3,1) = xi2(1:3)-xi1(1:3)
      Dxidt(1:3,2) = xi3(1:3)-xi1(1:3)
      Xi(1:3) = xi1(1:3) + T(1)*Dxidt(1:3,1) + T(2)*Dxidt(1:3,2)
!
   end subroutine face_param
!
!----------------------------------------------------------------------
!
!  ...return normal signs corresponding to face parametrizations......
      function Nsign_param(Ntype,Iface)
!
      integer Nsign_param
      integer Ntype,Iface
!
      select case(Ntype)
      case(MDLP,PRIS)
        select case(Iface)
        case(1,5)
          Nsign_param=-1
        case(2,3,4)
          Nsign_param= 1
        case default
          write(*,*) 'Nsign_param'; stop
        end select
      case(MDLB,BRIC)
        select case(Iface)
        case(1,5,6)
          Nsign_param=-1
        case(2,3,4)
          Nsign_param= 1
        case default
          write(*,*) 'Nsign_param'; stop
        end select
      case(MDLN,TETR)
        select case(Iface)
        case(1,4)
          Nsign_param=-1
        case(2,3)
          Nsign_param= 1
        case default
          write(*,*) 'Nsign_param'; stop
        end select
      case(MDLD,PYRA)
        select case(Iface)
        case(1,4,5)
          Nsign_param=-1
        case(2,3)
          Nsign_param= 1
        case default
          write(*,*) 'Nsign_param'; stop
        end select
      case default
        write(*,*) 'Nsign_param'; stop
      end select
!
      end function Nsign_param
!
!----------------------------------------------------------------------
!
!  ...check if element face is triangle or rectangle
      function Face_type(Ntype,Iface)
!
      integer :: Face_type
      integer :: Ntype,Iface
!
      select case(Ntype)
        case(MDLP,PRIS)
          select case(Iface)
            case(1,2)  ; Face_type = TRIA
            case(3,4,5); Face_type = RECT
            case default
               write(*,*) 'Face_type'; stop
          end select
        case(MDLB,BRIC)
          Face_type = RECT
        case(MDLN,TETR)
          Face_type = TRIA
        case(MDLD,PYRA)
          select case(Iface)
            case(1)      ; Face_type = RECT
            case(2,3,4,5); Face_type = TRIA
            case default
               write(*,*) 'Face_type'; stop
          end select
        case default
          write(*,*) 'Face_type'; stop
      end select
!
      end function Face_type
!
!----------------------------------------------------------------------
!
!  ...return local face node numbers
   subroutine face_nodes(Ntype,Iface, Nface_nodes,Nrfn)
!
      integer, intent(in)  :: Ntype,Iface
      integer, intent(out) :: Nface_nodes(9),Nrfn
!
      Nface_nodes(1:9) = 0
      select case(Ntype)
      case(MDLB,BRIC)
        Nface_nodes(1:4) = BRICK_FACE_TO_VERT(1:4,Iface)
        Nface_nodes(5:8) = 8+BRICK_FACE_TO_EDGE(1:4,Iface)
        Nface_nodes(9) = 20+Iface
        Nrfn=9
      case(MDLP,PRIS)
        select case(Iface)
        case(1,2)
          Nface_nodes(1:3) = PRISM_FACE_TO_VERT(1:3,Iface)
          Nface_nodes(4:6) = 6+PRISM_FACE_TO_EDGE(1:3,Iface)
          Nface_nodes(7) = 15+Iface
          Nrfn=7
        case(3,4,5)
          Nface_nodes(1:4) = PRISM_FACE_TO_VERT(1:4,Iface)
          Nface_nodes(5:8) = 6+PRISM_FACE_TO_EDGE(1:4,Iface)
          Nface_nodes(9) = 15+Iface
          Nrfn=9
        end select
      case(MDLN,TETR)
        Nface_nodes(1:3) = TETRA_FACE_TO_VERT(1:3,Iface)
        Nface_nodes(4:6) = 4+TETRA_FACE_TO_EDGE(1:3,Iface)
        Nface_nodes(7) = 10+Iface
        Nrfn=7
      case(MDLD,PYRA)
        select case(Iface)
        case(1)
          Nface_nodes(1:4) = PYRAM_FACE_TO_VERT(1:4,Iface)
          Nface_nodes(5:8) = 5+PYRAM_FACE_TO_EDGE(1:4,Iface)
          Nface_nodes(9) = 13+Iface
          Nrfn=9
        case(2,3,4,5)
          Nface_nodes(1:3) = PYRAM_FACE_TO_VERT(1:3,Iface)
          Nface_nodes(4:6) = 5+PYRAM_FACE_TO_EDGE(1:3,Iface)
          Nface_nodes(7) = 13+Iface
          Nrfn=7
        end select
      end select
!
   end subroutine face_nodes
!
!----------------------------------------------------------------------
!
!  ...return local face order of approximation
   subroutine face_order(Ntype,Iface,Norder, Nordf)
!
      integer, intent(in)  :: Ntype,Iface,Norder(19)
      integer, intent(out) :: Nordf(5)
!
      integer :: i,k
!
      Nordf(1:5) = 0
      select case(Ntype)
      case(MDLP,PRIS)
        select case(Iface)
        case(1,2)
          do i=1,3
            k = PRISM_FACE_TO_EDGE(i,Iface)
            Nordf(i) = Norder(k)
          enddo
          Nordf(4) = Norder(9+Iface)
        case(3,4,5)
          do i=1,4
            k = PRISM_FACE_TO_EDGE(i,Iface)
            Nordf(i) = Norder(k)
          enddo
          Nordf(5) = Norder(9+Iface)
        end select
      case(MDLB,BRIC)
        do i=1,4
          k = BRICK_FACE_TO_EDGE(i,Iface)
          Nordf(i) = Norder(k)
        enddo
        Nordf(5) = Norder(12+Iface)
      case(MDLN,TETR)
        do i=1,3
          k = TETRA_FACE_TO_EDGE(i,Iface)
          Nordf(i) = Norder(k)
        enddo
        Nordf(4) = Norder(6+Iface)
      case(MDLD,PYRA)
        select case(Iface)
        case(1)
          do i=1,4
            k = PYRAM_FACE_TO_EDGE(i,Iface)
            Nordf(i) = Norder(k)
          enddo
          Nordf(5) = Norder(8+Iface)
        case(2,3,4,5)
          do i=1,3
            k = PYRAM_FACE_TO_EDGE(i,Iface)
            Nordf(i) = Norder(k)
          enddo
          Nordf(4) = Norder(8+Iface)
        end select
      end select
!
   end subroutine face_order
!
!----------------------------------------------------------------------
!
!  ...return number of dof for a single component
   subroutine ndof_nod(Ntype,Nord, NdofH,NdofE,NdofV,NdofQ)
!
      integer, intent(in)  :: Ntype,Nord
      integer, intent(out) :: NdofH,NdofE,NdofV,NdofQ
      integer :: nordx,nordy,nordz,naux,n
!
      select case(Ntype)
!  ...VERTEX
      case(VERT)
        NdofH=1; NdofE=0; NdofV=0; NdofQ=0
!  ...EDGE
      case(MEDG)
        NdofH=Nord-1; NdofE=Nord; NdofV=0; NdofQ=0
!  ...TRIANGLE
      case(TRIA,MDLT)
        NdofH=(Nord-2)*(Nord-1)/2
        NdofE=(Nord-1)* Nord
        NdofV= Nord   *(Nord+1)/2
        NdofQ=0
!  ...QUAD
      case(RECT,MDLQ)
        call decode(Nord, nordx,nordy)
        NdofH=(nordx-1)*(nordy-1)
        NdofE= nordx   *(nordy-1) + &
              (nordx-1)* nordy
        NdofV= nordx   * nordy
        NdofQ=0
!  ...PRISM
      case(PRIS,MDLP)
        call decode(Nord, nordx,nordz)
        NdofH=(nordx-2)*(nordx-1)/2*(nordz-1)
        NdofE=(nordx-1)* nordx     *(nordz-1) + &
              (nordx-2)*(nordx-1)/2* nordz
        NdofV=(nordx-1)* nordx     * nordz    + &
               nordx   *(nordx+1)/2*(nordz-1)
        NdofQ=(nordx+1)* nordx   /2* nordz
!  ...HEXA
      case(BRIC,MDLB)
        call decode(Nord, naux,nordz)
        call decode(naux, nordx,nordy)
        NdofH=(nordx-1)*(nordy-1)*(nordz-1)
        NdofE= nordx   *(nordy-1)*(nordz-1) + &
              (nordx-1)* nordy   *(nordz-1) + &
              (nordx-1)*(nordy-1)* nordz
        NdofV=(nordx-1)* nordy   * nordz    + &
               nordx   *(nordy-1)* nordz    + &
               nordx   * nordy   *(nordz-1)
        NdofQ= nordx   * nordy   * nordz
!  ...TETRA
      case(TETR,MDLN)
        NdofH=(Nord-3)*(Nord-2)*(Nord-1)/6
        NdofE=(Nord-2)*(Nord-1)* Nord   /2
        NdofV=(Nord-1)* Nord   *(Nord+1)/2
        NdofQ= Nord   *(Nord+1)*(Nord+2)/6
!  ...PYRAMID
      case(PYRA,MDLD)
        n = Nord-1
!  LD: value for the original pyramid of Paolo
!!!        NdofH=(2*n**3+3*n**2+n)/6
        NdofH=  (Nord-1)**3
        NdofE=3*(Nord-1)**2*Nord
!        NdofV=2*(Nord-1)*Nord**2+ Nord**3
        NdofV=3*(Nord-1)   *Nord**2
!        NdofQ=Nord**2*(Nord+1)
        NdofQ=   Nord**3

      case default
        write(*,*) 'Error!!! ndof_nod: Type = ', S_Type(Ntype)
        NdofH=0; NdofE=0; NdofV=0; NdofQ=0
        stop
      end select
!!!      write(*,*) 'ndof_nod: S_Type(NType),Nord, NdofH,NdofE,NdofV,NdofQ = ', &
!!!                            S_Type(NType),Nord, NdofH,NdofE,NdofV,NdofQ
!!!      call pause
!
   end subroutine ndof_nod
!
!
   end module element_data
!
!
!======================================================================
!  REMARK : the following subroutines are somehow related to module   |
!           "element_data". They are (meant to be) robust wrt the     |
!           info contained in "element_data".                         |
!======================================================================
!
!
!
!----------------------------------------------------------------------
!> @brief      given a point in a 2D master element having coordinates
!              T(1:2) wrt a LOCAL system of coordinates, routine
!              returns coordinates Eta(1:2) wrt a GLOBAL system having
!              orientation Norient wrt the LOCAL system. Orientations
!              are defined by TRIAN,QUADR_L2G.
!
!> @param[in ] Ntype   - 2D element type
!> @param[in ] T       - coordinates in the LOCAL system
!> @param[in ] Norient - orientation of GLOBAL system wrt LOCAL system
!> @param[out] Eta     - coordinates in the GLOBAL system
!> @param[out] Detadt  - derivatives
!
!> @date       Feb 2023
!----------------------------------------------------------------------
   subroutine local2global(Ntype,T,Norient, Eta,Detadt)
!
      use element_data , only : TRIAN_L2G,QUADR_L2G,QUADR_EDGE_TO_VERT
      use node_types
!
      implicit none
      integer               ,intent(in ) :: Ntype
      real(8),dimension(2  ),intent(in ) :: T
      integer               ,intent(in ) :: Norient
      real(8),dimension(2  ),intent(out) :: Eta
      real(8),dimension(2,2),intent(out) :: Detadt
!
      real(8),dimension(  4) :: rlam
      real(8),dimension(2,4) :: drlam
      integer :: iv,iv1,iv2,iv3,iv4,iedge,i
!
!----------------------------------------------------------------------
!
!  ...initialize
      Eta(1:2)=0.d0 ; Detadt(1:2,1:2)=0.d0
!
!  ...vertex shape functions
      call vshape2(Ntype,T, rlam,drlam)
!
      select case(Ntype)
!
!  ...TRIANGLE.........................................................
      case(TRIA)
        iv2=TRIAN_L2G(2,Norient) ; iv3=TRIAN_L2G(3,Norient)
        Eta(1)=rlam(iv2) ; Detadt(1,1:2)=drlam(1:2,iv2)
        Eta(2)=rlam(iv3) ; Detadt(2,1:2)=drlam(1:2,iv3)
!
!  ...QUAD.............................................................
      case(QUAD)
        iv1=QUADR_L2G(1,Norient)
        iv2=QUADR_L2G(2,Norient)
        iv4=QUADR_L2G(4,Norient)
!
!  .....projection on (iv1,iv2) edge
        call quad_aux(iv1,iv2, iedge)
        do i=1,2
          iv=QUADR_EDGE_TO_VERT(i,iedge)
          Eta   (1)    =Eta   (1)    +rlam     (iv)
          Detadt(1,1:2)=Detadt(1,1:2)+drlam(1:2,iv)
        enddo
!  .....projection on (iv1,iv4) edge
        call quad_aux(iv1,iv4, iedge)
        do i=1,2
          iv=QUADR_EDGE_TO_VERT(i,iedge)
          Eta(   2    )=Eta(   2    )+rlam(     iv)
          Detadt(2,1:2)=Detadt(2,1:2)+drlam(1:2,iv)
        enddo
!
      case default
        write(*,7001) S_Type(Ntype)
 7001   format(' local2global: invalid Type = ',a10)
        stop
      end select
!
   end subroutine local2global
!
!----------------------------------------------------------------------
!> @brief      auxiliary routine needed by "hp3gen"
!
!> @param[in ] Ntype   - element type
!> @param[in ] Ivert   - vertex number
!> @param[out] Nrfaces - number of adjacent faces
!> @param[out] Nofaces - adjacent face numbers
!
!> @date       Feb 2023
!----------------------------------------------------------------------
!
   subroutine vert_to_faces(Ntype,Ivert, Nrfaces,Nofaces)
!
      use element_data , only : PRISM_VERT_TO_FACE, &
                                BRICK_VERT_TO_FACE, &
                                TETRA_VERT_TO_FACE, &
                                PYRAM_VERT_TO_FACE
      use node_types
!
      implicit none
      integer, intent(in ) :: Ntype
      integer, intent(in ) :: Ivert
      integer, intent(out) :: Nrfaces
      integer, intent(out) :: Nofaces(4)
!
!  ...number of connected faces
      Nrfaces=3
      select case(Ntype)
      case(BRIC,MDLB)
        Nofaces(1:3)=BRICK_VERT_TO_FACE(1:3,Ivert)
      case(TETR,MDLN)
        Nofaces(1:3)=TETRA_VERT_TO_FACE(1:3,Ivert)
      case(PRIS,MDLP)
        Nofaces(1:3)=PRISM_VERT_TO_FACE(1:3,Ivert)
      case(PYRA,MDLD)
        select case(Ivert)
!  .....lateral triangular faces
        case(1,2,3,4)
          Nofaces(1:3)=PYRAM_VERT_TO_FACE(1:3,Ivert)
!  .....rectangular base
        case(5)
          Nrfaces=4
          Nofaces(1:4)=PYRAM_VERT_TO_FACE(1:4,Ivert)
        end select
      end select
!
   end subroutine vert_to_faces
!
!----------------------------------------------------------------------
!> @brief      auxiliary routine needed by "hp3gen"
!
!> @param[in ] Type    - element type
!> @param[in ] Iedge   - edge number
!> @param[out] Nofaces - adjacent face numbers
!
!> @date       Feb 2023
!----------------------------------------------------------------------
!
   subroutine edge_to_faces(Ntype,Iedge, Nofaces)
!
      use element_data , only : PRISM_EDGE_TO_FACE, &
                                BRICK_EDGE_TO_FACE, &
                                TETRA_EDGE_TO_FACE, &
                                PYRAM_EDGE_TO_FACE
      use node_types
!
      implicit none
      integer             ,intent(in ) :: Ntype
      integer             ,intent(in ) :: Iedge
      integer,dimension(2),intent(out) :: Nofaces
!
      select case(Ntype)
      case(BRIC,MDLB)
        Nofaces(1:2)=BRICK_EDGE_TO_FACE(1:2,Iedge)
      case(TETR,MDLN)
        Nofaces(1:2)=TETRA_EDGE_TO_FACE(1:2,Iedge)
      case(PRIS,MDLP)
        Nofaces(1:2)=PRISM_EDGE_TO_FACE(1:2,Iedge)
      case(PYRA,MDLD)
        Nofaces(1:2)=PYRAM_EDGE_TO_FACE(1:2,Iedge)
      end select
!
   end subroutine edge_to_faces
!
!----------------------------------------------------------------------
!
   subroutine quad_aux(Iv1,Iv2, Iedge)
!
      use element_data , only : QUADR_EDGE_TO_VERT
!
      implicit none
      integer,intent(in ) :: Iv1,Iv2
      integer,intent(out) :: Iedge
!
      integer,save,dimension(4,4) :: AUX
!$OMP THREADPRIVATE(AUX)
      logical,save                :: INITIALIZED=.false.
!$OMP THREADPRIVATE(INITIALIZED)
      integer,       dimension(2) :: iedge1,iedge2,list
      integer                     :: i,j,i1,i2,ifound
!----------------------------------------------------------------------
!
      if (.not.INITIALIZED) then
!
!  ...initialize
      AUX(1:4,1:4)=0
!
!  ...double loop over edges
      do i=1,4
!
!  .....1st edge endpoints
        i1=QUADR_EDGE_TO_VERT(1,i)
        i2=QUADR_EDGE_TO_VERT(2,i)
        iedge1(1)=i1 ; iedge1(2)=i2
!
        do j=1,4
!
!  .......2nd edge endpoints
          iedge2(1)=QUADR_EDGE_TO_VERT(1,j)
          iedge2(2)=QUADR_EDGE_TO_VERT(2,j)
!
!  .......compare and record
          call list_minus(iedge1(1:2),iedge2(1:2),2, list,ifound)
          if (ifound.ne.1) cycle
!
          if (list(1).eq.i1)  AUX(i1,i2)=j
          if (list(1).eq.i2)  AUX(i2,i1)=j
        enddo
      enddo
!
!!!  ...printing
!!      do i=1,4 ; do j=1,4
!!!  .....print only relevant entries
!!        if (AUX(i,j).eq.0) cycle
!!        write(*,1000) i,j,AUX(i,j)
!! 1000   format(' iv1,2, edge = '2(i1,2x),2x,i1)
!!      enddo    ; enddo
!
!  ...update initialization flag
      INITIALIZED=.true.
!
      endif
!
!  ...retrieve edge number
      Iedge=AUX(Iv1,Iv2)
!
   end subroutine quad_aux
!
!----------------------------------------------------------------------
!> @brief      Returns {List3} = {List1} - {List2}
!
!> @param[in ] List1,2 - lists of integers
!> @param[in ] N       - length of List1,2
!> @param[out] List3   - difference
!> @param[out] J       - length of List3
!
!> @date       Feb 2023
!----------------------------------------------------------------------
   subroutine list_minus(List1,List2,N, List3,J)
!
      implicit none
      integer             ,intent(in ) :: N
      integer,dimension(N),intent(in ) :: List1,List2
      integer,dimension(N),intent(out) :: List3
      integer             ,intent(out) :: J
      integer :: i,ifound
!
!  ...initialize
      J=0 ; List3(1:N)=0
!
!  ...loop over entries of List1
      do i=1,N
!  .....locate entry of List1 on List2
        call locate(List1(i),List2(1:N),N, ifound)
!  .....if entry was not found, add it to List3
        if (ifound.eq.0) then
          J=J+1 ; List3(J)=List1(i)
        endif
      enddo
!
   end subroutine list_minus
!
!----------------------------------------------------------------------
!
   subroutine hexa_aux(Iv1,Iv2, Iface)
!
      use element_data , only : BRICK_EDGE_TO_VERT, &
                                BRICK_VERT_TO_FACE
!
      implicit none
      integer,intent(in ) :: Iv1,Iv2
      integer,intent(out) :: Iface
!
      integer,save,dimension(8,8) :: AUX
!$OMP THREADPRIVATE(AUX)
      logical,save                :: INITIALIZED=.false.
!$OMP THREADPRIVATE(INITIALIZED)
      integer,       dimension(3) :: iface1,iface2,list
      integer                     :: i,i1,i2,ifound
#if DEBUG_FLAG
      integer :: j,iprint=0
#endif
!
      if (.not.INITIALIZED) then
!
!  ...initialize
      AUX(1:8,1:8)=0
!
!  ...loop over edges
      do i=1,12
        i1=BRICK_EDGE_TO_VERT(1,i)
        i2=BRICK_EDGE_TO_VERT(2,i)
!
!  .....faces connected to 1st and 2nd vertex
        iface1(1:3)=BRICK_VERT_TO_FACE(1:3,i1)
        iface2(1:3)=BRICK_VERT_TO_FACE(1:3,i2)
!
!  .....determine non-common face, and store
        call list_minus(iface2(1:3),iface1(1:3),3, list,ifound)
        if (ifound.eq.0) then
          write(*,*)'hexa_aux: no face found!'
          stop
        endif
        AUX(i1,i2)=list(1)
!
        call list_minus(iface1(1:3),iface2(1:3),3, list,ifound)
        if (ifound.eq.0) then
          write(*,*)'hexa_aux: no face found!'
          stop
        endif
        AUX(i2,i1)=list(1)
      enddo
!
#if DEBUG_FLAG
!  ...printing
      do i=1,8 ; do j=1,8
!  .....print only relevant entries
        if (AUX(i,j).eq.0) cycle
        write(*,1000) i,j,AUX(i,j)
 1000   format(' iv1,2, face = '2(i1,2x),2x,i1)
      enddo    ; enddo
#endif
!
!  ...update initialization flag
      INITIALIZED=.true.
!
      endif
!
!  ...retrieve face number
      Iface=AUX(Iv1,Iv2)
!
   end subroutine hexa_aux
