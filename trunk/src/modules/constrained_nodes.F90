!----------------------------------------------------------------------
!
!   module name        - constrained_nodes
!
!----------------------------------------------------------------------
!> @brief   module contains information on constrained
!!          nodes for an element (a local data base)
!> @date    Feb 2023
!----------------------------------------------------------------------
module constrained_nodes
!
   implicit none
!
!..flag indicating whether the info on constraints is to be collected
   integer :: INFO_CONSTRAINTS
!$OMP THREADPRIVATE (INFO_CONSTRAINTS)
!
!..information saved from routine elem_nodes
   integer :: FATH_NODES(27), FATH_ORIENT(27), FATH_TYPE, SON_NUM
!$OMP THREADPRIVATE (FATH_NODES, FATH_ORIENT, FATH_TYPE, SON_NUM)
!
!..information on constrained nodes
   integer :: NODES_CONSTR(27)
!$OMP THREADPRIVATE (NODES_CONSTR)
!
!
!     a constrained node is identified with a nickname
!     nick = constraining edge or face number*100 + case number
!
!  ...nodes constrained by an edge
!     case 11: first mid-edge node constrained by an edge
!     case 12: second mid-edge node constrained by an edge
!     case 13: vertex node constrained by an edge
!
!  ...nodes constrained by an h4-refined quad face
!     case 21: south-west mid-face node
!     case 22: south-east mid-face node
!     case 23: north-east mid-face node
!     case 24: north-west mid-face node
!     case 25: south mid-edge node
!     case 26: east mid-edge node
!     case 27: north mid-edge node
!     case 28: west mid-edge node
!     case 29: vertex node
!
!  ...nodes constrained by the west half of an h2/h2 refined quad face
!     case 31: south mid-face node
!     case 32: north mid-face node
!     case 33: horizontal mid-edge node
!
!  ...nodes constrained by the east half of an h2/h2 refined quad face
!     case 34: south mid-face node
!     case 35: north mid-face node
!     case 36: horizontal mid-edge node
!
!  ...nodes constrained by the vertical edge of an h2/h2 refined quad face
!     case 37: south mid-edge node
!     case 38: north mid-edge node
!     case 39: vertex node
!
!  ...nodes constrained by the south half of an h2/h2 refined quad face
!     case 41: west mid-face node
!     case 42: east mid-face node
!     case 43: vertical mid-edge node
!
!  ...nodes constrained by the north half of an h2/h2 refined quad face
!     case 44: west mid-face node
!     case 45: east mid-face node
!     case 46: vertical mid-edge node
!
!  ...nodes constrained by the vertical edge of an h2/h2 refined quad face
!     case 47: west mid-edge node
!     case 48: east mid-edge node
!     case 49: vertex node
!
!  ...nodes constrained by a vertically h2-refined quad face
!     case 51: west mid-face node
!     case 52: east mid-face node
!     case 53: vertical mid-edge node
!
!  ...nodes constrained by a horizontally h2-refined quad face
!     case 61: south mid-face node
!     case 62: north mid-face node
!     case 63: horizontal mid-edge node
!
!  ...nodes constrained by a triangular face
!     case 71-74: mdlt nodes
!     case 75-77: medg nodes
!
!  ...mdlq nodes constrained by a refined triangular face
!     case 82: kref=2
!     case 83: kref=3
!     case 84: kref=4
!
!..data base for constraining edges
   integer, parameter :: MAXNRE=20
   integer :: NR_EDGES,NEDGC(MAXNRE),NEDG_CONS(2,MAXNRE)
!$OMP THREADPRIVATE (NR_EDGES, NEDGC, NEDG_CONS)
!
!
!     NR_EDGES - number of constraining edges in the data base
!     NEDGC(*) - constraining mid-edge node number
!     NEDG_CONS(1:2,*) - vertex nodes on the edge
!
!  ...data base for constraining faces
   integer, parameter :: MAXNRF=20
   integer :: NR_FACES,NFACEC(MAXNRF),NFACE_CONS(8,MAXNRF)
!$OMP THREADPRIVATE (NR_FACES, NFACEC, NFACE_CONS)
!
!
!     NR_FACES - number of constraining faces in the data base
!     NFACEC(*) - constraining mid-face node number
!     NFACE_CONS(1:4,*)   - mid-edge nodes on the face with sign
!                           indicating orientation
!     NFACE_CONS(5:8,*)   - vertex nodes on the face
!
   contains
!
!----------------------------------------------------------------------
!>@brief rotate vertex nodes on a constraining edge to fit edge
!!       global coordinate
!>@date  Feb 2023
      subroutine rotate_edge_nodes(Norient,J)
         integer,intent(in) :: Norient,J
         integer :: nloc(2)
         integer :: i,i1
         integer, external :: imod
         select case(Norient)
         case(1)
           nloc(1:2) = NEDG_CONS(1:2,J)
           do i=1,2
             i1 = imod(i+1,2)
             NEDG_CONS(i1,J) = nloc(i)
           enddo
         case(0)
         end select
      end subroutine rotate_edge_nodes
!
!----------------------------------------------------------------------
!>@brief rotate edge and vertex nodes on a constraining triangular face
!!       to fit face global coordinates
!>@date  Feb 2023
      subroutine rotate_trian_nodes(Norient,J)
         integer, intent(in) :: Norient,J
         integer :: nloc(8),nedg(3,0:5),nvrt(3,0:5)
         integer :: i,i1
         integer, external :: isgn
         data nedg/ 1, 2, 3,  2,-3,-1,  -3, 1,-2, &
                    3,-2, 1, -1, 3, 2,  -2,-1,-3 /
         data nvrt/ 1,2,3, 2,3,1, 3,1,2, 1,3,2, 2,1,3, 3,2,1 /
         nloc(1:8) = NFACE_CONS(1:8,J)
         do i=1,3
           i1 = iabs(nedg(i,Norient))
           NFACE_CONS(i,J) = nloc(i1)*isgn(nedg(i,Norient))
           i1 = nvrt(i,Norient)
           NFACE_CONS(4+i,J) = nloc(4+i1)
         enddo
      end subroutine rotate_trian_nodes
!
!----------------------------------------------------------------------
!>@brief rotate edge and vertex nodes on a constraining rectangular face
!!       to fit face global coordinates
!>@date  Feb 2023
      subroutine rotate_quadr_nodes(Norient,J)
         integer, intent(in) :: Norient,J
         integer :: nloc(8),nedg(4,0:7),nvrq(4,0:7)
         integer :: i,i1
         integer, external :: isgn
         data nedg/ 1, 2, 3, 4,  2,-3, 4,-1,  -3,-4,-1,-2, -4, 1,-2, 3, &
                    4, 3, 2, 1, -1, 4,-3, 2,  -2,-1,-4,-3,  3,-2, 1,-4 /
         data nvrq/ 1,2,3,4, 2,3,4,1, 3,4,1,2, 4,1,2,3, &
                    1,4,3,2, 2,1,4,3, 3,2,1,4, 4,3,2,1 /
         nloc(1:8) = NFACE_CONS(1:8,J)
         do i=1,4
           i1 = iabs(nedg(i,Norient))
           NFACE_CONS(i,J) = nloc(i1)*isgn(nedg(i,Norient))
           i1 = nvrq(i,Norient)
           NFACE_CONS(4+i,J) = nloc(4+i1)
         enddo
      end subroutine rotate_quadr_nodes
!
end module constrained_nodes
