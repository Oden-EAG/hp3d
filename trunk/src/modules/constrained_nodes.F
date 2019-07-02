c----------------------------------------------------------------------
c
c   module name        - constrained_nodes
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - module contains information on constrained
c                        nodes for an element (a local data base)
c
c
c----------------------------------------------------------------------
c
      module constrained_nodes
!      save
c
c  ...a flag indicating whether the info on constraints is to be
c     collected
      integer :: INFO_CONSTRAINTS
!$OMP THREADPRIVATE (INFO_CONSTRAINTS)
c
c  ...information saved from routine elem_nodes
      integer :: FATH_NODES(27), FATH_ORIENT(27), SON_NUM
!$OMP THREADPRIVATE (FATH_NODES, FATH_ORIENT, SON_NUM)
      character(len=4) :: FATH_TYPE
!$OMP THREADPRIVATE (FATH_TYPE)
c
c  ...information on constrained nodes
      integer :: NODES_CONSTR(27)
!$OMP THREADPRIVATE (NODES_CONSTR)
c
c
c     a constrained node is identified with a nickname
c     nick = constraining edge or face number*100 + case number
c
c  ...nodes constrained by an edge
c     case 11: first mid-edge node constrained by an edge
c     case 12: second mid-edge node constrained by an edge
c     case 13: vertex node constrained by an edge
c
c  ...nodes constrained by an h4-refined quad face
c     case 21: south-west mid-face node
c     case 22: south-east mid-face node
c     case 23: north-east mid-face node
c     case 24: north-west mid-face node
c     case 25: south mid-edge node
c     case 26: east mid-edge node
c     case 27: north mid-edge node
c     case 28: west mid-edge node
c     case 29: vertex node
c
c  ...nodes constrained by the west half of an h2/h2 refined quad face
c     case 31: south mid-face node
c     case 32: north mid-face node
c     case 33: horizontal mid-edge node
c
c  ...nodes constrained by the east half of an h2/h2 refined quad face
c     case 34: south mid-face node
c     case 35: north mid-face node
c     case 36: horizontal mid-edge node
c
c  ...nodes constrained by the vertical edge of an h2/h2 refined quad face
c     case 37: south mid-edge node
c     case 38: north mid-edge node
c     case 39: vertex node
c
c  ...nodes constrained by the south half of an h2/h2 refined quad face
c     case 41: west mid-face node
c     case 42: east mid-face node
c     case 43: vertical mid-edge node
c
c  ...nodes constrained by the north half of an h2/h2 refined quad face
c     case 44: west mid-face node
c     case 45: east mid-face node
c     case 46: vertical mid-edge node
c
c  ...nodes constrained by the vertical edge of an h2/h2 refined quad face
c     case 47: west mid-edge node
c     case 48: east mid-edge node
c     case 49: vertex node
c
c  ...nodes constrained by a vertically h2-refined quad face
c     case 51: west mid-face node
c     case 52: east mid-face node
c     case 53: vertical mid-edge node
c
c  ...nodes constrained by a horizontally h2-refined quad face
c     case 61: south mid-face node
c     case 62: north mid-face node
c     case 63: horizontal mid-edge node
c
c  ...nodes constrained by a triangular face
c     case 71-74: mdlt nodes
c     case 75-77: medg nodes
c
c  ...mdlq nodes constrained by a refined triangular face
c     case 82: kref=2
c     case 83: kref=3
c     case 84: kref=4
c
c  ...data base for constraining edges
      integer, parameter :: MAXNRE=20
      integer :: NR_EDGES,NEDGC(MAXNRE),NEDG_CONS(2,MAXNRE)
!$OMP THREADPRIVATE (NR_EDGES, NEDGC, NEDG_CONS)
c
c
c     NR_EDGES - number of constraining edges in the data base
c     NEDGC(*) - constraining mid-edge node number
c     NEDG_CONS(1:2,*) - vertex nodes on the edge
c
c  ...data base for constraining faces
      integer, parameter :: MAXNRF=20
      integer :: NR_FACES,NFACEC(MAXNRF),NFACE_CONS(8,MAXNRF)
!$OMP THREADPRIVATE (NR_FACES, NFACEC, NFACE_CONS)
c
c
c     NR_FACES - number of constraining faces in the data base
c     NFACEC(*) - constraining mid-face node number
c     NFACE_CONS(1:4,*)   - mid-edge nodes on the face with sign
c                           indicating orientation
c     NFACE_CONS(5:8,*)   - vertex nodes on the face
c
      contains
c
c  ...rotate vertex nodes on a constraining edge to fit edge global
c     coordinate
      subroutine rotate_edge_nodes(Norient,J)
      integer Norient,J,nloc(2)
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
      select case(Norient)
      case(0)
      case(1)
        nloc(1:2) = NEDG_CONS(1:2,J)
        do i=1,2
          i1 = imod(i+1,2)
          NEDG_CONS(i1,J) = nloc(i)
        enddo
      end select
      end subroutine rotate_edge_nodes
c
c  ...rotate edge and vertex nodes on a constraining triangular face
c     to fit face global coordinates
      subroutine rotate_trian_nodes(Norient,J)
      integer Norient,J,nloc(8),nedg(3,0:5),nvrt(3,0:5)
      data nedg/ 1, 2, 3,  2,-3,-1,  -3, 1,-2,
     .           3,-2, 1, -1, 3, 2,  -2,-1,-3/
      data nvrt/1,2,3, 2,3,1, 3,1,2, 1,3,2, 2,1,3, 3,2,1/
      nloc(1:8) = NFACE_CONS(1:8,J)
      do i=1,3
        i1 = iabs(nedg(i,Norient))
        NFACE_CONS(i,J) = nloc(i1)*isgn(nedg(i,Norient))
        i1 = nvrt(i,Norient)
        NFACE_CONS(4+i,J) = nloc(4+i1)
      enddo
      end subroutine rotate_trian_nodes
c
c  ...rotate edge and vertex nodes on a constraining rectangular face
c     to fit face global coordinates
      subroutine rotate_quadr_nodes(Norient,J)
      integer Norient,J,nloc(8),nedg(4,0:7),nvrq(4,0:7)
      data nedg/ 1, 2, 3, 4,  2,-3, 4,-1,  -3,-4,-1,-2, -4, 1,-2, 3,
     .           4, 3, 2, 1, -1, 4,-3, 2,  -2,-1,-4,-3,  3,-2, 1,-4/
      data nvrq/1,2,3,4, 2,3,4,1, 3,4,1,2, 4,1,2,3,
     .          1,4,3,2, 2,1,4,3, 3,2,1,4, 4,3,2,1/
      nloc(1:8) = NFACE_CONS(1:8,J)
      do i=1,4
        i1 = iabs(nedg(i,Norient))
        NFACE_CONS(i,J) = nloc(i1)*isgn(nedg(i,Norient))
        i1 = nvrq(i,Norient)
        NFACE_CONS(4+i,J) = nloc(4+i1)
      enddo
      end subroutine rotate_quadr_nodes
c
      end module constrained_nodes
