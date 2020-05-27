!----------------------------------------------------------------------
!
!   module name        - module 2Drefinements
!
!----------------------------------------------------------------------
!
!   latest revision    - Jul 01
!
!   purpose            - module contains logical information
!                        related to 2D h-refinements
!
!
!----------------------------------------------------------------------
!
      module refinements_2D
!
!  ...denumeration of sons of an edge, in a local, edge system
!     of coordinates
!
!     --1--3--2-->
!
!  ...denumeration of sons of an h11-refined quad middle node,
!     in the local system of coordinates
!
!     -----------------
!     |       |       |
!     |   4   7   3   |
!     |       |       |
!     ----8---9---6----
!     |       |       |
!     |   1   5   2   |
!     |       |       |
!     -----------------
!
!
!  ...denumeration of sons of an h10 and h01-refined quad middle
!     nodes in the local system of coordinates
!
!     -----------------     -----------------
!     |       |       |     |               |
!     |       |       |     |       2       |
!     |       |       |     |               |
!     |   1   3   2   |     --------3--------
!     |       |       |     |               |
!     |       |       |     |       1       |
!     |       |       |     |               |
!     -----------------     -----------------
!
!
!  ...denumeration of sons of a triangle middle node, in the local
!     system of coordinates
!
!     .
!     . .
!     .   .
!     . 3   .
!     ....7....
!     . .     . .
!     .   5 4 6   .
!     . 1   . . 2   .
!     .................
!
!  ...h4 refinement of a triangle...
      integer  :: PARENT_NODES_T4(2,7,4) = reshape( &
        (/ 0,0,  4,3,  6,3,  4,1,  7,5,  6,1,  7,1, &
           4,3,  0,0,  5,3,  4,2,  5,1,  7,6,  7,2, &
           6,3,  5,3,  0,0,  7,7,  5,2,  6,2,  7,3, &
           5,3,  6,3,  4,3,  7,7,  7,5,  7,6,  7,4 /),  &
        (/2,7,4/))
!
!  ...h4 refinement of a quad...
      integer  :: PARENT_NODES_Q11(2,9,4) = reshape( &
        (/ 0,0,  5,3,  9,9,  8,3,  5,1,  9,5,  9,8,  8,1,  9,1,  &
           5,3,  0,0,  6,3,  9,9,  5,2,  6,1,  9,6,  9,5,  9,2, &
            9,9,  6,3,  0,0,  7,3,  9,6,  6,2,  7,2,  9,7,  9,3, &
            8,3,  9,9,  7,3,  0,0,  9,8,  9,7,  7,1,  8,2,  9,4/), &
         (/2,9,4/))
!
!  ...h2 refinements of a quad...
      integer  :: PARENT_NODES_Q10(2,9,2) = reshape( &
         (/ 0,0,  5,3,  7,3,  0,0,  5,1,  9,3,  7,1,  0,0,  9,1, &
            5,3,  0,0,  0,0,  7,3,  5,2,  0,0,  7,2,  9,3,  9,2/), &
         (/2,9,2/))
      integer  :: PARENT_NODES_Q01(2,9,2) = reshape( &
         (/ 0,0,  0,0,  6,3,  8,3,  0,0,  6,1,  9,3,  8,1,  9,1, &
            8,3,  6,3,  0,0,  0,0,  9,3,  6,2,  0,0,  8,2,  9,2/), &
         (/2,9,2/))
!
!  ...for j-th node of i-th son
!     PARENT_NODES_???(1,j,i) = number of the parent node of the father
!     PARENT_NODES_???(2,j,i) = the corresponding son's node number
!
      integer  :: T4_ORIENT(3,4) = reshape( &
        (/0,0,0, 0,0,0, 0,0,0, 1,1,1/), (/3,4/))
!
!  ...for i-th edge of j-th son of a refined triangle
!     T4_ORIENT(i,j) = orientation of the mid-edge node
!                      meaningful for mid-edge nodes that are sons
!                      of the middle node
!
!
      contains
!
!-----------------------------------------------------------------------
!
!  ...In:  middle node father type, node number, son number, ref kind
!  ...Out: element node number of the father node
      function Npar_ref(Type,J,Nson,Nref)
      character(len=4) :: Type
      integer :: J,Nson,Nref, Npar_ref
!
      select case(Type)
      case('mdlt')
        Npar_ref   = PARENT_NODES_T4(1,J,Nson)
      case('mdlq')
        select case(Nref)
        case(11)
          Npar_ref   = PARENT_NODES_Q11(1,J,Nson)
        case(10)
          Npar_ref   = PARENT_NODES_Q10(1,J,Nson)
        case(01)
          Npar_ref   = PARENT_NODES_Q01(1,J,Nson)
        case default
          write(*,7001) Type,J,Nson,Nref; stop 1
 7001     format('Npar_ref: Type,J,Nson,Nref = ',a4,3i4)
        end select
      case default
        write(*,7001) Type,J,Nson,Nref; stop 1
      end select
!
      end function Npar_ref
!
!
      function Nson_ref(Type,J,Nson,Nref)
      character(len=4) :: Type
      integer :: J,Nson,Nref, Nson_ref
!
      select case(Type)
      case('mdlt')
        Nson_ref   = PARENT_NODES_T4(2,J,Nson)
      case('mdlq')
        select case(Nref)
        case(11)
          Nson_ref   = PARENT_NODES_Q11(2,J,Nson)
        case(10)
          Nson_ref   = PARENT_NODES_Q10(2,J,Nson)
        case(01)
          Nson_ref   = PARENT_NODES_Q01(2,J,Nson)
        case default
          write(*,7001) Type,J,Nson,Nref; stop 1
 7001     format('Nson_ref: Type,J,Nson,Nref = ',a4,3i4)
        end select
      case default
        write(*,7001) Type,J,Nson,Nref; stop 1
      end select
!
      end function Nson_ref
!
!
      function Nort_ref(Type,J,Nson,Nref)
      character(len=4) :: Type
      integer :: J,Nson,Nref, Nort_ref
!
      select case(Type)
      case('mdlt')
        select case(J)
        case(1,2,3,7)
          Nort_ref   = 0
        case(4,5,6)
          Nort_ref = T4_ORIENT(J-3,Nson)
        case default
          write(*,7001) Type,J,Nson,Nref; stop 1
        end select
      case('mdlq')
        Nort_ref   = 0
      case default
        write(*,7001) Type,J,Nson,Nref; stop 1
 7001   format('Nort_ref: Type,J,Nson,Nref = ',a4,3i4)
      end select
!
      end function Nort_ref
!-----------------------------------------------------------------------
!
!  ...return number of middle node sons for a specifi! refinement
      subroutine nr_mdle_sons(Type,Kref, Nrsons)
      character(len=4) :: Type
      integer Kref,Nrsons
!
      select case(Type)
      case('mdlt')
        select case(Kref)
        case(1)
          Nrsons=4
        case default
          write(*,7001) Type,Kref; stop 1
 7001     format('nr_mdle_sons: Type,Kref = ',a4,i3)
        end select
      case('mdlq')
        select case(Kref)
        case(11)
          Nrsons=4
        case(10,01)
          Nrsons=2
        case default
          write(*,7001) Type,Kref; stop 1
        end select
      case default
        write(*,*) Type,Kref; stop 1
      end select
!
      end subroutine nr_mdle_sons
!-----------------------------------------------------------------------
!
!  ...return number of sons for different refinements of different nodes
      subroutine nr_sons(Type,Kref, Nrsons)
!
      character(len=4) :: Type
      integer :: Kref,Nrsons
!
      select case(Type)
      case('medg')
        select case(Kref)
        case(1); Nrsons=3
        case default; go to 10
        end select
      case('mdlt')
        select case(Kref)
        case(1); Nrsons=7
        case default; go to 10
        end select
      case('mdlq')
        select case(Kref)
        case(11); Nrsons=9
        case(10,01); Nrsons=3
        case default; go to 10
        end select
      end select
!
!
      return
 10   write(*,7001) Type,Kref
 7001 format('nr_sons: Type,Kref = ',a4,2x,i3)
      stop 1
!
      end subroutine
!
!-----------------------------------------------------------------------
!
!  ...modify son number 'Is' and node orientation 'Nort' for a son of
!     a mid-edge node according to the orientation 'Norient' of
!     the father
      subroutine rotate_edge(Norient,Is,Nort)
      integer Norient,Is,Nort,imod,j,m
      imod(j,m) = j-(j-1)/m*m
      select case(Norient)
      case(0)
      case(1)
        select case(Is)
        case(1,2)
          Is = imod(Is+1,2)
          Nort = mod(Nort+1,2)
        case(3)
        end select
      end select
      end subroutine rotate_edge
!
!
      end module refinements_2D



