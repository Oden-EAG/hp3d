module node_types
   implicit none
!  integer, parameter :: NODE_TYPE_EMPTY = -1, &
!                        NODE_TYPE_VERTEX = 0, &
!                        NODE_TYPE_EDGE = 1, &
!                        NODE_TYPE_TRI = 2, &
!                        NODE_TYPE_QUAD = 3, &
!                        NODE_TYPE_TETRA = 4, &
!                        NODE_TYPE_PRISM = 5, &
!                        NODE_TYPE_PYRAMID = 6, &
!                        NODE_TYPE_HEXA = 7
!  ...node types
      integer, parameter :: MDLB = 1
      integer, parameter :: MDLN = 2
      integer, parameter :: MDLP = 3
      integer, parameter :: MDLD = 4
      integer, parameter :: MDLQ = 5
      integer, parameter :: MDLT = 6
      integer, parameter :: MEDG = 7
      integer, parameter :: VERT = 8
!  ...element types
      integer, parameter :: BRIC = 9
      integer, parameter :: TETR = 10
      integer, parameter :: PRIS = 11
      integer, parameter :: PYRA = 12
!  ...face types
      integer, parameter :: TRIA = 13
      integer, parameter :: QUAD = 14
      integer, parameter :: RECT = 15
!  ...edge type
      integer, parameter :: SEGM = 16
!
!
   contains
!
!-----------------------------------------------------------------------
!
!@brief  return string representation of a type
!@date   Feb 2023
      function S_Type(Itype)
         Integer Itype
         character(4) S_Type
         select case(Itype)
            ! node types
            case(MDLB); S_Type = 'mdlb'
            case(MDLN); S_Type = 'mdln'
            case(MDLP); S_Type = 'mdlp'
            case(MDLD); S_Type = 'mdld'
            case(MDLQ); S_Type = 'mdlq'
            case(MDLT); S_Type = 'mdlt'
            case(MEDG); S_Type = 'medg'
            case(VERT); S_Type = 'vert'
            ! element types
            case(BRIC); S_Type = 'bric'
            case(TETR); S_Type = 'tetr'
            case(PRIS); S_Type = 'pris'
            case(PYRA); S_Type = 'pyra'
            ! face types
            case(TRIA); S_Type = 'tria'
            case(QUAD); S_Type = 'quad'
            case(RECT); S_Type = 'rect'
            ! edge type
            case(SEGM); S_Type = 'segm'
            case default
               write(*,*) 'S_Type: Itype = ', Itype
               stop
         end select
      end function S_Type
!
!-----------------------------------------------------------------------
!
!@brief  return integer representation of a type
!@date   Feb 2023
      function I_Type(Stype)
         character(4) Stype
         integer I_Type
         select case(Stype)
            ! node types
            case('mdlb'); I_Type = MDLB
            case('mdln'); I_Type = MDLN
            case('mdlp'); I_Type = MDLP
            case('mdld'); I_Type = MDLD
            case('mdlq'); I_Type = MDLQ
            case('mdlt'); I_Type = MDLT
            case('medg'); I_Type = MEDG
            case('vert'); I_Type = VERT
            ! element types
            case('bric'); I_Type = BRIC
            case('tetr'); I_Type = TETR
            case('pris'); I_Type = PRIS
            case('pyra'); I_Type = PYRA
            ! face types
            case('tria'); I_Type = TRIA
            case('quad'); I_Type = QUAD
            case('rect'); I_Type = RECT
            ! edge type
            case('segm'); I_Type = SEGM
            case default
               write(*,*) 'I_Type: Stype = ', Stype
               stop
         end select
      end function I_Type
!
!-----------------------------------------------------------------------
!
end module node_types
