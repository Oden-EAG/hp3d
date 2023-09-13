! @brief  module defines node type parameters
! @date   Feb 2023
module node_types
!
   implicit none
!
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
      integer, parameter :: QUAD = 13
      integer, parameter :: TRIA = 14
      integer, parameter :: RECT = 15
!  ...edge type
      integer, parameter :: SEGM = 16
!
!  ...number of vertices for node types
      integer, parameter, dimension(16) :: NVERT =    &
       (/   8,   4,   6,   5,   4,   3,   2,   0,     &
!     (/ MDLB,MDLN,MDLP,MDLD,MDLQ,MDLT,MEDG,VERT )
            8,   4,   6,   5,   4,   3,   4,   2 /)
!     (/ BRIC,TETR,PRIS,PYRA,QUAD,TRIA,RECT,SEGM )
!
!  ...number of edges for node types
      integer, parameter, dimension(16) :: NEDGE =    &
       (/  12,   6,   9,   8,   4,   3,   0,   0,     &
!     (/ MDLB,MDLN,MDLP,MDLD,MDLQ,MDLT,MEDG,VERT )
           12,   6,   9,   8,   4,   3,   4,   0 /)
!     (/ BRIC,TETR,PRIS,PYRA,QUAD,TRIA,RECT,SEGM )
!
!  ...number of face for node types
      integer, parameter, dimension(16) :: NFACE =    &
       (/   6,   4,   5,   5,   0,   0,   0,   0,     &
!     (/ MDLB,MDLN,MDLP,MDLD,MDLQ,MDLT,MEDG,VERT )
            6,   4,   5,   5,   0,   0,   0,   0 /)
!     (/ BRIC,TETR,PRIS,PYRA,QUAD,TRIA,RECT,SEGM )
!
!  ...node types for each element type
      integer, parameter, dimension(27,4) :: TYPE_NOD = reshape(  &
      ! MDLB: 8 VERT, 12 MEDG, 6 MDLQ, 1 MDLB (27 nodes)
      (/ VERT,VERT, VERT,VERT, VERT,VERT, VERT,VERT,                       &
         MEDG,MEDG, MEDG,MEDG, MEDG,MEDG, MEDG,MEDG, MEDG,MEDG, MEDG,MEDG, &
         MDLQ,MDLQ, MDLQ,MDLQ, MDLQ,MDLQ,                                  &
         MDLB,                                                             &
      ! MDLN: 4 VERT, 6 MEDG, 4 MDLT, 1 MDLN (15 nodes)
         VERT,VERT, VERT,VERT,                                             &
         MEDG,MEDG, MEDG,MEDG, MEDG,MEDG,                                  &
         MDLT,MDLT, MDLT,MDLT,                                             &
         MDLN,                                                             &
         0,0,0, 0,0,0, 0,0,0, 0,0,0,                                       &
      ! MDLP: 6 VERT, 9 MEDG, 2 MDLT, 3 MDLQ, 1 MDLP (21 nodes)
         VERT,VERT, VERT,VERT, VERT,VERT,                                  &
         MEDG,MEDG, MEDG,MEDG, MEDG,MEDG, MEDG,MEDG, MEDG,                 &
         MDLT,MDLT, MDLQ,MDLQ, MDLQ,                                       &
         MDLP,                                                             &
         0,0,0, 0,0,0,                                                     &
      ! MDLD: 5 VERT, 8 MEDG, 1 MDLQ, 4 MDLT, 1 MDLD (19 nodes)
         VERT,VERT, VERT,VERT, VERT,                                       &
         MEDG,MEDG, MEDG,MEDG, MEDG,MEDG, MEDG,MEDG,                       &
         MDLQ,      MDLT,MDLT, MDLT,MDLT,                                  &
         MDLD,                                                             &
         0,0,0, 0,0,0, 0,0                                                 &
      /),(/27,4/))
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
            case(QUAD); S_Type = 'quad'
            case(TRIA); S_Type = 'tria'
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
            case('quad'); I_Type = QUAD
            case('tria'); I_Type = TRIA
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
