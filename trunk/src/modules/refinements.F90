!----------------------------------------------------------------------
!
!   module name        - module refinements
!
!----------------------------------------------------------------------
!
!   latest revision    - May 10
!
!   purpose            - module defines an information necessary
!                        to reconstruct element to nodes connectivities
!
!
!----------------------------------------------------------------------
!
module refinements
!  save
  !     nod / nson / iref1/ iref 2 ...
  integer, dimension(1:21,1:8,1:3,1:6)    :: TETRA_PAR,TETRA_SON,TETRA_ORT
  integer, dimension(1:21,1:8,0:1,0:1)    :: PRISM_PAR,PRISM_SON,PRISM_ORT
  integer, dimension(1:21,1:4,0:1,0:1)    :: PYRAM_PAR,PYRAM_SON,PYRAM_ORT
  integer, dimension(1:27,1:8,0:1,0:1,0:1):: BRICK_PAR,BRICK_SON,BRICK_ORT
  !

  ! Testing unimplemented refinements crash the code
  ! Even dummy arguments should be avoided
  !integer, dimension(1:13), parameter :: TETRA_REF = &
  !     (/11,12,13, 21,22,23,24,25,26, 31,32,33,34/)
  integer, dimension(1:5), parameter :: TETRA_REF = &
       (/11,12,13, 24,32/)   !  the last two are not operational !
  integer, dimension(1:3),  parameter :: PRISM_REF = (/11,10,01/)
  !integer, dimension(1:2),  parameter :: PYRAM_REF = (/10,01/)
  integer, dimension(1:1),  parameter :: PYRAM_REF = (/10/)  ! inoperational
  integer, dimension(1:7),  parameter :: BRICK_REF = &
       (/111,110,101,011,100,010,001/)
  !
  logical :: ISO_ONLY = .FALSE.
  !
  interface elem_show
     subroutine elem_show_var1(Mdle)
       integer, intent(in)    :: Mdle
     end subroutine elem_show_var1

     subroutine elem_show_var2(Mdle, Type, Nodesl, Norientl)
       character(len=4),       intent(in) :: Type
       integer,                intent(in) :: Mdle
       integer, dimension(27), intent(in) :: nodesl, norientl
     end subroutine elem_show_var2
  end interface
  !
contains
  !
  !-----------------------------------------------------------------------
  subroutine disable_iso_only
    ISO_ONLY = .FALSE.
  end subroutine disable_iso_only
  !-----------------------------------------------------------------------
  subroutine enable_iso_only
    ISO_ONLY = .TRUE.
  end subroutine enable_iso_only
  !-----------------------------------------------------------------------
  subroutine init_refinements_pyhp3d
    character(len=*), parameter :: vis_file_dir = "files/ref"
    call init_refinements(vis_file_dir)
  end subroutine init_refinements_pyhp3d

  subroutine init_refinements(Fp)
    implicit none
    character(len=*), intent(in) :: Fp
    !
    TETRA_PAR = 0;      TETRA_SON = 0;      TETRA_ORT = 0
    PRISM_PAR = 0;      PRISM_SON = 0;      PRISM_ORT = 0
    PYRAM_PAR = 0;      PYRAM_SON = 0;      PYRAM_ORT = 0
    BRICK_PAR = 0;      BRICK_SON = 0;      BRICK_ORT = 0
    !
    call init_ref_2(TETRA_PAR,TETRA_SON,TETRA_ORT,1,1, 0,Fp//'/tetra_11')
    call init_ref_2(TETRA_PAR,TETRA_SON,TETRA_ORT,1,2, 0,Fp//'/tetra_12')
    call init_ref_2(TETRA_PAR,TETRA_SON,TETRA_ORT,1,3, 0,Fp//'/tetra_13')
    !
    call init_ref_2(TETRA_PAR,TETRA_SON,TETRA_ORT,2,4, 0,Fp//'/tetra_24')
    call init_ref_2(TETRA_PAR,TETRA_SON,TETRA_ORT,3,2, 0,Fp//'/tetra_32')
    !
    call init_ref_2(PRISM_PAR,PRISM_SON,PRISM_ORT,1,1, 1,Fp//'/prism_11')
    call init_ref_2(PRISM_PAR,PRISM_SON,PRISM_ORT,1,0, 1,Fp//'/prism_10')
    call init_ref_2(PRISM_PAR,PRISM_SON,PRISM_ORT,0,1, 1,Fp//'/prism_01')
    !
    call init_ref_2(PYRAM_PAR,PYRAM_SON,PYRAM_ORT,1,0, 1,Fp//'/pyram_10')
    !
    call init_ref_3(BRICK_PAR,BRICK_SON,BRICK_ORT,1,1,1, 1,Fp//'/brick_111')
    call init_ref_3(BRICK_PAR,BRICK_SON,BRICK_ORT,1,0,0, 1,Fp//'/brick_100')
    call init_ref_3(BRICK_PAR,BRICK_SON,BRICK_ORT,0,1,0, 1,Fp//'/brick_010')
    call init_ref_3(BRICK_PAR,BRICK_SON,BRICK_ORT,0,0,1, 1,Fp//'/brick_001')
    call init_ref_3(BRICK_PAR,BRICK_SON,BRICK_ORT,1,1,0, 1,Fp//'/brick_110')
    call init_ref_3(BRICK_PAR,BRICK_SON,BRICK_ORT,1,0,1, 1,Fp//'/brick_101')
    call init_ref_3(BRICK_PAR,BRICK_SON,BRICK_ORT,0,1,1, 1,Fp//'/brick_011')
    !
  end subroutine init_refinements
  !
  !-----------------------------------------------------------------------
  subroutine init_ref_2(Par,Son,Ort, Kref1,Kref2, Ibegin, Fp)
    implicit none
    ! input arguments
    integer, dimension(:,:,:,:), intent(out) :: &
         Par, Son, Ort
    integer, intent(in)          :: Kref1, Kref2, Ibegin
    character(len=*), intent(in) :: Fp
    ! local variables
    integer, parameter :: nin = 30
    integer            :: i, nr, j, nson, iref1, iref2, n
    !
    open(unit=nin, file=Fp,form='formatted', &
         access='sequential',status='old',action='read')
    read(nin,*) nr
    do i=1,nr
       read(nin,*) j, nson, iref1, iref2, n
       Par(j, nson, Kref1+Ibegin, Kref2+Ibegin) = n
    end do
    read(nin,*) nr
    do i=1,nr
       read(nin,*) j, nson, iref1, iref2, n
       Son(j, nson, Kref1+Ibegin, Kref2+Ibegin) = n
    end do
    read(nin,*) nr
    do i=1,nr
       read(nin,*) j, nson, iref1, iref2, n
       Ort(j, nson, Kref1+Ibegin, Kref2+Ibegin) = n
    end do
    !
    close(nin)
  end subroutine init_ref_2
  !
  !-----------------------------------------------------------------------
  subroutine init_ref_3(Par,Son,Ort, Kref1,Kref2,Kref3, Ibegin, Fp)
    implicit none
    ! input arguments
    integer, dimension(:,:,:,:,:), intent(inout) :: &
         Par, Son, Ort
    integer, intent(in)          :: Kref1, Kref2, Kref3, Ibegin
    character(len=*), intent(in) :: Fp
    ! local variables
    integer, parameter :: nin = 30
    integer            :: i, nr, j, nson, iref1, iref2, iref3, n
    !
    open(unit=nin, file=Fp,form='formatted', &
         access='sequential',status='old',action='read')
    read(nin,*) nr
    do i=1,nr
       read(nin,*) j, nson, iref1, iref2, iref3, n
       Par(j, nson, Kref1+Ibegin, Kref2+Ibegin, Kref3+Ibegin) = n
    end do
    read(nin,*) nr
    do i=1,nr
       read(nin,*) j, nson, iref1, iref2, iref3, n
       Son(j, nson, Kref1+Ibegin, Kref2+Ibegin, Kref3+Ibegin) = n
    end do
    read(nin,*) nr
    do i=1,nr
       read(nin,*) j, nson, iref1, iref2, iref3, n
       Ort(j, nson, Kref1+Ibegin, Kref2+Ibegin, Kref3+Ibegin) = n
    end do
    !
    close(nin)
  end subroutine init_ref_3
  !
  !-----------------------------------------------------------------------
  subroutine decode_ref(Type,Nref_kind,  Iref1,Iref2,Iref3)
    implicit none
    character(len=4), intent(in) :: Type
    integer, intent(in)  :: Nref_kind
    integer, intent(out) :: Iref1, Iref2, Iref3

    ! local variable
    integer :: iref0
    !
    select case(Type)
    case('mdln','mdlp','mdld')
       call decode(Nref_kind, Iref1,Iref2)
       Iref3 = 0
    case('mdlb')
       call decode(Nref_kind, iref0,Iref3)
       call decode(iref0, Iref1,Iref2)
    end select
    !
  end subroutine decode_ref
  !
  !-----------------------------------------------------------------------
  logical function is_iso_only()
    is_iso_only = ISO_ONLY
  end function is_iso_only
  !-----------------------------------------------------------------------
  integer function kref_kind(I, Type)
    character(len=4) :: Type
    integer :: I
    select case(Type)
     case('mdln'); kref_kind = TETRA_REF(I)
     case('mdlp'); kref_kind = PRISM_REF(I)
     case('mdld'); kref_kind = PYRAM_REF(I)
     case('mdlb'); kref_kind = BRICK_REF(I)
     end select
   end function kref_kind
   !-----------------------------------------------------------------------
  integer function nr_ref(Type)
    character(len=4) :: Type
    ! tetrahedra aniso refinement is not considered to be used in closing
    select case(Type)
    case('mdln'); nr_ref = 3 !nr_ref=13
    case('mdlp'); nr_ref = 3
    case('mdld'); nr_ref = 0
    case('mdlb'); nr_ref = 7
    case default
       nr_ref = 0
    end select
  end function nr_ref
  !-----------------------------------------------------------------------
  integer function npar_ref(Type, J, Nson, Iref1, Iref2, Iref3)
    implicit none
    character(len=4), intent(in) :: Type
    integer, intent(in)          :: &
         J, Nson, Iref1, Iref2, Iref3
    !
    select case(Type)
    case('mdln')
       npar_ref   = TETRA_PAR(J,Nson,Iref1,Iref2)
    case('mdlp')
       npar_ref   = PRISM_PAR(J,Nson,Iref1,Iref2)
    case('mdld')
       npar_ref   = PYRAM_PAR(J,Nson,Iref1,Iref2)
    case('mdlb')
       npar_ref   = BRICK_PAR(J,Nson,Iref1,Iref2,Iref3)
    end select
    !
  end function npar_ref
  !
  !-----------------------------------------------------------------------
  integer function nson_ref(Type, J, Nson, Iref1, Iref2, Iref3)
    implicit none
    character(len=4), intent(in) :: Type
    integer, intent(in)          ::  &
         J, Nson, Iref1, Iref2, Iref3
    !
    select case(Type)
    case('mdln')
       nson_ref   = TETRA_SON(J,Nson,Iref1,Iref2)
    case('mdlp')
       nson_ref   = PRISM_SON(J,Nson,Iref1,Iref2)
    case('mdld')
       nson_ref   = PYRAM_SON(J,Nson,Iref1,Iref2)
    case('mdlb')
       nson_ref   = BRICK_SON(J,Nson,Iref1,Iref2,Iref3)
    end select
    !
  end function nson_ref
  !
  !-----------------------------------------------------------------------
  integer function nort_ref(Type, J, Nson, Iref1, Iref2, Iref3)
    implicit none
    character(len=4), intent(in) :: Type
    integer, intent(in)          :: &
         J, Nson, Iref1, Iref2, Iref3
    !
    select case(Type)
    case('mdln')
       nort_ref   = TETRA_ORT(J,Nson,Iref1,Iref2)
    case('mdlp')
       nort_ref   = PRISM_ORT(J,Nson,Iref1,Iref2)
    case('mdld')
       nort_ref   = PYRAM_ORT(J,Nson,Iref1,Iref2)
    case('mdlb')
       nort_ref   = BRICK_ORT(J,Nson,Iref1,Iref2,Iref3)
    end select
    !
  end function nort_ref
  !
  !-----------------------------------------------------------------------
  !  ...return number of mid-face sons of a mid-face node
  subroutine nr_face_sons(Type,Kref, Nrsons)
    implicit none
    character(len=4), intent(in) :: Type
    integer, intent(in)  ::  Kref
    integer, intent(out) ::  Nrsons
    !
    Nrsons=0
    select case(Type)
    case('mdlt')
       select case(Kref)
       case(1)
          Nrsons=4
       case(2,3,4)
          Nrsons=2
       end select
    case('mdlq')
       select case(Kref)
       case(11)
          Nrsons=4
       case(10, 01)
          Nrsons=2
       end select
    end select
    !
  end subroutine nr_face_sons
  !-----------------------------------------------------------------------
  !
  !  ...return number of middle node sons for a specific refinement
  subroutine nr_mdle_sons(Type,Kref, Nrsons)
    implicit none
    character(len=4), intent(in) :: Type
    integer, intent(in)  ::  Kref
    integer, intent(out) ::  Nrsons
    !
    Nrsons=0
    select case(Type)
    case('mdln')
       select case(Kref)
       case(11,12,13)
          Nrsons=8
       case(21,22,23,24,25,26)
          Nrsons=4
       case(31,32,33,34)
          Nrsons=2
       end select
    case('mdlp')
       select case(Kref)
       case(11)
          Nrsons=8
       case(10)
          Nrsons=4
       case(01)
          Nrsons=2
       end select
    case('mdld')
       select case(Kref)
       case(10,01)
          Nrsons=4
       end select
    case('mdlb')
       select case(Kref)
       case(111)
          Nrsons=8
       case(110,101,011)
          Nrsons=4
       case(100,010,001)
          Nrsons=2
       end select
    end select
    !
  end subroutine nr_mdle_sons
  !-----------------------------------------------------------------------
  !
  !  ...return number of sons for different refinements of different nodes
  subroutine nr_sons(Type,Kref, Nrsons)
    implicit none
    character(len=4), intent(in) :: Type
    integer, intent(in)  ::  Kref
    integer, intent(out) ::  Nrsons
    !
    Nrsons=0
    select case(Type)
    case('medg')
       select case(Kref)
       case(1); Nrsons=3
       case default; go to 10
       end select
    case('mdlt')
       select case(Kref)
       case(1); Nrsons=7
       case(2,3,4); Nrsons=3
       case default; go to 10
       end select
    case('mdlq')
       select case(Kref)
       case(11); Nrsons=9
       case(10,01); Nrsons=3
       case default; go to 10
       end select
    case('mdlp')
       select case(Kref)
       case(11); Nrsons=21
       case(10); Nrsons=7
       case(01); Nrsons=3
       case default; go to 10
       end select
    case('mdln')
       select case(Kref)
       case(11,12,13); Nrsons=17
       case(21,22,23,24,25,26); Nrsons=7
       case(31,32,33,34); Nrsons=3
       case default; go to 10
       end select
    case('mdld')
       select case(Kref)
       case(10,01); Nrsons=7
       case default; go to 10
       end select
    case('mdlb')
       select case(Kref)
       case(111); Nrsons=27
       case(110,101,011); Nrsons=9
       case(100,010,001); Nrsons=9
       case default; go to 10
       end select
    end select
    !
    !
    return
10  write(*,7001) Type,Kref
7001 format('nr_sons: Type,Kref = ',a5,2x,i3)
    stop 1
    !
  end subroutine nr_sons
  !-----------------------------------------------------------------------
  !
  !  ...modify son number 'Is' and node orientation 'Nort' for a son of
  !     a mid-edge node according to the orientation 'Norient' of
  !     the father
  subroutine rotate_edge(Norient,Is,Nort)
    implicit none
    integer, external      :: imod
    integer, intent(in)    :: Norient
    integer, intent(inout) :: Is,Nort

    select case(Norient)
    case(0)
    case(1)
       select case(Is)
       case(1,2)
          Is   = imod(Is+1,2)
          Nort =  mod(Nort+1,2)
       case(3)
       end select
    end select
  end subroutine rotate_edge
  !-----------------------------------------------------------------------
  !
  !  ...modify son number 'Is' and node orientation 'Nort' for a son of
  !     a mid-triangle node according to the orientation 'Norient' of
  !     the father, here Iref is the local and Ireff is the global
  !     refinement flag for the face
  subroutine rotate_trian(Iref,Ireff,Norient, Is,Nort)
    implicit none

    ! Input arguments
    integer, intent(in)    :: Iref,Ireff,Norient
    integer, intent(inout) :: Is,Nort

    ! Local variables
    integer :: nr_rotat,nr_flips,i,nrr,nrf


    integer, parameter, dimension(1:3,0:5) :: nson &
         = reshape( (/1,2,3, 3,1,2, 2,3,1, 1,3,2, 2,1,3, 3,2,1/), &
                    (/3,6/) )
    integer, parameter, dimension(1:3,0:5) :: nsgn &
         = reshape( (/0,0,0, 0,1,1, 1,1,0, 1,0,0, 0,0,1, 1,1,1/), &
                    (/3,6/) )
    integer, parameter, dimension(0:2,0:1) :: nrotat_trig &
         = reshape( (/0,1,2, 0,2,1/), &
                    (/3,2/) )
    integer, parameter, dimension(0:5,2:4) :: norient_quad &
         = reshape( (/0,1,2, 7,5,6,   0,2,3, 4,5,7,   0,1,3, 4,5,6/), &
                    (/6,3/) )
    integer, parameter, dimension(0:5,2:4) :: edge_orient &
         = reshape( (/0,0,1,1,0,1,    0,1,1,0,0,1,    0,1,0,0,1,1/), &
                    (/6,3/) )
    integer, parameter, dimension(0:5,2:4) :: irefg  &
         = reshape( (/2,4,3, 2,3,4,   3,2,4, 4,2,3,   4,3,2, 3,4,2/), &
                    (/6,3/) )

!   iwork(i,j) = orientation of a quad system of coordinates having orientation
!                i (0 <= i <= 7) wrt a system of coordinates having orientation
!                j (0 <= j <= 7)
!
    integer, parameter, dimension(0:7,0:7) :: iwork =      &
           reshape( (/0,1,2,3,4,5,6,7,                     &
                      1,2,3,0,5,6,7,4,                     &
                      2,3,0,1,6,7,4,5,                     &
                      3,0,1,2,7,4,5,6,                     &
                      4,7,6,5,0,3,2,1,                     &
                      5,4,7,6,1,0,3,2,                     &
                      6,5,4,7,2,1,0,3,                     &
                      7,6,5,4,3,2,1,0 /) ,(/8,8/))
!
    select case(Iref)
       !
       !  ...h4 refinement
    case(1)
       !
       !  .....son to father
       nrr = mod(Nort,3)
       nrf = Nort/3
       !
       !  .....father to global
       nr_rotat = mod(Norient,3)
       nr_flips = Norient/3
       !
       !  .....calculation of number of rotations has to be done in the
       !       son reference frame
       nr_rotat = nrotat_trig(nr_rotat, nrf)
       !
       select case(Is)
       case(1,2,3)
          Is = nson(Is,Norient)
          Nort = mod(nr_flips+nrf,2)*3 + mod(nr_rotat+nrr,3)
       case(4)
          Nort = mod(nr_flips+nrf,2)*3 + mod(nr_rotat+nrr,3)
       case(5,6,7)
          i=Is-4; Is = 4+nson(i,Norient)
          Nort = mod(Nort+nsgn(i,Norient),2)
       end select
!
!   refinement into a quad and triangle
    case(2,3,4)
       if (Ireff.ne.irefg(Norient,Iref)) then
          write(*,7001) Norient,Iref,Ireff
7001      format('rotate_trian: Norient,Iref,Ireff = ',3i3)
          Nort = 100
          return
       endif
!
!      select son number [ son number is invariant : 1 - tria ; 2 - quad ; 3 edge ]
       select case(Is)
       case(1) ; Nort=Norient
       case(2) ; Nort=iwork(norient_quad(Norient,Iref),Nort)

         !!! !  .......from parent Triangle
         !!! nr_rotat = mod(norient_quad(Norient, Iref),4)
         !!! nr_flips = norient_quad(Norient, Iref)/4
         !!! !
         !!! !  .......from Quad
         !!! nrr      = mod(Nort,4)
         !!! nrf      = Nort/4
         !!! !
         !!! Nort     = mod(nr_flips+nrf,2)*4 + mod(nr_rotat+nrr,4)

       case(3) ; Nort=edge_orient(Norient,Iref)
       endselect
!
    endselect
!
!
endsubroutine rotate_trian
!
!
!
  !-----------------------------------------------------------------------
  !
  !  ...modify son number 'Is' and node orientation 'Nort' for a son of
  !     a mid-triangle node according to the orientation 'Norient' of
  !     the father
  !
  !  ...in:  Iref - local refinement flag for a quad node
  !          Ireff - actual refinement flag for the quad node
  !          Norient - orientation for the quad node
  !          Is - local son number
  !     out: Is - actual son number
  !          Is1 - actual son number for the second refinement
  !          Nort - orientation of the son node
  subroutine rotate_quad(Iref,Ireff,Norient, Is,Is1,Nort)
    implicit none
    integer, intent(in)    :: Iref,Ireff,Norient
    integer, intent(inout) :: Is
    integer, intent(out)   :: Is1,Nort

    integer, parameter, dimension(1:4,0:7) :: nson &
         = reshape( (/1,2,3,4, 4,1,2,3, 3,4,1,2, 2,3,4,1, &
                      1,4,3,2, 2,1,4,3, 3,2,1,4, 4,3,2,1/), &
                    (/4,8/) )
    integer, parameter, dimension(1:4) :: nsgn1 = (/0,0,1,1/)
    integer, parameter, dimension(1:4) :: nsgn2 = (/0,1,1,0/)

    integer, parameter, dimension(1:2,1:9) :: h11_to_h10 &
         = reshape( (/1,1, 2,1, 2,2, 1,2, 3,1, 2,3, 3,2, 1,3, 3,3/), &
                    (/2,9/))
    integer, parameter, dimension(1:2,1:9) :: h11_to_h01 &
         = reshape( (/1,1, 1,2, 2,2, 2,1, 1,3, 3,2, 2,3, 3,1, 3,3/), &
                    (/2,9/) )

    !  ...for j-th son on h11-refined face which has been h??-refined
    !     H11_TO_H??(1,j) = son number for the father node
    !     H11_TO_H??(2,j) = son number for the node
    !
    integer,parameter, dimension(1:2,1:9,0:7) :: face_orient_h11    &
         = reshape( (/ 1,0, 2,0, 3,0, 4,0, 5,0, 6,0, 7,0, 8,0, 9,0, &
                       4,1, 1,1, 2,1, 3,1, 8,0, 5,1, 6,0, 7,1, 9,0, &
                       3,2, 4,2, 1,2, 2,2, 7,1, 8,1, 5,1, 6,1, 9,0, &
                       2,3, 3,3, 4,3, 1,3, 6,1, 7,0, 8,1, 5,0, 9,0, &
                       1,4, 4,4, 3,4, 2,4, 8,0, 7,0, 6,0, 5,0, 9,0, &
                       2,5, 1,5, 4,5, 3,5, 5,0, 8,1, 7,0, 6,1, 9,0, &
                       3,6, 2,6, 1,6, 4,6, 6,1, 5,1, 8,1, 7,1, 9,0, &
                       4,7, 3,7, 2,7, 1,7, 7,1, 6,0, 5,1, 8,0, 9,0/), &
                     (/2,9,8/) )
    integer, parameter, dimension(1:2,1:3,0:7) :: face_orient_h10 &
         = reshape( (/ 1,0, 2,0, 3,0, &
                       2,1, 1,1, 3,0, &
                       2,2, 1,2, 3,1, &
                       1,3, 2,3, 3,1, &
                       1,4, 2,4, 3,0, &
                       2,5, 1,5, 3,0, &
                       2,6, 1,6, 3,1, &
                       1,7, 2,7, 3,1/), &
                     (/2,3,8/) )
    integer, parameter, dimension(1:2,1:3,0:7) :: face_orient_h01 &
         = reshape( (/ 1,0, 2,0, 3,0, &
                       1,1, 2,1, 3,1, &
                       2,2, 1,2, 3,1, &
                       2,3, 1,3, 3,0, &
                       1,4, 2,4, 3,0, &
                       1,5, 2,5, 3,1, &
                       2,6, 1,6, 3,1, &
                       2,7, 1,7, 3,0/), &
                    (/2,3,8/) )

    !  ...for local son 'Is' of face with orientation 'Norient',
    !     face_orient_??(1,Is,Norient) = actual son no
    !     face_orient_??(2,Is,Norient) = its orientation
    select case(Iref)
    case(11)
       Nort = face_orient_h11(2,Is,Norient)
       Is   = face_orient_h11(1,Is,Norient)
       !
       !  .....'Is' is the son number in face coordinates
       select case(Ireff)
       case(11)
          Is1=0
       case(10)
          Is1 = h11_to_h10(2,Is); Is = h11_to_h10(1,Is)
       case(01)
          Is1 = h11_to_h01(2,Is); Is = h11_to_h01(1,Is)
       case default
          write(*,*) "rotate_quad : ERROR, not a valid Ireff type"
          stop 1
       end select
    case(10)
       Nort = face_orient_h10(2,Is,Norient)
       Is   = face_orient_h10(1,Is,Norient)
       Is1=0
    case(01)
       Nort = face_orient_h01(2,Is,Norient)
       Is   = face_orient_h01(1,Is,Norient)
       Is1=0
    case default
       write(*,*) "rotate_quad : ERROR, not a valid Iref type"
       stop 1
    end select
    !
  end subroutine rotate_quad
  !
end module refinements

