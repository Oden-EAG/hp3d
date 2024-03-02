!----------------------------------------------------------------------
!
!   module name        - module refinements
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - module defines information necessary to
!                        reconstruct element-to-nodes connectivities
!
!----------------------------------------------------------------------
!
module refinements
!
  use node_types
!
  implicit none
!
!          nod / nson / iref1 / iref2 / ...
  integer, dimension(1:21,1:8,1:3,1:6)    :: TETRA_PAR,TETRA_SON,TETRA_ORT
  integer, dimension(1:21,1:8,0:1,0:1)    :: PRISM_PAR,PRISM_SON,PRISM_ORT
  integer, dimension(1:21,1:4,0:1,0:1)    :: PYRAM_PAR,PYRAM_SON,PYRAM_ORT
  integer, dimension(1:27,1:8,0:1,0:1,0:1):: BRICK_PAR,BRICK_SON,BRICK_ORT
  !

  ! Testing unimplemented refinements crash the code
  ! Even dummy arguments should be avoided
  !integer, dimension(1:13), parameter :: TETRA_REF = &
  !     (/11,12,13, 21,22,23,24,25,26, 31,32,33,34/)
  integer, parameter :: TETRA_REF(5) = (/11,12,13, 24,32/)
  integer, parameter :: PRISM_REF(3) = (/11,10,01/)
  !integer, parameter :: PYRAM_REF(2) = (/10,01/)
  integer, parameter :: PYRAM_REF(1) = (/10/)
  integer, parameter :: BRICK_REF(7) = (/111,110,101,011,100,010,001/)
!
  logical :: ISO_ONLY = .false.
!
!
#if HP3D_DEBUG
!-----------------------------------------------------------------------
!< @date Mar 2023
  interface elem_show
     !
     !< @date Mar 2023
      subroutine elem_show_var1(Mdle)
       integer, intent(in) :: Mdle
      end subroutine elem_show_var1
     !
     !< @date Mar 2023
      subroutine elem_show_var2(Mdle, Ntype, Nodesl, Norientl)
       integer, intent(in) :: Mdle, Ntype
       integer, intent(in) :: Nodesl(27), Norientl(27)
      end subroutine elem_show_var2
     !
  end interface elem_show
#endif
!
!
  contains
!
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  subroutine disable_iso_only
    ISO_ONLY = .false.
  end subroutine disable_iso_only
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  subroutine enable_iso_only
    ISO_ONLY = .true.
  end subroutine enable_iso_only
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  subroutine init_refinements(Fp)
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
!< @date Mar 2023
  subroutine init_ref_2(Par,Son,Ort, Kref1,Kref2, Ibegin, Fp)
    integer, dimension(:,:,:,:), intent(out) :: Par, Son, Ort
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
!< @date Mar 2023
  subroutine init_ref_3(Par,Son,Ort, Kref1,Kref2,Kref3, Ibegin, Fp)
    integer, dimension(:,:,:,:,:), intent(inout) :: Par, Son, Ort
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
!< @date Mar 2023
  subroutine decode_ref(Ntype,Nref_kind,  Iref1,Iref2,Iref3)
    integer, intent(in)  :: Ntype, Nref_kind
    integer, intent(out) :: Iref1, Iref2, Iref3
    ! local variable
    integer :: iref0
    !
    select case(Ntype)
    case(MDLB)
       call decode(Nref_kind, iref0,Iref3)
       call decode(iref0, Iref1,Iref2)
    case(MDLN,MDLP,MDLD)
       call decode(Nref_kind, Iref1,Iref2)
       Iref3 = 0
    end select
  end subroutine decode_ref
!
!-----------------------------------------------------------------------
  logical function is_iso_only()
    is_iso_only = ISO_ONLY
  end function is_iso_only
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  integer function kref_kind(I, Ntype)
    integer :: Ntype
    integer :: I
    select case(Ntype)
      case(MDLB); kref_kind = BRICK_REF(I)
      case(MDLN); kref_kind = TETRA_REF(I)
      case(MDLP); kref_kind = PRISM_REF(I)
      case(MDLD); kref_kind = PYRAM_REF(I)
      case default
        write(*,*) 'kref_kind'; stop
    end select
  end function kref_kind
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  integer function nr_ref(Ntype)
    integer :: Ntype
    ! tetrahedra aniso refinement is not considered to be used in closing
    select case(Ntype)
      case(MDLB); nr_ref = 7
      case(MDLN); nr_ref = 3 !nr_ref=13
      case(MDLP); nr_ref = 3
      case(MDLD); nr_ref = 0
      case default; nr_ref = 0
    end select
  end function nr_ref
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  integer function npar_ref(Ntype,J,Nson,Iref1,Iref2,Iref3)
    integer, intent(in) :: Ntype,J,Nson,Iref1,Iref2,Iref3
    !
    select case(Ntype)
      case(MDLB); npar_ref = BRICK_PAR(J,Nson,Iref1,Iref2,Iref3)
      case(MDLN); npar_ref = TETRA_PAR(J,Nson,Iref1,Iref2)
      case(MDLP); npar_ref = PRISM_PAR(J,Nson,Iref1,Iref2)
      case(MDLD); npar_ref = PYRAM_PAR(J,Nson,Iref1,Iref2)
      case default
        write(*,*) 'npar_ref: Ntype = ',Ntype
        stop
    end select
    !
  end function npar_ref
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  subroutine npar_ref_all(Ntype,Nson,Iref1,Iref2,Iref3, Npar_refs)
    integer, intent(in)  :: Ntype,Nson,Iref1,Iref2,Iref3
    integer, intent(out) :: Npar_refs(27)
    !
    Npar_refs(1:27) = 0
    select case(Ntype)
      case(MDLB); Npar_refs(1:27) = BRICK_PAR(1:27,Nson,Iref1,Iref2,Iref3)
      case(MDLN); Npar_refs(1:21) = TETRA_PAR(1:21,Nson,Iref1,Iref2)
      case(MDLP); Npar_refs(1:21) = PRISM_PAR(1:21,Nson,Iref1,Iref2)
      case(MDLD); Npar_refs(1:21) = PYRAM_PAR(1:21,Nson,Iref1,Iref2)
      case default
        write(*,*) 'npar_ref_all: Ntype = ',Ntype
        stop
    end select
    !
  end subroutine npar_ref_all
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  integer function nson_ref(Ntype,J,Nson,Iref1,Iref2,Iref3)
    integer, intent(in) :: Ntype,J,Nson,Iref1,Iref2,Iref3
    !
    select case(Ntype)
      case(MDLB); nson_ref = BRICK_SON(J,Nson,Iref1,Iref2,Iref3)
      case(MDLN); nson_ref = TETRA_SON(J,Nson,Iref1,Iref2)
      case(MDLP); nson_ref = PRISM_SON(J,Nson,Iref1,Iref2)
      case(MDLD); nson_ref = PYRAM_SON(J,Nson,Iref1,Iref2)
      case default
        write(*,*) 'nson_ref: Ntype = ',Ntype
        stop
    end select
    !
  end function nson_ref
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  subroutine nson_ref_all(Ntype,Nson,Iref1,Iref2,Iref3, Nson_refs)
    integer, intent(in)  :: Ntype,Nson,Iref1,Iref2,Iref3
    integer, intent(out) :: Nson_refs(27)
    !
    Nson_refs(1:27) = 0
    select case(Ntype)
      case(MDLB); Nson_refs(1:27) = BRICK_SON(1:27,Nson,Iref1,Iref2,Iref3)
      case(MDLN); Nson_refs(1:21) = TETRA_SON(1:21,Nson,Iref1,Iref2)
      case(MDLP); Nson_refs(1:21) = PRISM_SON(1:21,Nson,Iref1,Iref2)
      case(MDLD); Nson_refs(1:21) = PYRAM_SON(1:21,Nson,Iref1,Iref2)
      case default
        write(*,*) 'nson_ref_all: Ntype = ',Ntype
        stop
    end select
    !
  end subroutine nson_ref_all
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  integer function nort_ref(Ntype,J,Nson,Iref1,Iref2,Iref3)
    integer, intent(in) :: Ntype,J,Nson,Iref1,Iref2,Iref3
    !
    select case(Ntype)
      case(MDLB); nort_ref = BRICK_ORT(J,Nson,Iref1,Iref2,Iref3)
      case(MDLN); nort_ref = TETRA_ORT(J,Nson,Iref1,Iref2)
      case(MDLP); nort_ref = PRISM_ORT(J,Nson,Iref1,Iref2)
      case(MDLD); nort_ref = PYRAM_ORT(J,Nson,Iref1,Iref2)
      case default
        write(*,*) 'nort_ref: Ntype = ',Ntype
        stop
    end select
    !
  end function nort_ref
!
!-----------------------------------------------------------------------
!< @date Mar 2023
  subroutine nort_ref_all(Ntype,Nson,Iref1,Iref2,Iref3, Nort_refs)
    integer, intent(in)  :: Ntype,Nson,Iref1,Iref2,Iref3
    integer, intent(out) :: Nort_refs(27)
    !
    Nort_refs(1:27) = 0
    select case(Ntype)
      case(MDLB); Nort_refs(1:27) = BRICK_ORT(1:27,Nson,Iref1,Iref2,Iref3)
      case(MDLN); Nort_refs(1:21) = TETRA_ORT(1:21,Nson,Iref1,Iref2)
      case(MDLP); Nort_refs(1:21) = PRISM_ORT(1:21,Nson,Iref1,Iref2)
      case(MDLD); Nort_refs(1:21) = PYRAM_ORT(1:21,Nson,Iref1,Iref2)
      case default
        write(*,*) 'nort_ref_all: Ntype = ',Ntype
        stop
    end select
    !
  end subroutine nort_ref_all
!
!-----------------------------------------------------------------------
!> @brief Returns number of mid-face sons of a mid-face node
!> @date  Feb 2023
  subroutine nr_face_sons(Ntype,Kref, Nrsons)
    integer, intent(in)  :: Ntype,Kref
    integer, intent(out) :: Nrsons
    !
    Nrsons=0
    select case(Ntype)
      case(MDLT)
        select case(Kref)
          case(1)    ; Nrsons=4
          case(2,3,4); Nrsons=2
        end select
      case(MDLQ)
        select case(Kref)
          case(11)   ; Nrsons=4
          case(10,01); Nrsons=2
        end select
    end select
    !
  end subroutine nr_face_sons
!
!-----------------------------------------------------------------------
!> @brief Returns number of middle node sons for a specific refinement
!> @date  Feb 2023
  subroutine nr_mdle_sons(Ntype,Kref, Nrsons)
    integer, intent(in)  :: Ntype,Kref
    integer, intent(out) :: Nrsons
    !
    Nrsons=0
    select case(Ntype)
      case(MDLB)
        select case(Kref)
          case(111)        ; Nrsons=8
          case(110,101,011); Nrsons=4
          case(100,010,001); Nrsons=2
        end select
      case(MDLN)
        select case(Kref)
          case(11,12,13)         ; Nrsons=8
          case(21,22,23,24,25,26); Nrsons=4
          case(31,32,33,34)      ; Nrsons=2
        end select
      case(MDLP)
        select case(Kref)
          case(11); Nrsons=8
          case(10); Nrsons=4
          case(01); Nrsons=2
        end select
      case(MDLD)
        select case(Kref)
          case(10,01); Nrsons=4
        end select
    end select
    !
  end subroutine nr_mdle_sons
!
!-----------------------------------------------------------------------
!> @brief Returns number of sons for different refinements of different
!!        nodes
!> @date  Feb 2023
  subroutine nr_sons(Ntype,Kref, Nrsons)
    integer, intent(in)  :: Ntype,Kref
    integer, intent(out) :: Nrsons
    !
    Nrsons=0
    select case(Ntype)
    case(MEDG)
       select case(Kref)
         case(1)     ; Nrsons=3
         case default; goto 10
       end select
    case(MDLT)
       select case(Kref)
         case(1)    ; Nrsons=7
         case(2,3,4); Nrsons=3
         case default; goto 10
       end select
    case(MDLQ)
       select case(Kref)
         case(11)   ; Nrsons=9
         case(10,01); Nrsons=3
         case default; goto 10
       end select
    case(MDLB)
       select case(Kref)
         case(111)        ; Nrsons=27
         case(110,101,011); Nrsons=9
         case(100,010,001); Nrsons=9
         case default; goto 10
       end select
    case(MDLP)
       select case(Kref)
         case(11); Nrsons=21
         case(10); Nrsons=7
         case(01); Nrsons=3
         case default; goto 10
       end select
    case(MDLN)
       select case(Kref)
         case(11,12,13)         ; Nrsons=17
         case(21,22,23,24,25,26); Nrsons=7
         case(31,32,33,34)      ; Nrsons=3
         case default; goto 10
       end select
    case(MDLD)
       select case(Kref)
         case(10,01) ; Nrsons=7
         case default; goto 10
       end select
    end select
    !
    return
10  write(*,7001) S_Type(Ntype),Kref
7001 format('nr_sons: Type,Kref = ',a5,2x,i3)
    stop 1
    !
  end subroutine nr_sons
!
!-----------------------------------------------------------------------
!> @brief Modifies son number 'Is' and node orientation 'Nort' for a son
!!        of a mid-edge node according to the orientation 'Norient' of
!!        the father
!> @date  Feb 2023
  subroutine rotate_edge(Norient,Is,Nort)
    integer, external      :: imod
    integer, intent(in)    :: Norient
    integer, intent(inout) :: Is,Nort
    !
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
!
!-----------------------------------------------------------------------
!> @brief Modifies son number 'Is' and node orientation 'Nort' for a son
!!        of a mid-triangle node according to the orientation 'Norient'
!!        of the father; here, Iref is the local and Ireff is the global
!!        refinement flag for the face
!> @date  Feb 2023
  subroutine rotate_trian(Iref,Ireff,Norient, Is,Nort)
    integer, intent(in)    :: Iref,Ireff,Norient
    integer, intent(inout) :: Is,Nort
!
    integer :: nr_rotat,nr_flips,i,nrr,nrf
!
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
!
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
       end select
!
    end select
!
end subroutine rotate_trian
!
!-----------------------------------------------------------------------
!> @brief Modify son number 'Is' and node orientation 'Nort' for a son
!!        of a mid-triangle node according to the orientation 'Norient'
!!        of the father
!> @param[in]      Iref    - local refinement flag for a quad node
!> @param[in]      Ireff   - actual refinement flag for the quad node
!> @param[in]      Norient - orientation for the quad node
!> @param[in,out]  Is      - local son number -> actual son number
!> @param[out]     Is1     - actual son number for the second refinement
!> @param[out]     Nort    - orientation of the son node
!> @date  Feb 2023
  subroutine rotate_quad(Iref,Ireff,Norient, Is,Is1,Nort)
    integer, intent(in)    :: Iref,Ireff,Norient
    integer, intent(inout) :: Is
    integer, intent(out)   :: Is1,Nort
!
!    integer, parameter, dimension(1:4,0:7) :: nson &
!         = reshape( (/1,2,3,4, 4,1,2,3, 3,4,1,2, 2,3,4,1, &
!                      1,4,3,2, 2,1,4,3, 3,2,1,4, 4,3,2,1/), &
!                    (/4,8/) )
!
    integer, parameter, dimension(1:2,1:9) :: h11_to_h10 &
         = reshape( (/1,1, 2,1, 2,2, 1,2, 3,1, 2,3, 3,2, 1,3, 3,3/), &
                    (/2,9/))
    integer, parameter, dimension(1:2,1:9) :: h11_to_h01 &
         = reshape( (/1,1, 1,2, 2,2, 2,1, 1,3, 3,2, 2,3, 3,1, 3,3/), &
                    (/2,9/) )
!
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
!
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
            Is1 = 0
          case(10)
            Is1 = h11_to_h10(2,Is)
            Is  = h11_to_h10(1,Is)
          case(01)
            Is1 = h11_to_h01(2,Is)
            Is  = h11_to_h01(1,Is)
          case default
            write(*,*) "rotate_quad : ERROR, not a valid Ireff type"
            stop 1
        end select
      case(10)
        Nort = face_orient_h10(2,Is,Norient)
        Is   = face_orient_h10(1,Is,Norient)
        Is1  = 0
      case(01)
        Nort = face_orient_h01(2,Is,Norient)
        Is   = face_orient_h01(1,Is,Norient)
        Is1  = 0
      case default
        write(*,*) "rotate_quad : ERROR, not a valid Iref type"
        stop 1
    end select
    !
  end subroutine rotate_quad
  !
end module refinements
