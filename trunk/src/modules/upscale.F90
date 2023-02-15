!----------------------------------------------------------------------
!
!   module name        - module upscale
!
!----------------------------------------------------------------------
!
!   latest revision    - Oct 2019
!
!   purpose            - define data point for upscaling of VTK output
!
!----------------------------------------------------------------------
module upscale

  save
!
! define generic type for a visualization object
  type :: vis
!
!    number of vertices required by visualization level (0 - 3)
     integer :: NR_VERT
!
!    list of vertices' coordinates
     real(8), pointer :: VERT(:,:)
!
!    number of elements required by visualization level (0 - 3)
     integer :: NR_ELEM
!
!    list of elements' vertices (indexing starts at 0)
     integer, pointer :: ELEM(:,:)
  endtype vis
!
  type(vis) :: TETR_VIS, PRIS_VIS, HEXA_VIS
!
contains
!
! return cell type (XDMF 2) for vis object
  integer function ivis_type(Type)
    character(len=4) :: Type
    select case(Type)
    case('tetr','mdln'); ivis_type = 6
    case('pris','mdlp'); ivis_type = 8
    case('hexa','mdlb'); ivis_type = 9
    end select
  end function ivis_type
!
! return num points to describe vis object
  integer function nobj_conf(Type)
     character(len=4) :: Type
     select case(Type)
        case('tetr','mdln'); nobj_conf = 4
        case('pris','mdlp'); nobj_conf = 6
        case('hexa','mdlb'); nobj_conf = 8
        case default
           write(*,*) 'nobj_conf: unexpected type. stop.'
           stop
     end select
  end function nobj_conf
!
! return preloaded vis element configuration
  type(vis) function vis_on_type(Type)
    character(len=4) :: Type
    select case(Type)
    case('tetr','mdln'); vis_on_type = TETR_VIS
    case('pris','mdlp'); vis_on_type = PRIS_VIS
    case('hexa','mdlb'); vis_on_type = HEXA_VIS
    end select
  end function vis_on_type
!
! deallocate vis object (configuration)
  subroutine clear_vis(V)
    implicit none
    type(vis), intent(inout) :: V
    V%NR_VERT = 0
    if (associated(V%VERT)) then
       deallocate(V%VERT)
    end if

    V%NR_ELEM = 0
    if (associated(V%ELEM)) then
       deallocate(V%ELEM)
    end if
  end subroutine clear_vis
!
! set up vis element configuration
  subroutine load_vis(V,Fp,Type)
    implicit none
    type(vis)       , intent(inout) :: V
    character(len=*), intent(in)    :: Fp, Type

    integer, parameter :: nin = 29
    integer            :: i, istat, ierr

    ierr = 0
!
!   clear_vis nullifies vis%elem and vis%vert pointers
    call clear_vis(V)

    open(unit=nin, file=Fp,form='formatted', &
         access='sequential',status='unknown')

    read(nin,*) V%NR_VERT
    allocate(V%VERT(0:V%NR_VERT, 3), STAT=istat)
    V%VERT = 0.d0
    ierr = ierr + istat

    do i=0,(V%NR_VERT-1)
       read(nin,*) V%VERT(i,1), V%VERT(i,2), V%VERT(i,3)
    end do

    read(nin,*) V%NR_ELEM
    allocate(V%ELEM(1:V%NR_ELEM, 0:8), STAT=istat)
    V%ELEM = 0
    ierr = ierr + istat

    do i=1,V%NR_ELEM
       select case (Type)
       case('tetr','mdln')
          read(nin,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3), V%ELEM(i,4)
       case('pris','mdlp')
          read(nin,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3),&
               V%ELEM(i,4), V%ELEM(i,5), V%ELEM(i,6)
       case('hexa','mdlb')
          read(nin,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3), V%ELEM(i,4), &
               V%ELEM(i,5), V%ELEM(i,6), V%ELEM(i,7), V%ELEM(i,8)
       end select
    end do

    close(nin)
    ierr = ierr + istat
    if (ierr.ne.0) then
       write(*,*) 'module upscale::Something wrong', ierr
    end if
  end subroutine load_vis
!
! return point coordinates of a point associated with a vis object
  subroutine get_vis_point(V,Idx, Pt)
    implicit none
    type(vis), intent(in) :: V
    integer  , intent(in) :: Idx
    real(8)  , intent(out):: Pt(3)
    Pt(1:3) = V%VERT(Idx,1:3)
  end subroutine get_vis_point
!
! return list of vertices for vis object
  subroutine get_vis_elem(V,Idx,Ioffs, Iverl)
    implicit none
    type(vis), intent(in)  :: V
    integer  , intent(in)  :: Idx, Ioffs
    integer  , intent(out) :: Iverl(8)
    Iverl(1:8) = V%ELEM(Idx,1:8) + Ioffs
  end subroutine get_vis_elem
!
! return number of vis elements of the vis object
  subroutine get_vis_nrelem(Type, Nrelem)
    implicit none
    character(len=4), intent(in)  :: Type
    integer         , intent(out) :: Nrelem
    select case(Type)
       case('tetr','mdln'); Nrelem = TETR_VIS%NR_ELEM
       case('pris','mdlp'); Nrelem = PRIS_VIS%NR_ELEM
       case('hexa','mdlb'); Nrelem = HEXA_VIS%NR_ELEM
       case default
          write(*,*) 'get_vis_nrelem: unexpected type. stop.'
          stop
    end select
  end subroutine get_vis_nrelem
!
!
  subroutine disp_vis(Nout, V, Type)
    implicit none
    type(vis), intent(in):: V
    integer,   intent(in):: Nout
    character(len=*), intent(in):: Type

    integer :: i

    write(Nout,*) Type, ' NR_VERT = ', V%NR_VERT
    do i=0,(V%NR_VERT-1)
       write(Nout,*) V%VERT(i,1:3)
    end do

    write(Nout,*) Type, ' NR_ELEM = ', V%NR_ELEM
    do i=1,V%NR_ELEM
       select case (Type)
       case('tetr','mdln')
          write(Nout,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3), V%ELEM(i,4)
       case('pris','mdlp')
          write(Nout,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3),&
               V%ELEM(i,4), V%ELEM(i,5), V%ELEM(i,6)
       case('hexa','mdlb')
          write(Nout,*) &
               V%ELEM(i,0), &
               V%ELEM(i,1), V%ELEM(i,2), V%ELEM(i,3), V%ELEM(i,4), &
               V%ELEM(i,5), V%ELEM(i,6), V%ELEM(i,7), V%ELEM(i,8)
       end select
    end do
  end subroutine disp_vis
!
!
end module upscale
