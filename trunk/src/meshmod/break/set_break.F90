!---------------------------------------------------------------------------------------------
!> Purpose : determine the infos necessary for node breaking
!!
!! @param[in ] Type_nod  - node type
!! @param[in ] Kref      - refinement kind
!! @param[in ] Kfilter   - refinement filter
!! @param[in ] Nord      - order of approximation for the node
!! @param[in ] Nbc       - BC flags
!! @param[out] Nrsons    - number of sons
!! @param[out] Type_sons - son type
!! @param[out] Norder    - order of approximation for sons
!! @param[out] Nfilter   - refinement filter for sons
!! @param[out] Nbcond    - BC flags for sons
!!
!! rev@Dec 12
!---------------------------------------------------------------------------------------------
subroutine set_break(Type_nod,Kref,Kfilter,Nord,Nbc, Nrsons,Type_sons,Norder,Nfilter,Nbcond)
  implicit none
! 
  character(len=4),                intent(in)  :: Type_nod
  integer,                         intent(in)  :: Kref, Kfilter, Nord, Nbc
!
  integer,          dimension(27), intent(out) :: Norder, Nbcond, Nfilter
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type_sons

! initialize refinement filter for sons
  Nfilter=0
  select case (Type_nod)
! EDGE  
  case('medg') ; call set_edge_break(Kref,Kfilter,Nord, Nrsons,Type_sons,Norder,Nfilter)
! TRIANGLE          
  case('mdlt') ; call set_trig_break(Kref,        Nord, Nrsons,Type_sons,Norder        )
! QUAD          
  case('mdlq') ; call set_quad_break(Kref,        Nord, Nrsons,Type_sons,Norder        )
! BRICK          
  case('mdlb') ; call set_bric_break(Kref,Kfilter,Nord, Nrsons,Type_sons,Norder,Nfilter)
! TET          
  case('mdln') ; call set_tetr_break(Kref,Kfilter,Nord, Nrsons,Type_sons,Norder,Nfilter)
! PRISM          
  case('mdlp') ; call set_pris_break(Kref,Kfilter,Nord, Nrsons,Type_sons,Norder,Nfilter)
! PYRAMID          
  case('mdld') ; call set_pyra_break(Kref,Kfilter,Nord, Nrsons,Type_sons,Norder,Nfilter)
  case default
     write(*,*) 'set_break: NOT SUPPORTED TYPE ', Type_nod
     stop
  endselect
!
! BC flags are inherited directly from the father node
  Nbcond(1:Nrsons)=Nbc
!
!
endsubroutine set_break
