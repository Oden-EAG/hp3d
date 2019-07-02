!---------------------------------------------------------------------------------
!> Purpose : determines info for breaking a pyramid's middle node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Kfilter - node refinement filter
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Type    - sons' types
!> @param[out] Norder  - sons' orders of approximation
!> @param[out] Nfilter - sons' refinement filters
!!
!> rev@Dec 12
!---------------------------------------------------------------------------------
subroutine set_pyra_break(Kref,Kfilter,Nord, Nrsons,Type,Norder,Nfilter)
  implicit none
  integer,                         intent(in)  :: Kref, Kfilter, Nord
  integer,          dimension(27), intent(out) :: Norder, Nfilter
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type

  integer :: nordp
!---------------------------------------------------------------------------------
!
  nordp=nord*10+nord
!  
! initialize
  Type(1:27)='none'
  Norder(1:27)=0 ; Nfilter(1:27)=0
!  
! select refinement kind
  select case(Kref)
  case(10)
     Nrsons = 7
!                  | INTERIOR NODES            | FACE NODES         | 
     Type(  1:7)=(/ 'mdld','mdlp','mdlp','mdlp','mdlq','mdlq','mdlq' /)
     Norder(1:7)=(/   nord, nordp, nordp, nordp, nordp, nordp, nordp /)
!
!    set up refinement filters
     if (Kfilter.ne.0) then
        Nfilter(1:7)=(/10,10,10,10,0,0,0/)
     endif
  case default
     Nrsons=0
  endselect
!
!
endsubroutine set_pyra_break
