!---------------------------------------------------------------------------------
!> Purpose : determines info for breaking a pyramid's middle node
!!
!> @param[in ] Kref    - node refinement kind
!> @param[in ] Nord    - node order of approximation
!> @param[out] Nrsons  - number of sons
!> @param[out] Type    - sons' types
!> @param[out] Norder  - sons' orders of approximation
!!
!> rev@Aug 2019
!---------------------------------------------------------------------------------
subroutine set_pyra_break(Kref,Nord, Nrsons,Type,Norder)
  implicit none
  integer,                         intent(in)  :: Kref,Nord
  integer,          dimension(27), intent(out) :: Norder
  integer,                         intent(out) :: Nrsons
  character(len=4), dimension(27), intent(out) :: Type

  integer :: nordp
!---------------------------------------------------------------------------------
!
  nordp=nord*10+nord
!  
! initialize
  Norder = 0; Type(1:27) = 'none'
!
! select refinement kind
  select case(Kref)
  case(10)
     Nrsons = 7
!                  | INTERIOR NODES            | FACE NODES         | 
     Type(  1:7)=(/ 'mdld','mdlp','mdlp','mdlp','mdlq','mdlq','mdlq' /)
     Norder(1:7)=(/   nord, nordp, nordp, nordp, nordp, nordp, nordp /)
  case default
     Nrsons=0
  end select
!
end subroutine set_pyra_break

