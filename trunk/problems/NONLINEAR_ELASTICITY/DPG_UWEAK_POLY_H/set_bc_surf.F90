! set_bc_surf
! Subroutine to set boundary condition index according to boundary surface number
!
! Problem and geometry dependent
! Current problem and geometry:
! Linear Elasticity + Cubic matrix and a solid spherical inclusion inside
!
subroutine set_bc_surf(NSurf,Ibc)
implicit none
integer, intent(in ) :: NSurf
integer, intent(out) :: Ibc

Ibc = 0

select case(Nsurf)
! Set Prescribed displacement b.c. (1) on all 6 faces of the cube
case(1,2,3,4,5)
  Ibc = 1
case(6)
  Ibc = 1
! The remaining surface is an interface but needs no b.c.
end select

end subroutine