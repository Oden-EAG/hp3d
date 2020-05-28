!
!-----------------------------------------------------------------------
!
!   routine name       - find_elem_type
!
!-----------------------------------------------------------------------
!
!   latest revision:   - June 2018
!
!   purpose:           - routine returns the type of the initial mesh
!                        ancestor for an element
!                        (i.e. 'Linear', 'TraHex', ...)
!
!   arguments:
!           in         - Mdle: middle node of an element
!           out        - Type: type of the initial mesh element (ancestor)
!
!-----------------------------------------------------------------------
!
!
subroutine find_elem_type(Mdle, Type)
!
   use GMP
   use data_structure3D
!
   implicit none
!
   integer         , intent(in)  :: Mdle
   character(len=6), intent(out) :: Type
!
   integer :: nfath, nel,no,nick
!
!-----------------------------------------------------------------------------
!
!..determine initial mesh ancestor
   nfath = NODES(Mdle)%father
   do while(nfath.gt.0)
      nfath = NODES(nfath)%father
   enddo
!
!..determine father type
   nel = abs(nfath)
!
   call decode(ELEMS(nel)%GMPblock, no,nick)
!
   select case(nick)
!
!  ...prism
      case(1)
         Type = PRISMS(no)%type
!
!  ...hexahedron
      case(2)
         Type = HEXAS(no)%type
!
!  ...tetrahedron
      case(3)
         Type = TETRAS(no)%type
!
!  ...pyramid
      case(4)
         Type = PYRAMIDS(no)%type
!
   end select
!
end subroutine find_elem_type
