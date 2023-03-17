!-----------------------------------------------------------------------
!
!   routine name       - find_elem_type
!
!-----------------------------------------------------------------------
!> @brief Returns the type of the initial mesh ancestor for an element
!!        (e.g. 'Linear', 'TraHex', ...)
!> params[in]  Mdle : middle node of an element
!> params[out] Type : type of the initial mesh element (ancestor)
!> @date Feb 2023
!-----------------------------------------------------------------------
subroutine find_elem_type(Mdle, Etype)
!
   use GMP
   use data_structure3D
!
   implicit none
!
   integer         , intent(in)  :: Mdle
   character(len=6), intent(out) :: Etype
!
   integer :: nfath,nel,no,nick
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
         Etype = PRISMS(no)%type(1:6)
!
!  ...hexahedron
      case(2)
         Etype = HEXAS(no)%type(1:6)
!
!  ...tetrahedron
      case(3)
         Etype = TETRAS(no)%type(1:6)
!
!  ...pyramid
      case(4)
         Etype = PYRAMIDS(no)%type(1:6)
!
      case default
         write(*,*) 'find_elem_type'; stop
!
   end select
!
end subroutine find_elem_type
