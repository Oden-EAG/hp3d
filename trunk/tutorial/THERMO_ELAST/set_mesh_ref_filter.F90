subroutine set_mesh_ref_filter
  use data_structure3D
  implicit none
  integer :: nel, ndom

  do nel=1,NRELIS
     select case(NODES(nel)%type)
     case('mdlb'); NODES(nel)%ref_filter = 101
     end select
  enddo
end subroutine set_mesh_ref_filter
