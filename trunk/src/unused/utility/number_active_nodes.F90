!> Purpose : Write points with numbering on visit flag
subroutine number_active_nodes(N)
  use data_structure3D
  use element_data
  implicit none
  character(len=4) :: type
  integer, intent(out) :: N
  integer, dimension(27) :: nodesl, norientl
  integer :: mdle, iel, nod, ino, n_nodes
  !
  mdle = 0; N = 0
  do iel=1, NRELES
     call nelcon(mdle, mdle)
     call elem_nodes(mdle, nodesl,norientl)
     !
     type = NODES(mdle)%type
     n_nodes = nvert(type) + nedge(type) + nface(type) + 1
     !
     do ino=1, n_nodes
        nod = nodesl(ino)
        if (NODES(nod)%visit.eq.0) then
           N = N + 1
           NODES(nod)%visit = 1
        end if
     end do
  end do
end subroutine number_active_nodes
