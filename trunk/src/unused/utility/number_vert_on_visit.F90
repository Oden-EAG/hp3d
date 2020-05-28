!> Purpose : Write points with numbering on visit flag
subroutine number_vert_on_visit(Nr_vert)
  use data_structure3D
  use element_data
  implicit none
  integer, intent(out) :: Nr_vert
  integer, dimension(27) :: nodesl, norientl
  integer :: mdle, iel, nod, ino, nver
  !
  mdle = 0; nver = 0
  do iel=1, NRELES
     call nelcon(mdle, mdle)
     call elem_nodes(mdle, nodesl,norientl)
     !
     do ino=1, nvert(NODES(mdle)%type)
        nod = nodesl(ino)
        if (NODES(nod)%visit.eq.0) then
           nver = nver + 1
           NODES(nod)%visit = nver
        end if
     end do
  end do
  Nr_vert = nver
end subroutine number_vert_on_visit
