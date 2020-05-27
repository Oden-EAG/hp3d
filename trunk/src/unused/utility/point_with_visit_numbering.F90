!> Purpose : Write points with numbering on visit flag
subroutine point_with_visit_numbering(Nout, Iflag)
  use data_structure3D
  use element_data
  implicit none
  integer, intent(in) :: Nout, Iflag
  integer :: iel, ino, mdle, nod, nodesl(27), norientl(27)
  !
  mdle = 0
  do iel=1, NRELES
     call nelcon(mdle, mdle)
     call elem_nodes(mdle, nodesl, norientl)
     !
     do ino=1, nvert(NODES(mdle)%type)
        nod = nodesl(ino)
        if (NODES(nod)%visit.gt.0) then 
           if (Iflag.eq.1) then
              write(Nout,7000) NODES(nod)%visit, NODES(ino)%dof%coord(1:NDIMEN,1)
           else
              write(Nout,7001) NODES(nod)%dof%coord(1:NDIMEN,1)
           end if
           NODES(nod)%visit = -NODES(nod)%visit
        end if
     end do
  end do
  !
7000 format(i10, '   ', 3(f20.15,2x))
7001 format('   ', 3(f20.15,2x))
  !
end subroutine point_with_visit_numbering
