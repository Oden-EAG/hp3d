subroutine mesh2hp3d(Nout, Nr_vert)

  use data_structure3D
  use element_data
  implicit none

  integer, intent(in)  :: Nout
  integer, intent(out) :: Nr_vert
  integer :: i, nr_cell, mdle, ndom, nodesl(27), norientl(27)

  ! Step 0: Numbering
  call number_vert_on_visit(Nr_vert)

  ! Step 1: Header
  write(Nout,6000) Nr_vert
6000 format(i10, ' # NR POINTS')

  ! Step 2: Points
  call point_with_visit_numbering(Nout, 0)

  ! Step 4: Cell header
  write(Nout,6005) NRELES
6005 format(i10)

  ! Step 5: Visit number start with 1
  call mdle_with_visit_numbering(Nout, 0, 1)

end subroutine mesh2hp3d
