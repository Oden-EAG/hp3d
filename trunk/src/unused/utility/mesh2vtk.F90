subroutine mesh2vtk(Nout, Nr_vert)

  use data_structure3D
  use element_data
  implicit none

  integer, intent(in)  :: Nout
  integer, intent(out) :: Nr_vert
  integer :: i, nr_cell, mdle, nodesl(27), norientl(27)

  ! Step 0: Numbering
  call number_vert_on_visit(Nr_vert)

  ! Step 1: Header
  write(Nout,6000)
  write(Nout,6001)
  write(Nout,6002)
  write(Nout,6003)
  write(Nout,6004) Nr_vert

6000 format('# vtk DataFile Version 2.0')
6001 format('hp3d export mesh ')
6002 format('ASCII')
6003 format('DATASET UNSTRUCTURED_GRID')
6004 format('POINTS ',i10, '  double')

  ! Step 2: Points
  call point_with_visit_numbering(Nout, 0)

  ! Step 3: Count cells
  mdle = 0; nr_cell = NRELES
  do i=1,NRELES
     call nelcon(mdle, mdle)
     call elem_nodes(mdle, nodesl,norientl)

     nr_cell = nr_cell + nvert(NODES(mdle)%type)
  enddo

  ! Step 4: Cell header
  write(Nout,6005) NRELES, nr_cell
6005 format('CELLS ', i10, i10)

  ! Step 5: Write cell with begin number is 1
  call mdle_with_visit_numbering(Nout, 1, 1)

  ! Step 6: VTK cell type
  write(Nout,6006) NRELES
6006 format('CELL_TYPES ', i10)
  mdle = 0
  do i=1,NRELES
     call nelcon(mdle, mdle)
     write(Nout,*) ivtk_type(NODES(mdle)%type)
  enddo

end subroutine mesh2vtk
