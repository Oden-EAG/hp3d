#if DEBUG_MODE

!> Purpose : show the mesh data structure
subroutine mesh_show
  use data_structure3D
  use refinements
  implicit none
  integer :: mdle, i
  !
  write(*,*) 'mesh_show: NRELES = ', NRELES
  mdle = 0
  do i=1, NRELES
     call nelcon(mdle, mdle)
     call elem_show(mdle)
  enddo
  !
end subroutine mesh_show

#endif
