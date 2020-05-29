!> Purpose : show the mesh data structure
subroutine neig_show
  use data_structure3D
  use refinements
  implicit none
  character(len=4)        :: type
  integer                 :: mdle, i
  integer, dimension(4,6) :: neig
  !
  mdle = 0
  do i=1, NRELES
     call nelcon(mdle, mdle)
     call find_neig(mdle, neig)
     type = NODES(mdle)%type
     write(*,*) 'neig_show: mdle, type, neig= ', &
          mdle, type, neig(1,1:nface(type))
  enddo

end subroutine neig_show
