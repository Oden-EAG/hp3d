function delta_tensor(i,j)
  implicit none
  integer, intent(in) :: i, j
  integer :: delta_tensor

  if (i == j) then
     delta_tensor = 1
  else
     delta_tensor = 0
  endif
end function delta_tensor
