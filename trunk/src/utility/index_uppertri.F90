!---------------------------------------------------------------------------
!> Purpose : Given the row and column indexes of a full square matrix it
!> returns the index pointing to the same coefficient in a vector used
!> to store the upper triangular part of the matrix
!>
!>  In other words, 
!>
!>    matrix_uppertri(index_uppertri(row,column))=matrix_full(row,column)
!
!! @param[in]  row       - row index of full matrix
!! @param[in]  column    - column index of full matrix
!---------------------------------------------------------------------------
function index_uppertri(row,column)
  implicit none
  integer :: row,column
  integer :: index_uppertri

  ! row index runs first (upper triangular)
  index_uppertri = ((column-1)*column)/2+row

end function index_uppertri

