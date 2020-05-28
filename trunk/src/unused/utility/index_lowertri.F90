!---------------------------------------------------------------------------
!> Purpose : Given the row and column indexes of a full square matrix it
!> returns the index pointing to the same coefficient in a vector used
!> to store the lower triangular part of the matrix
!>
!>  In other words,
!>
!>    matrix_lowertri(index_lowertri(row,column))=matrix_full(row,column)
!
!! @param[in]  row          - row index of full matrix
!! @param[in]  column       - column index of full matrix
!! @param[in]  matrix_size  - number of rows=columns of the full square matrix
!---------------------------------------------------------------------------
function index_lowertri(row,column,matrix_size)
  implicit none
  integer :: row,column,matrix_size
  integer :: index_lowertri

  ! row index runs first (lower triangular)
  index_lowertri =  row+((2*matrix_size-column)*(column-1))/2


end function index_lowertri
