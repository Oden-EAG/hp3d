!> @param[in]  str - String containing number
!> @param[out] int - Integer to be converted
subroutine str2int(str, int)
    implicit none
    character(len=*), intent(in)  :: str
    integer,          intent(out) :: int
    read(str, '(i10)') int
end subroutine str2int
