  subroutine str2int(str, int)
    implicit none
    character(len=*), intent(in)  :: str ! String containing number
    integer,          intent(out) :: int ! Integer to be converted

    read(str, '(i10)') int
  end subroutine str2int
