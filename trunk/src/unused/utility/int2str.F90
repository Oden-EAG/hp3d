  subroutine int2str(int, str)
    implicit none
    integer,          intent(in)  :: int ! Integer to be converted
    character(len=*), intent(out) :: str ! String containing number
    character(len=11) :: longstr         ! Longest is -2147483647

    write(longstr,'(I11)') int
    longstr = adjustl(longstr)
    if (len_trim(longstr) > len(str)) then
       write(*,7001) 'int2str: WARNING: can''t fit '//trim(longstr)// &
            ' into a ',len(str),'-character variable'
7001   format(A,I3,A)
    endif
    str = longstr
  end subroutine int2str
