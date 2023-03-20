!> @date Mar 2023
module error
!
   implicit none
!
   integer, parameter :: ERR_INVALID_VALUE     = -10
   integer, parameter :: ERR_ALLOC_FAILURE     = -11
   integer, parameter :: ERR_OUT_OF_RANGE      = -12
!
   integer, parameter :: SUCCESS               =  0
   integer, parameter :: FAILURE               = -1
!
   contains
!
   subroutine logic_error(Ierr, File, Line)
!
      integer     , intent(in) :: Ierr, Line
      character(*), intent(in) :: File
!
      select case(Ierr)
!
         case(FAILURE)
            write(*,*) 'Error :: Fail to process, at ',FILE,',',LINE
!
         case(ERR_INVALID_VALUE)
            write(*,*) 'Error :: Invalid value, at '  ,FILE,',',LINE
!
         case(ERR_ALLOC_FAILURE)
            write(*,*) 'Error :: Alloc failure, at '  ,FILE,',',LINE
!
         case(ERR_OUT_OF_RANGE)
            write(*,*) 'Error :: Out of range, at '   ,FILE,',',LINE
!
         case default
            write(*,*) 'Error :: Unknown Ierr, at '   ,FILE,',',LINE
!
      end select
!
      call pause
      stop 1
!
   end subroutine logic_error

end module error
