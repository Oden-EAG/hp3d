!--------------------------------------------------------------------
!> Purpose : routine performs a global h-refinement
!
!> @date July 2019
!--------------------------------------------------------------------
subroutine global_href
!
   use error
   use data_structure3D
   use environment , only : QUIET_MODE
   use mpi_param   , only : RANK,ROOT
!
   implicit none
!
   integer :: mdle_list(NRELES)
   integer :: nr_elements_to_refine,mdle,i,kref,ierr
!
!..collect elements
   nr_elements_to_refine = NRELES
!
   mdle=0
   do i=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(i) = mdle
   enddo
!
!..break the elements
   do i=1,nr_elements_to_refine
      mdle = mdle_list(i)
      if (is_leaf(mdle)) then
         call get_isoref(mdle, kref)
         call break(mdle,kref)
      endif
   enddo
!
   call refresh
!
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,100) ' # of elements broken         = ', nr_elements_to_refine
      write(*,200) ' # of current elements, nodes = ', NRELES, NRNODS
   endif
 100 format(A,I8)
 200 format(A,I8,', ',I9)
!
end subroutine global_href
