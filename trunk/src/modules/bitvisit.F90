
!----------------------------------------------------------------------
!
!   module name        - bitvisit
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - Stores visitation flag in an array (bitpacked
!                        or not), simplifies reducing visitation
!
!----------------------------------------------------------------------
!
module bitvisit
!
!..number of integers required to store all visitation
   integer :: nword
!
!..flags stored per word (32 if bitpacked; 1 otherwise)
   integer :: bpw
!
!..visitation flags
   integer, allocatable :: visitation(:)
!
   logical :: is_init
!
!..whether to bit pack visitation flags (packed: fast MPI, slow OpenMP)
   logical :: bitvisit_packed = .false


   contains


!--------------------------------------------------------------------
!> @brief initialize visitation array
!> @date  Feb 2023
!--------------------------------------------------------------------
   subroutine bitvisit_init(n)
!
      implicit none
!
      integer, intent(in) :: n
!
!----------------------------------------------------------------------
!
      if (bitvisit_packed) then
         bpw = 32
      else
         bpw = 1
      endif
!
      nword = (n-1)/bpw + 1
      allocate(visitation(nword))
      visitation(:) = 0
!
      is_init = .true.
!
   end subroutine bitvisit_init



!--------------------------------------------------------------------
!> @brief destroy visitation array
!> @date  Feb 2023
!--------------------------------------------------------------------
   subroutine bitvisit_finalize
!
      deallocate(visitation)
!
      is_init = .false.
!
   end subroutine bitvisit_finalize



!--------------------------------------------------------------------
!> @brief reset visitation flags to 0
!> @date  Feb 2023
!--------------------------------------------------------------------
   subroutine bitvisit_reset
!
      visitation(:) = 0
!
   end subroutine bitvisit_reset



!--------------------------------------------------------------------
!> @brief set visitation array entry
!> @date  Feb 2023
!--------------------------------------------------------------------
   subroutine set_bit(i,val)
!
      integer, intent(in) :: i, val
      integer :: iword
!
      if (bitvisit_packed) then
!     ...get word number for entry
         iword = (i - 1)/bpw + 1
!$omp critical
         if (val .eq. 1) then
            visitation(iword) = ibset(visitation(iword),mod(i-1,bpw))
         else
            visitation(iword) = ibclr(visitation(iword),mod(i-1,bpw))
         endif
!$omp end critical
      else
         visitation(i) = val
      endif
!
   end subroutine set_bit



!--------------------------------------------------------------------
!> @brief set visitation array entry to 1
!> @date  Feb 2023
!--------------------------------------------------------------------
   subroutine visit(i)
!
      integer, intent(in) :: i
!
      if (.not. is_init) then
         write(*,*) 'bitvisit: Need to initialize before visiting'
      endif

      call set_bit(i,1)
!
   end subroutine visit



!--------------------------------------------------------------------
!> @brief set visitation array entry to 0
!> @date  Feb 2023
!--------------------------------------------------------------------
   subroutine unvisit(i)
!
      integer, intent(in) :: i
!
      if (.not. is_init) then
         write(*,*) 'bitvisit: Need to initialize before visiting'
      endif
!
      call set_bit(i,0)
!
   end subroutine unvisit



!--------------------------------------------------------------------
!> @brief set visitation array entry to 0
!> @date  Feb 2023
!--------------------------------------------------------------------
   function get_bit(i)
!
      integer, intent(in) :: i
      logical :: get_bit
      integer :: iword
!
      if (bitvisit_packed) then
!     ...get word number for entry
         iword = (i - 1)/bpw + 1
         get_bit = btest(visitation(iword),mod(i-1,32))
      else
         get_bit = visitation(i).ne.0
      endif
!
   end function get_bit



!--------------------------------------------------------------------
!> @brief test whether value has been visited
!> @date  Feb 2023
!--------------------------------------------------------------------
   function visited(i)
!
      integer, intent(in) :: i
      logical :: visited
!
      visited = get_bit(i)
!
   end function visited



!--------------------------------------------------------------------
!> @brief MPI_AllReduce visitation aray
!> @date  Feb 2023
!--------------------------------------------------------------------
   subroutine reduce_visit
!
      use mpi_param,       only: RANK,ROOT
      use MPI,             only: MPI_INTEGER, MPI_BOR, MPI_IN_PLACE, &
                                 MPI_COMM_WORLD
      use par_mesh,        only: DISTRIBUTED
!
      implicit none
!
      integer :: ierr
!
!----------------------------------------------------------------------
!
      if (.not. DISTRIBUTED) return
!
      call MPI_Allreduce(MPI_IN_PLACE,visitation, nword,MPI_INTEGER,MPI_BOR,MPI_COMM_WORLD,ierr)
!
   end subroutine reduce_visit
!
!
end module bitvisit
