!----------------------------------------------------------------------
!
!   module name      - index_map
!
!----------------------------------------------------------------------
!
!   latest revision  - Feb 2023
!
!   purpose          - Maps item to location in list
!
!----------------------------------------------------------------------
!
   module index_map
!
   implicit none
!
      type Indmap
         integer, allocatable :: sorted_list(:)
         integer, allocatable :: index(:)
         integer :: nlist
      end type Indmap


   contains


!-----------------------------------------------------------------------
!> @brief instantiate Indmap, mapping items in list to index in list
!! @param[in]  list - list of items to map
!! @param[in]  n    - number of items in list
!> @date  Feb 2023
!-----------------------------------------------------------------------
      function Index_map_init(list,n)
!
         use sorts
!
         implicit none
!
         integer, intent(in) :: n,list(n)
         integer :: i,j
         type(indmap) :: Index_map_init
!
!-----------------------------------------------------------------------
!
         Index_map_init%nlist = n
         allocate(Index_map_init%sorted_list(n))
         allocate(Index_map_init%index(n))
!
!     ...check if list sorted
         i = 1
         do 
            if (i .ge. n) exit
            if (list(i) .gt. list(i+1)) exit
            i = i + 1
         enddo
!
!     ...list already sorted
         if (i.eq.n) then
            Index_map_init%sorted_list(1:n) = list(1:n)
            do i=1,n
               Index_map_init%index(i) = i
            enddo
         else
            call sortIndex(n,list,Index_map_init%index)
!
!        ...sort list using permutation given by sortIndex
            do i=1,n
               j = Index_map_init%index(i)
               Index_map_init%sorted_list(i) = list(j)
            enddo
         endif
!
#if DEBUG_MODE
!     ...check that list is sorted
         if(n.gt.1) then
           do i=1,n-1
              if (Index_map_init%sorted_list(i) .gt.     &
                  Index_map_init%sorted_list(i+1)) then
                 write(*,*) 'index_map: bad sort', i, n, Index_map_init%sorted_list(i), Index_map_init%sorted_list(i+1)
                 stop 1
              endif
           enddo
         endif
#endif
!
      end function Index_map_init



!-----------------------------------------------------------------------
!> @brief destroy an instance of Indmap
!! @param[in]  this - Indmap instance to destroy
!> @date  Feb 2023
!-----------------------------------------------------------------------
      subroutine Index_map_finalize(this)
!
         type(Indmap) :: this
!
         this%nlist = 0
         if (allocated(this%sorted_list)) then
            deallocate(this%sorted_list, this%index)
         endif
!
      end subroutine Index_map_finalize



!-----------------------------------------------------------------------
!> @brief      gets index of item in sorted array (using binary search)
!! @param[in]  this - Indmap instance to search
!! @param[in]  item - item to search for
!> @date  Feb 2023
!-----------------------------------------------------------------------
      function map_index(this,item)
!
         use mpi_param
!
         implicit none
!
         type(Indmap), intent(in) :: this
         integer, intent(in) :: item
         integer :: map_index, i,j,k
!
!-----------------------------------------------------------------------
!
!     ...if list not initialized, return
         if (this%nlist .eq. 0) then
            map_index = 0
            return
         endif
!
!     ...if item outside of list bounds, return
         if (item .lt. this%sorted_list(1) .or. &
             item .gt. this%sorted_list(this%nlist)) then
            map_index = 0
            return
         endif
!
!     ...binary search for item in list
         i = 1
         j = this%nlist
         do while(i.lt.j)
            k = (i+j)/2
            if (item .eq. this%sorted_list(k)) then
               map_index = this%index(k)
               return
            elseif (item < this%sorted_list(k)) then
               j = k - 1
            else
               i = k + 1
            endif
         enddo
         if (this%sorted_list(i) .eq. item) then
            map_index = this%index(i)
            return
         elseif (this%sorted_list(j) .eq. item) then
            map_index = this%index(j)
            return
         else
            map_index = 0
            return
         endif
!
      end function map_index
!
   end module index_map
