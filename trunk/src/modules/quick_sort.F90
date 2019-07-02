module quick_sort

  implicit none
  public  :: qsort_real, qsort_real_index
  private :: partition_real, partition_real_index


  interface qsort
     module procedure qsort_real, qsort_real_index
  end interface

contains

  recursive subroutine qsort_real(A)
    implicit none
    real*8, intent(inout), dimension(:) :: A
    integer :: marker

    if ( size(A) > 1 ) then
       call partition_real( A, marker )
       call qsort_real( A(:marker-1) )
       call qsort_real( A(marker:)   )
    endif

  end subroutine qsort_real

  subroutine partition_real(A, marker)
    implicit none
    real*8,  intent(inout), dimension(:) :: A
    integer, intent(out)                 :: marker
    integer :: i, j
    real :: temp, x

    x = A(1);    i = 0;    j = size(A) + 1

    do
       j = j-1
       do
          if ( A(j).le.x ) exit
          j = j-1
       end do
       i = i + 1
       do
          if ( A(i).ge.x ) exit
          i = i + 1
       end do
       if ( i.lt.j ) then
          ! exchange A(i) and A(j)
          temp = A(i);          A(i) = A(j);          A(j) = temp
       else if ( i.eq.j ) then
          marker = i + 1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine partition_real

  recursive subroutine qsort_real_index(A, N)
    implicit none
    real*8, intent(inout), dimension(:) :: A, N
    integer :: marker

    ! error check
    if ( size(A).ne.size(N) ) then
       write(*,*) 'size of A and N should match each other'
       stop 1
    endif

    if ( size(A) > 1 ) then
       call partition_real_index( A, N, marker )
       call qsort_real_index( A(:marker-1), N(:marker-1) )
       call qsort_real_index( A(marker:),   N(marker:)   )
    endif

  end subroutine qsort_real_index

  subroutine partition_real_index(A, N, marker)
    implicit none
    real*8,  intent(inout), dimension(:) :: A, N
    integer, intent(out)                 :: marker
    integer :: i, j, itemp;
    real :: temp, x

    x = A(1);    i = 0;    j = size(A) + 1
    do
       j = j-1
       do
          if ( A(j).le.x ) exit
          j = j-1
       end do
       i = i+1
       do
          if ( A(i).ge.x ) exit
          i = i+1
       end do

       if ( i.lt.j ) then
          ! exchange A(i) and A(j)
          temp = A(i);          A(i) = A(j);          A(j) = temp
          ! exchange N(i) and N(j)
          itemp = N(i);         N(i) = N(j);          N(j) = itemp

       else if ( i.eq.j ) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    enddo

  end subroutine partition_real_index

end module quick_sort
