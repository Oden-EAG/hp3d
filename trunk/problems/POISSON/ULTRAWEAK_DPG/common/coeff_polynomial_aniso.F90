

subroutine coeff_polynomial_aniso(Nord_mod,Nextract,etype,Bwork,coeff_sum)


    implicit none

    integer,    intent(in)  :: Nord_mod
    integer,    intent(in) :: Nextract(*)
    character(len=4), intent(in)  :: etype
    real(8),    intent(in)  :: Bwork(*)
    real(8),    intent(out) :: coeff_sum(*)
    integer :: nrdofmQ,px,py,pz
    integer :: i,j
    real(8),    allocatable :: data_sort(:)
    integer,    allocatable :: index_sort(:)

    if(etype .eq. 'mdlb') then

        coeff_sum(1:3) = 0.d0

        call ddecode(Nord_mod,px,py,pz)
        nrdofmQ = px * py * pz
        allocate(data_sort(nrdofmQ))
        allocate(index_sort(nrdofmQ))

        do i = 1,nrdofmQ
            data_sort(i) = real(Nextract(i),8)
            index_sort(i) = i
        enddo

        call qsort_duplet_increasing(index_sort,data_sort,nrdofmQ,1,nrdofmQ)

        ! now we can compute coeff_sum

        do j = 2,px

            coeff_sum(1) = coeff_sum(1) + abs(Bwork(index_sort(j)))
    
        enddo
    
        do j = 1,py-1
    
            coeff_sum(2) = coeff_sum(2) + abs(Bwork(index_sort(j*px + 1)))
    
        enddo
    
        do j = 1,pz-1
    
            coeff_sum(3) = coeff_sum(3) + abs(Bwork(index_sort(j * px * py + 1)))
    
        enddo

    endif

end subroutine coeff_polynomial_aniso

recursive subroutine qsort_duplet_increasing(Iel_array,data_sort,N,First,Last)

implicit none 

integer , intent(in)    :: N,First,Last
integer , intent(inout) :: Iel_array(N)
real(8) , intent(inout) :: data_sort(N)


real(8) :: x,pivot
integer :: i,j,l

pivot = data_sort((First+Last) / 2)
i = First
j = Last

do
    !  ...find first element from the left that needs to be swapped
          do while ((data_sort(i) < pivot))
             i = i + 1
          end do
    !  ...find first element from the right that needs to be swapped
          do while ((pivot < data_sort(j)))
             j = j - 1
          end do
    !  ...end loop if no elements need to be swapped
          if (i >= j) exit
    !  ...swap the elements
          l = Iel_array(i); Iel_array(i) = Iel_array(j); Iel_array(j) = l
          x = data_sort(i); data_sort(i) = data_sort(j); data_sort(j) = x
          i = i + 1
          j = j - 1
enddo

if (First < i-1) call qsort_duplet_increasing(Iel_array,data_sort,N,First,i-1 )
if (j+1 < Last)  call qsort_duplet_increasing(Iel_array,data_sort,N,j+1,  Last)

end subroutine qsort_duplet_increasing