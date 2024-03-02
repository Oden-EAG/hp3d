!----------------------------------------------------------------------
!
!   module name      - sorts
!
!----------------------------------------------------------------------
!
!   latest revision  - Feb 2023
!
!   purpose          - Collection of sorting routines
!
!----------------------------------------------------------------------
!
module sorts
!
   implicit none
!
   contains
!
!> @brief - permute index such that Array(Index) is ascending
!> @param[in]   n     - size of Array
!> @param[in]   Array - list of integers to sort
!> @param[out]  Index - permutation of sort
!-----------------------------------------------------------------------
   subroutine sortIndex(n,Array,Index)
!
      implicit none
!
      integer, intent(in)  :: n
      integer, intent(in)  :: Array(n)
      integer, intent(out) :: Index(n)
!
      integer, parameter   :: nn=15, nstack=50
      integer              :: k,i,j,indext,jstack,l,r
      integer              :: istack(nstack)
      integer              :: a
!
!----------------------------------------------------------------------
!
      do j = 1,n
         Index(j) = j
      end do
!
      jstack=0
      l=1
      r=n
      do
         if (r-l < nn) then
            do j=l+1,r
               indext=Index(j)
               a=Array(indext)
               do i=j-1,l,-1
                  if (Array(Index(i)) <= a) exit
                  Index(i+1)=Index(i)
               end do
               Index(i+1)=indext
            end do
            if (jstack == 0) return
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
         else
            k=(l+r)/2
            call swap(Index(k),Index(l+1))
            call exchangeIndex(Index(l),Index(r))
            call exchangeIndex(Index(l+1),Index(r))
            call exchangeIndex(Index(l),Index(l+1))
            i=l+1
            j=r
            indext=Index(l+1)
            a=Array(indext)
!
            do
               do
                  i=i+1
                  if (Array(Index(i)) >= a) exit
               end do
               do
                  j=j-1
                  if (Array(Index(j)) <= a) exit
               end do
               if (j < i) exit
               call swap(Index(i),Index(j))
            end do
!
            Index(l+1)=Index(j)
            Index(j)=indext
            jstack=jstack+2
!
            if (jstack > nstack) then
               write(*,*) 'sortIndex: NSTACK too small in indexArray()'
               stop
            end if
!
            if (r-i+1 >= j-l) then
               istack(jstack)=r
               istack(jstack-1)=i
               r=j-1
            else
               istack(jstack)=j-1
               istack(jstack-1)=l
               l=i
            end if
!
         end if
!
      end do


   contains


!> @brief - sorts two elements
!> @param[inout]  i - first index
!> @param[inout]  j - second index
!-----------------------------------------------------------------------
      subroutine exchangeIndex(i,j)
!
         implicit none
!
         integer, intent(inout) :: i,j
         integer                :: swp
!
         if (Array(j) .lt. Array(i)) then
             swp=i
             i=j
             j=swp
         end if
!
      end subroutine exchangeIndex



!> @brief - swaps two indices
!> @param[inout]  a - first index
!> @param[inout]  b - second index
!-----------------------------------------------------------------------
     pure elemental subroutine swap(a,b)
!
         implicit none
!
         integer, intent(inout) :: a,b
         integer :: dum
!
         dum=a
         a=b
         b=dum
!
      end subroutine swap
!
   end subroutine sortIndex




! Note this routine has similar function to sortIndex but also reorders
! real array (as opposed to integers).
!
!> @brief - sorts an array of duplets (iel,residual)
!> @param[inout] Iel_array - Index array for permutation
!> @param[inout] Residuals - array to sort
!> @param[in]    N         - list of integers to sort
!> @param[in]    First     - first index of current partition (needs to be 1 on initial call)
!> @param[in]    Last      - last index of current partition (needs to be N on initial call)
   recursive subroutine qsort_duplet(Iel_array,Residuals,N,First,Last)
!
      implicit none
!
!  ...declare variables
      integer , intent(in)    :: N,First,Last
      integer , intent(inout) :: Iel_array(N)
      real(8) , intent(inout) :: Residuals(N)
!
      real(8) :: x,pivot
      integer :: i,j,l
!
      pivot = Residuals((First+Last) / 2)
      i = First
      j = Last
!  ...iterate through the array to be sorted
      do
!     ...find first element from the left that needs to be swapped
         do while ((Residuals(i) > pivot))
            i = i + 1
         end do
!     ...find first element from the right that needs to be swapped
         do while ((pivot > Residuals(j)))
            j = j - 1
         end do
!     ...end loop if no elements need to be swapped
         if (i >= j) exit
!     ...swap the elements
         l = Iel_array(i); Iel_array(i) = Iel_array(j); Iel_array(j) = l
         x = Residuals(i); Residuals(i) = Residuals(j); Residuals(j) = x
         i = i + 1
         j = j - 1
      end do
      if (First < i-1) call qsort_duplet(Iel_array,Residuals,N,First,i-1 )
      if (j+1 < Last)  call qsort_duplet(Iel_array,Residuals,N,j+1,  Last)
!
   end subroutine qsort_duplet




! Note: qsort_triplet functions used to assemble compressed sparse formats
!
!> @brief - sorts an array of triplets (IA,JA,Val) with
!!           IA (major), JA (minor), in ascending order
!> @param[inout] ia    - integer array (row indices)
!> @param[inout] ja    - integer array (colun indices)
!> @param[inout] xa    - array of values
!> @param[in]    n     - list of integers to sort
!> @param[in]    first - first index of current partition (1 on initial call)
!> @param[in]    last  - last index of current partition (N on initial call)
   recursive subroutine qsort_triplet(ia,ja,xa,n,first,last)
!
      implicit none
!
!..declare variables
      integer     , intent(in)    :: n,first,last
      integer     , intent(inout) :: ia(n), ja(n)
#if HP3D_COMPLEX
      complex(8)  , intent(inout) :: xa(n)
      complex(8)                  :: x
#else
      real(8)     , intent(inout) :: xa(n)
      real(8)                     :: x
#endif
      integer                     :: i,j,ki,kj,l
!
      ki = ia((first+last) / 2)
      kj = ja((first+last) / 2)
      i = first
      j = last
      do
         do while ((ia(i) < ki) .or. ((ia(i) .eq. ki) .and. (ja(i) < kj)))
            i = i + 1
         end do
         do while ((ki < ia(j)) .or. ((ki .eq. ia(j)) .and. (kj < ja(j))))
            j = j - 1
         end do
         if (i >= j) exit
         l = ia(i); ia(i) = ia(j); ia(j) = l
         l = ja(i); ja(i) = ja(j); ja(j) = l
         x = xa(i); xa(i) = xa(j); xa(j) = x
         i = i + 1
         j = j - 1
      end do
      if (first < i-1) call qsort_triplet(ia,ja,xa,n,first,i-1 )
      if (j+1 < last)  call qsort_triplet(ia,ja,xa,n,j+1,  last)
!
   end subroutine qsort_triplet




!> @brief - partitions triplet array around a pivot
!> @param[in]    ip    - first index of current partition
!> @param[in]    jp    - last index of current partition
!> @param[in]    n     - list of integers to sort
!> @param[out]   k     - pivot
!> @param[inout] ia    - integer array (row indices)
!> @param[inout] ja    - integer array (colun indices)
!> @param[inout] xa    - array of values
   recursive subroutine partition_triplet_omp(ip,jp,n,k,ia,ja,xa)
!
      implicit none
!
!  ...declare variables
      integer    , intent(in)     :: ip,jp,n
      integer    , intent(out)    :: k
      integer    , intent(inout)  :: ia(n), ja(n)
#if HP3D_COMPLEX
      complex(8) , intent(inout)  :: xa(n)
      complex(8)                  :: xaux
#else
      real(8)    , intent(inout)  :: xa(n)
      real(8)                     :: xaux
#endif
      integer                     :: i,ig,iaux,jaux,kl,kr,kto
!  ...tuning parameter
      ig = 1000000
!  ...counter
      k = 0
      if (n < ig) then
!  ...sequential partitioning
      do i = 1, n
         if ((ia(i) < ip) .or. ((ia(i) .eq. ip) .and. (ja(i) .le. jp))) then
            k = k + 1
            iaux = ia(k); ia(k) = ia(i); ia(i) = iaux
            jaux = ja(k); ja(k) = ja(i); ja(i) = jaux
            xaux = xa(k); xa(k) = xa(i); xa(i) = xaux
         endif
      enddo
   else
!  ...parallel partitioning
!  ...1. recursive calls: left and right block
!$OMP TASK DEFAULT(SHARED)
      call partition_triplet_omp(ip,jp,  n/2,kl,ia(1:n/2  ),ja(1:n/2  ),xa(1:n/2  ))
!$OMP END TASK
!$OMP TASK DEFAULT(SHARED)
      call partition_triplet_omp(ip,jp,n-n/2,kr,ia(n/2+1:n),ja(n/2+1:n),xa(n/2+1:n))
!$OMP END TASK
!$OMP TASKWAIT
!
!  ...2. set iterations for parallel loop
      if (n/2 - kl > kr) then
         kto = kr       ! left block has more elements to swap
      else
         kto = n/2 - kl ! right block has more elements to swap
      endif
!
!  ...3. swap elements
!$OMP PARALLEL DO
      do i = 1, kto
         iaux = ia(kl+i); ia(kl+i) = ia(n/2+kr-i+1); ia(n/2+kr-i+1) = iaux
         jaux = ja(kl+i); ja(kl+i) = ja(n/2+kr-i+1); ja(n/2+kr-i+1) = jaux
         xaux = xa(kl+i); xa(kl+i) = xa(n/2+kr-i+1); xa(n/2+kr-i+1) = xaux
      enddo
!$OMP END PARALLEL DO
      k = kl + kr
   endif
!
end subroutine partition_triplet_omp





!> @brief - Threaded version of qsort_triplet
recursive subroutine qsort_triplet_omp(n,ia,ja,xa)
!
   implicit none
!
!..declare variables
   integer    , intent(in)    :: n
   integer    , intent(inout) :: ia(n), ja(n)
#if HP3D_COMPLEX
   complex(8) , intent(inout) :: xa(n)
#else
   real(8)    , intent(inout) :: xa(n)
#endif
   integer                    :: i,ig,ip,jp,n_le
   real                       :: r
!..omp tuning parameter
!..note: deadlock if more than 'ig-1' elements are equal
   ig = 1000
!
   if (n < 2) then
!  ...do nothing
   elseif (n < ig) then
!  ...sequential sorting
         call qsort_triplet(ia,ja,xa,n,1,n)
      else
!
!        Parallel sorting
!
!     ...1. select pivot
         call random_number(r)
         i = INT(1+r*n)
         ip = ia(i)
         jp = ja(i)
!
!     ...2. partition around pivot
!           (n_le is #elements <= pivot)
         call partition_triplet_omp(ip,jp,n,n_le,ia,ja,xa)
!
!     ...3. sort left and right partitions separately and in parallel
         if (n_le > 1) then
!$OMP TASK DEFAULT(SHARED)
            call qsort_triplet_omp(  n_le, ia(1:n_le)  , ja(1:n_le)  , xa(1:n_le)  )
!$OMP END TASK
         endif
         if (n_le < n-1) then
!$OMP TASK DEFAULT(SHARED)
            call qsort_triplet_omp(n-n_le, ia(n_le+1:n), ja(n_le+1:n), xa(n_le+1:n))
!$OMP END TASK
         endif
!$OMP TASKWAIT
      endif
!
   end subroutine qsort_triplet_omp




!  Note: This array is similar to qsort_duplet but sorts by index, instead of value
!
!> @brief - sorts duplet (ia,val) by index
   recursive subroutine qsort_double(ia,val,n,first, last)
!
      implicit none
!
      integer     :: i, j, n, k, l
      integer     :: ia(n)
#if HP3D_COMPLEX
      complex(8)  :: val(n), x
#else
      real(8)     :: val(n), x
#endif
      integer :: first, last
!
      k = ia( (first+last) / 2 )
      i = first
      j = last
      do
         do while (ia(i) < k)
            i=i+1
         end do
         do while (k < ia(j))
            j=j-1
         end do
         if (i >= j) exit
         l = ia(i);   ia(i)  = ia(j);   ia(j)  = l
         x = val(i);  val(i) = val(j);  val(j) = x
         i=i+1
         j=j-1
      end do
      if (first < i-1) call qsort_double(ia,val,n, first, i-1)
      if (j+1 < last)  call qsort_double(ia,val,n, j+1,  last)
!
   end subroutine qsort_double


end module sorts
