#include "implicit_none.h"
!-----------------------------------------------------------------------
!
!    routine name:      - coo2csc
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Jul 2018
!
!    purpose:           - sparse coo format to csc
!
!    arguments:
!           in/out
!                       - IA  : matrix row indices/matrix row pointers
!                       - JA  : matrix column indices
!                       - XA  : matrix nonzero values
!           in
!                       - nz  : size of triplet arrays (nonzeros)
!           out
!                       - nnz : number of non-duplicate nonzeros
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine coo2csc(IA,JA,XA,nz, nnz)
!
   implicit none
!..declare variables
   integer    , intent(in)    :: nz
   integer    , intent(out)   :: nnz
   integer    , intent(inout) :: IA(nz),JA(nz)
   VTYPE      , intent(inout) :: XA(nz)
!
!..sort and remove duplicates
   call assemble_triplet(JA,IA,XA,nz,nnz)
!
!..create column pointers
   call get_pointers(JA(1:nnz), nnz)
   
!
end subroutine coo2csc
!
!
!-----------------------------------------------------------------------
!
!    routine name       - coo2csr
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Jul 2018
!
!    purpose:           - sparse coo format to csr
!
!    arguments:
!           in/out
!                       - IA  : matrix row indices/matrix row pointers
!                       - JA  : matrix column indices
!                       - XA  : matrix nonzero values
!           in
!                       - nz  : size of triplet arrays (nonzeros)
!           out
!                       - nnz : number of non-duplicate nonzeros
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine coo2csr(IA,JA,XA,nz, nnz)
!
   implicit none
!..declare variables
   integer    , intent(in)    :: nz
   integer    , intent(out)   :: nnz
   integer    , intent(inout) :: IA(nz),JA(nz)
   VTYPE      , intent(inout) :: XA(nz)
!
!..sort and remove duplicates
   call assemble_triplet(IA,JA,XA,nz,nnz)
!
!..create CSR row pointers
   call get_pointers(IA(1:nnz), nnz)
!   
end subroutine coo2csr
!
!
!-----------------------------------------------------------------------
!
!    routine name       - get_pointers
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Aug 2018
!
!    purpose:           - convert row vector indices to pointers 
!                         (compatible to CSR format)
!
!    arguments:
!           in/out
!                       - IA  : matrix row indices/matrix row pointers
!           in
!                       - nnz : number of non-duplicate nonzeros
!
!-----------------------------------------------------------------------
subroutine get_pointers(IA,nnz)

   implicit none

   integer, intent(in)    :: nnz
   integer, intent(inout) :: IA(nnz)
   integer :: i,j,p,t

!..create row pointers
   t = 1
   p = 2
   do i=2,nnz
      j = IA(i)
      if (IA(i) .ne. t) then
         IA(p) = i
         p = p + 1
      endif
      t = j
   enddo
   IA(p) = i

end subroutine get_pointers
!
!
!-----------------------------------------------------------------------
!
!    routine name:      - assemble_triplet
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Jul 2018
!
!    purpose:           - sorts an array of triplets (row,col,val) and
!                         subsequently removes (row,col) duplicates by
!                         reduction (sum) over val
!
!    arguments:
!           in/out
!                       - row : 1D integer array (row indices)
!                       - col : 1D integer array (column indices)
!                       - val : 1D real/complex array (nz values)
!           in
!                       - nz  : size of triplet arrays (nonzeros)
!                       - nnz : number of non-duplicate nonzeros
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine assemble_triplet(row,col,val,nz,nnz)
!
   implicit none
!..declare variables
   integer    , intent(inout) :: row(nz),col(nz)
#if C_MODE
   complex(8) , intent(inout) :: val(nz)
#else
   real(8)    , intent(inout) :: val(nz)
#endif
   integer    , intent(in)    :: nz
   integer    , intent(out)   :: nnz
   integer                    :: i,j,k
   integer    , external      :: OMP_GET_NUM_THREADS
!
!..sort triplets indices (sort key: row major, column minor)
!
!$OMP PARALLEL
!$OMP SINGLE
   k = OMP_GET_NUM_THREADS()
   if (k .eq. 1) then
!  ...sequential sort
      call qsort_triplet(row,col,val,nz,1,nz)
   else
!  ...parallel sort
      call qsort_triplet_omp(nz,row,col,val)
   endif
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
!
!..remove duplicates
   j = 1
   do i = 2,nz
      if (row(j) .ne. row(i) .or. col(j) .ne. col(i)) then
         j = j + 1
         row(j) = row(i)
         col(j) = col(i)
         val(j) = val(i)
      else
         val(j) = val(j) + val(i)
      endif
   enddo   
!
   nnz = j
!
end subroutine assemble_triplet
!
!
!-----------------------------------------------------------------------
!
!    routine name:      - qsort_triplet_old
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Jul 2018
!
!    purpose:           - sorts an array of triplets (IA,JA,XA) with
!                         sort key IA, in ascending order
!
!    arguments:
!           in/out
!                       - ia    : 1D integer array (row indices)
!                       - ja    : 1D integer array (column indices)
!                       - xa    : 1D real/complex array (nz values)
!           in
!                       - n     : size of array
!                       - first : first index of current partition
!                       - last  : last  index of current partition
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
recursive subroutine qsort_triplet_old(ia,ja,xa,n,first,last)
!
   implicit none
!..declare variables
   integer     , intent(in)    :: n,first,last
   integer     , intent(inout) :: ia(n), ja(n)
#if C_MODE   
   complex(8)  , intent(inout) :: xa(n)
   complex(8)                  :: x
#else
   real(8)     , intent(inout) :: xa(n)
   real(8)                     :: x
#endif
   integer                     :: i,j,k,l
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
      l = ia(i); ia(i) = ia(j); ia(j) = l
      l = ja(i); ja(i) = ja(j); ja(j) = l
      x = xa(i); xa(i) = xa(j); xa(j) = x
      i=i+1
      j=j-1
   end do
   if (first < i-1) call qsort_triplet_old(ia,ja,xa,n,first,i-1 )
   if (j+1 < last)  call qsort_triplet_old(ia,ja,xa,n,j+1,  last)
!
!  
end subroutine qsort_triplet_old
!
!
!-----------------------------------------------------------------------
!
!    routine name:      - qsort_triplet
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Jul 2018
!
!    purpose:           - sorts an array of triplets (IA,JA,Val) with
!                         IA (major), JA (minor), in ascending order
!
!    arguments:
!           in/out
!                       - ia    : 1D integer array (row indices)
!                       - ja    : 1D integer array (column indices)
!                       - xa    : 1D real/complex array (nz values)
!           in
!                       - n     : size of array
!                       - first : first index of current partition
!                       - last  : last  index of current partition
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
recursive subroutine qsort_triplet(ia,ja,xa,n,first,last)
!
   implicit none
!..declare variables
   integer     , intent(in)    :: n,first,last
   integer     , intent(inout) :: ia(n), ja(n)
#if C_MODE   
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
!
!
!-----------------------------------------------------------------------
!
!    routine name:      - partition_triplet_omp
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Jul 2018
!
!    purpose:           - partitions an array of triples around a pivot
!
!    arguments:         
!           in/out
!                       - ia : 1D integer array (row indices)
!                       - ja : 1D integer array (column indices)
!                       - xa : 1D real/complex array (nz values)
!           in
!                       - ip : pivot row index (major key)
!                       - jp : pivot column index (minor key)
!                       - n  : size of array
!           out
!                       - k  : partition counter (#elements <= pivot)
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
recursive subroutine partition_triplet_omp(ip,jp,n,k,ia,ja,xa)
!
   implicit none
!..declare variables
   integer    , intent(in)     :: ip,jp,n
   integer    , intent(out)    :: k
   integer    , intent(inout)  :: ia(n), ja(n)
#if C_MODE   
   complex(8) , intent(inout)  :: xa(n)
   complex(8)                  :: xaux
#else
   real(8)    , intent(inout)  :: xa(n)
   real(8)                     :: xaux
#endif
   integer                     :: i,ig,iaux,jaux,kl,kr,kto
!..tuning parameter
   ig = 1000000
!..counter
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
!
!
!-----------------------------------------------------------------------
!
!    routine name:      - qsort_triplet_omp
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Jul 2018
!
!    purpose:           - sorts an array of triplets (IA,JA,Val) with
!                         IA (major), JA (minor), in ascending order
!
!    arguments:
!           in/out
!                       - ia : 1D integer array (row indices)
!                       - ja : 1D integer array (column indices)
!                       - xa : 1D real/complex array (nz values)
!           in
!                       - n  : size of arrays
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
recursive subroutine qsort_triplet_omp(n,ia,ja,xa)
!
   implicit none
!..declare variables
   integer    , intent(in)    :: n
   integer    , intent(inout) :: ia(n), ja(n)
#if C_MODE
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
!  ...parallel sorting
!  ...1. select pivot
      call random_number(r)
      i = INT(1+r*n)
      ip = ia(i)
      jp = ja(i)
!  ...2. partition around pivot
!  ...   n_le is #elements <= pivot
      call partition_triplet_omp(ip,jp,n,n_le,ia,ja,xa)
!  ...3. sort left and right partitions separately and in parallel
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
!
!
!-----------------------------------------------------------------------
!
!    routine name:      - assemble_double
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Nov 2016
!
!    purpose:           - 
!
!    arguments:         - 
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine assemble_double(row,val,nz,nnz)
   implicit none

   integer       :: row(nz)
#if C_MODE
   complex(8)    :: Val(nz)
#else
   real(8)       :: Val(nz)
#endif
   integer       :: nz, nnz
   integer       :: first, last, k, i 

   first = 1
   last  = nz 

!..sort first
!
   call qsort_double(row,val,nz,first, last)
! 
!..remove duplicates
   k=1
   do i=2,nz
      if (row(k) .ne. row(i)) then
         k=k+1
         row(k) = row(i)
         val(k) = val(i)
      else
         val(k) = val(k) + val(i)
      endif
   enddo   
!
   nnz = k
!   
end subroutine assemble_double
!
!
!-----------------------------------------------------------------------
!
!    routine name       - qsort_double
!
!-----------------------------------------------------------------------
!
!    latest revision    - NOV 16
!
!    purpose            - sort the couple JA, Val for 
!                         sparse coordinate format
!
!   arguments:
!     in/out
!             ia       - 1D integer array (indices)
!             Val      - 1D real/complex array (nz values)
!     in
!             first    - first index of array
!             last     - last index of array
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
recursive subroutine qsort_double(ia,val,n,first, last)
!
   implicit none
   integer     :: i, j, n, k, l
   integer     :: ia(n) 
#if C_MODE   
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

