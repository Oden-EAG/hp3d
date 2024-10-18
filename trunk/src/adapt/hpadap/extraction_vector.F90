!
!-----------------------------------------------------------------------
!> @brief  Routine extracts vector for location of smaller symmetric submatrix in larger symmetric matrix
!!         Here larger matrix corresponds to a L2 gram matrix of larger space and smaller submatrix
!!         is the L2 gram matrix of the smaller space
!> @param[in]   Nord_prev       - previous order
!> @param[in]   Nord_current    - current order
!> @param[in]   Nord_glob       - highest possible order
!> @param[in]   NrdofmQ         - number of dofs corresponding to Nord_current
!> @param[in]   Nextract_prev   - previous extraction vector
!> @param[out]  Nextract        - current extraction vector
!> @date May 2024
!-----------------------------------------------------------------------
subroutine extraction_vector_new(Nord_prev,Nord_current,Nord_glob,NrdofmQ,Nextract_prev, Nextract)
!
   implicit none
!
   integer, intent(in) ::  Nord_prev
   integer, intent(in) ::  Nord_current
   integer, intent(in) ::  Nord_glob
   integer, intent(in) ::  NrdofmQ
!
   integer,    intent(in)                        :: Nextract_prev(*)
   integer,    dimension(NrdofmQ),   intent(out) :: Nextract
!
   integer :: nordxp,nordyp,nordzp
   integer :: nordx,nordy,nordz
   integer :: nordxg,nordyg,nordzg
   integer :: i,j,k,m
!
   call ddecode(Nord_current,nordx,nordy,nordz)
   call ddecode(Nord_prev,nordxp,nordyp,nordzp)
   call ddecode(Nord_glob,nordxg,nordyg,nordzg)
!
   m = 0
!..filling the current extraction vector with the previous extraction vector
!  as the previous extraction vector is subset of the current one.
   do i = 1,nordxp*nordyp*nordzp
      m = m + 1
      Nextract(m) = Nextract_prev(m)
   enddo
!..adding new indexes corresponding to increase along px
   do k = 0,nordzp-1
      do j = 0,nordyp-1
         do i = nordxp,nordx-1
               m = m + 1
               Nextract(m) =  k * (nordxg * nordyg) + j * nordxg + i + 1
         enddo
      enddo
   enddo
!..adding new indexes corresponding to increase along py
   do k = 0,nordzp-1
      do j = nordyp,nordy-1
         do i = 0,nordx-1
               m = m + 1
               Nextract(m) =  k * (nordxg * nordyg) + j * nordxg + i + 1
         enddo
      enddo
   enddo
!..adding new indexes corresponding to increase along pz
   do k = nordzp,nordz-1
      do j = 0,nordy-1
         do i = 0,nordx-1
               m = m + 1
               Nextract(m) =  k * (nordxg * nordyg) + j * nordxg + i + 1
         enddo
      enddo
   enddo
end subroutine extraction_vector_new