!----------------------------------------------------------------------
!
!   function name      - ij_upper_to_packed
!
!----------------------------------------------------------------------
!> @brief Given index (I,J) of upper triangle (1<=I<=J) of a symmetric
!!        matrix, returns the corresponding index in a packed
!!        columnwise linear array.
!!
!> @param[in] I - row index
!> @param[in] J - col index
!!
!> @date  Feb 2024
!----------------------------------------------------------------------
integer function ij_upper_to_packed(I,J)
!
   implicit none
!
   integer :: I,J
!
#if DEBUG_MODE
   if (I < 1 .or. I > J) then
      write(*,*) 'ij_upper_to_packed: I,J = ',I,J
      stop
   endif
#endif
!
   ij_upper_to_packed = I + (J-1)*J/2
!
end function ij_upper_to_packed
!
!
!----------------------------------------------------------------------
!
!   function name      - ij_lower_to_packed
!
!----------------------------------------------------------------------
!> @brief Given index (I,J) of lower triangle (J<=I<=N) of a symmetric
!!        matrix, returns the corresponding index in a packed
!!        columnwise linear array.
!!
!> @param[in] I - row index
!> @param[in] J - col index
!> @param[in] N - size of matrix
!!
!> @date  Feb 2024
!----------------------------------------------------------------------
integer function ij_lower_to_packed(I,J,N)
!
   implicit none
!
   integer :: I,J,N
!
#if DEBUG_MODE
   if (J > I .or. I > N) then
      write(*,*) 'ij_lower_to_packed: I,J,N = ',I,J,N
      stop
   endif
#endif
!
   ij_lower_to_packed = I + (J-1)*(2*N-J)/2
!
end function ij_lower_to_packed
!
!
!----------------------------------------------------------------------
!
!   function name      - ij_to_packed
!
!----------------------------------------------------------------------
!> @brief Given index (I,J) of lower or upper triangle of a symmetric
!!        matrix, returns the corresponding index in a packed
!!        columnwise linear array.
!!
!> @param[in] I - row index
!> @param[in] J - col index
!> @param[in] N - size of matrix
!> @param[in] N - size of matrix
!!
!> @date  Feb 2024
!----------------------------------------------------------------------
integer function ij_to_packed(I,J,N,UPLO)
!
   implicit none
!
   integer      :: I,J,N
   character(1) :: UPLO
!
   integer, external :: ij_upper_to_packed, ij_lower_to_packed
!
   if (UPLO .eq. 'U') then
      ij_to_packed = ij_upper_to_packed(I,J)
   elseif (UPLO .eq. 'L') then
      ij_to_packed = ij_lower_to_packed(I,J,N)
   else
      write(*,*) 'ij_to_packed: UPLO = ',UPLO
      stop
   endif
!
end function ij_to_packed
