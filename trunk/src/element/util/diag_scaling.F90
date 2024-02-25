
!--------------------------------------------------------------------
!
!     routine name      - diag_scaling
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          - routine returns performs diagonal scaling
!                         of the form
!                         GP ---> D^-1/2 * GP * D^-1/2
!                         B  ---> D^-1/2 * B
!
!     arguments:
!
!     in:
!          m            - number of test functions
!          n            - number of trial functions plus one
!                         (for the load vector)
!     in/out:
!          GP           - Gram matrix in packed form
!          B            - Stiffness matrix containing also
!                         the load vector in the last column
!
!-----------------------------------------------------------------------
!
   subroutine diag_scaling_p(m,n,GP,B)
!
      implicit none
!
      integer, intent(in)       :: m, n
      integer                   :: i,j,k
!
#if HP3D_COMPLEX
      complex(8), intent(inout) :: GP(m*(m+1)/2)
      complex(8), intent(inout) :: B(m,n)
#else
      real(8), intent(inout)    :: GP(m*(m+1)/2)
      real(8), intent(inout)    :: B(m,n)
#endif
      real(8)                   :: D(m),diag
      integer                   :: ip
!
      integer, external         :: ij_upper_to_packed
!
!-----------------------------------------------------------------------
!
!  ...preconditioning: GP ---> D^-1/2 * GP * D^-1/2
!  ...extract the diagonal of the matrix
      do i=1,m
         k = ij_upper_to_packed(i,i)
         diag = dsqrt(real(GP(k),8))
         ip = floor(log(diag)/log(2.d0))
         D(i) = 0.5d0**ip
      enddo
!  ...apply the scaling to the Gram matrix
      do j=1,m
         do i=1,j
            k = ij_upper_to_packed(i,j)
            GP(k) = GP(k)*D(i)*D(j)
         enddo
      enddo
!
!  ...preconditioning: B ---> D^-1/2 * B
!  ...apply the scaling to the Stiffness matrix
      do j=1,n
         do i=1,m
            B(i,j) = B(i,j)*D(i)
         enddo
      enddo
!
  end subroutine diag_scaling_p
!
!
!--------------------------------------------------------------------
!
!     routine name      - diag_scaling
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          - routine returns performs diagonal scaling
!                         of the form
!                         GP ---> D^-1/2 * GP * D^-1/2
!                         B  ---> D^-1/2 * B
!
!     arguments:
!
!     in:
!          m            - number of test functions
!          n            - number of trial functions plus one
!                         (for the load vector)
!     in/out:
!          GP           - Gram matrix in packed form
!          B            - Stiffness matrix containing also
!                         the load vector in the last column
!
!-----------------------------------------------------------------------
!
   subroutine diag_scaling(m,n,GP,B)
!
      implicit none

      integer, intent(in)       :: m, n
      integer                   :: i,j
#if HP3D_COMPLEX
      complex(8), intent(inout) :: GP(m,m)
      complex(8), intent(inout) :: B(m,n)
#else
      real(8), intent(inout)    :: GP(m,m)
      real(8), intent(inout)    :: B(m,n)
#endif
      real(8)                   :: D(m),diag
      integer                   :: ip
!
!-----------------------------------------------------------------------
!
!  ...preconditioning: GP ---> D^-1/2 * GP * D^-1/2
!  ...extract the diagonal of the matrix
      do i=1,m
         diag = dsqrt(real(GP(i,i),8))
         ip = floor(log(diag)/log(2.d0))
         D(i) = 0.5d0**ip
      enddo
!  ...apply the scaling to the Gram matrix
      do j=1,m
         do i=1,j
            GP(i,j) = GP(i,j)*D(i)*D(j)
         enddo
      enddo
!
!  ...preconditioning: B ---> D^-1/2 * B
!  ...apply the scaling to the Stiffness matrix
      do j=1,n
         do i=1,m
            B(i,j) = B(i,j)*D(i)
         enddo
      enddo
!
  end subroutine diag_scaling
!
!
!--------------------------------------------------------------------
!
!     routine name      - diag_scaling_s
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          - routine returns performs diagonal scaling
!                         of the form
!                         SA   ---> D^-1/2 * SA * D^-1/2
!                         Rhs  ---> D^-1/2 * Rhs
!
!     arguments:
!
!     in:
!          n            - dimension of the matrix
!       Nrhs            - number of right hand sides
!         Nz            - number of non-zeros
!         IA            - row indices
!         JA            - column indices
!         DA            - diagonal of the matrix
!
!     in/out:
!          SA           - Sparse matrix in coordinate format
!          Rhs          - Right hand side
!
!-----------------------------------------------------------------------
!
   subroutine diag_scaling_s(n,Nrhs,Nz,IA,JA,DA, SA,Rhs)
!
      implicit none

      integer, intent(in)       :: n, Nrhs, Nz, IA(Nz), JA(Nz)
      integer                   :: i,j,k
#if HP3D_COMPLEX
      complex(8), intent(in)    :: DA(n)
      complex(8), intent(inout) :: SA(Nz), Rhs(n,Nrhs)
#else
      real(8),    intent(in)    :: DA(n)
      real(8),    intent(inout) :: SA(Nz), Rhs(n,Nrhs)
#endif
!
!-----------------------------------------------------------------------
!
!  ...preconditioning: SA ---> D^-1/2 * SA * D^-1/2
   do k = 1, Nz
!  ...get the row number
      i = IA(k)
!  ...get the column number
      j = JA(k)
!  ...scale by the ith diagonal element
      SA(k) = SA(k)/sqrt(DA(i)*DA(j))
   enddo
   do i = 1, n
      Rhs(i,:) = Rhs(i,:)/sqrt(DA(i))
   enddo
!
!
   end subroutine diag_scaling_s
