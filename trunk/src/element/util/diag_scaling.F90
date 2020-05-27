
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

      integer, intent(in)       :: m, n
      integer                   :: i,j,k
      integer                   :: nk, n1, n2
#if C_MODE
      complex(8), intent(inout) :: GP(m*(m+1)/2)
      complex(8), intent(inout) :: B(m,n)
      complex(8)                :: D(m)
#else
      real(8), intent(inout)    :: GP(m*(m+1)/2)
      real(8), intent(inout)    :: B(m,n)
      real(8)                   :: D(m)
#endif
!  ...function for vector storage of Hermitian Gram matrix
      nk(n1,n2) = (n2-1)*n2/2+n1
!
!-----------------------------------------------------------------------
!
!  ...preconditioning: GP ---> D^-1/2 * GP * D^-1/2
!  ...extract the diagonal of the matrix
      do i=1,m
         k = nk(i,i)
         D(i) = GP(k)
      enddo
!  ...apply the scaling to the Gram matrix
      do i=1,m
         do j=i,m
            k = nk(i,j)
            GP(k) = GP(k)/sqrt(D(i)*D(j))
         enddo
      enddo
!
!  ...preconditioning: B ---> D^-1/2 * B
!  ...apply the scaling to the Stiffness matrix
      do j=1,n
         do i=1,m
            B(i,j) = B(i,j)/sqrt(D(i))
         enddo
      enddo

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
      integer                   :: i,j,k
      integer                   :: nk, n1, n2
#if C_MODE
      complex(8), intent(inout) :: GP(m,m)
      complex(8), intent(inout) :: B(m,n)
      complex(8)                :: D(m)
#else
      real(8), intent(inout)    :: GP(m,m)
      real(8), intent(inout)    :: B(m,n)
      real(8)                   :: D(m)
#endif
!
!-----------------------------------------------------------------------
!
!  ...preconditioning: GP ---> D^-1/2 * GP * D^-1/2
!  ...extract the diagonal of the matrix
      do i=1,m
         D(i) = GP(i,i)
      enddo
!  ...apply the scaling to the Gram matrix
      do j=1,m
         do i=1,j
            GP(i,j) = GP(i,j)/sqrt(D(i)*D(j))
         enddo
      enddo
!
!  ...preconditioning: B ---> D^-1/2 * B
!  ...apply the scaling to the Stiffness matrix
      do j=1,n
         do i=1,m
            B(i,j) = B(i,j)/sqrt(D(i))
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
#if C_MODE
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

