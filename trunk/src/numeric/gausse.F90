!-----------------------------------------------------------------------
!
!   routine name       - gausse
!
!-----------------------------------------------------------------------
!
!   purpose            - solution of a full matrix linear system
!                        of equations using gauss elimination
!                        without pivoting
!
!   arguments
!        in:      gk   - n by n matrix
!                 igk  - row dimension of matrix gk exactly as
!                        specified in the dimension statement in the
!                        calling program
!                 gf   - vector of length n (the right-hand side)
!                 n    - number of equations
!      out:       u    - vector of length n containing the solution
!
!   required routines  - rhsub,tri
!
!-----------------------------------------------------------------------
!> @date Feb 2023
!-----------------------------------------------------------------------
subroutine gausse(gk,igk,gf,u,n)
!
   implicit none
!
   integer :: igk,n
   real(8) :: gk(igk,*),gf(*),u(*)
!
   call tri(gk,igk,n)
   call rhsub(gk,u,gf,igk,n)
!
end subroutine gausse
