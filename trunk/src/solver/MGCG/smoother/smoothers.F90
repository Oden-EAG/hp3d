!-----------------------------------------------------------------------
!
!   routine name       - block_jacobi_smoother
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2018
!
!   purpose            - perform the smoothing step of the preconditioned
!                        conjugate gradient solver using element wise
!                        computations
!
!   arguments :
!
!         in
!                  r    - the residual vector
!                  n    - size of the residual vector
!              maxit    - maximum number of iterations
!        out
!                  z    - the correction vector
!
!-----------------------------------------------------------------------
!
!
#include "implicit_none.h"
!
   subroutine block_jacobi_smoother(Igrid,r,n, maxit, z)
!

   use mg_data_structure, only: GRID
   use data_structure3D,  only: ZERO, ZONE
   use patch_info
   use pcg_info,          only: global_stiff_multiply, THETA
!
   implicit none
!
   integer, intent(in)  :: Igrid,n, maxit
   VTYPE,   intent(in)  :: r(n)
   VTYPE,   intent(out) :: z(n)
   VTYPE                :: r_new(n), r_old(n), zaux(n)
   VTYPE, allocatable   :: zxp(:)
!..locals
   integer :: k, i, nrdof, ip,ig, info, kk, iel, ndof
   real*8  :: w, rnorm
   real*8, external :: dznrm2
!
!--------------------------------------------------------
!..damping factor
   w = THETA
!
   r_old = r ;  z = ZERO
!
!..smoother iterations
   do k = 1, maxit
!
!  ...initialize auxiliary vector
      zaux = ZERO
!
!  ...loop through the number of patches
!$omp parallel default(shared)   &
!$omp private(i,nrdof,ip,ig,zxp)
!$omp do reduction(+:zaux) schedule(dynamic)
      do i = 1,GRID(Igrid)%nrpatch
!     ...get the size of the block
         nrdof = GRID(Igrid)%ptch(i)%nrdof
!     ...allocate solution vector associated with the patch
         allocate(zxp(nrdof))
!     ...loop through the patch dof
         do ip = 1,nrdof
!        ...pick up the global dof number
            ig = GRID(Igrid)%ptch(i)%lcon(ip)
            zxp(ip) = r_old(ig)
!        ...end of loop through patch dof
         enddo
!     ...back substitution
         ! if (dznrm2(nrdof,zxp,1) .gt. 1e-5) then
#if C_MODE
            call ZPPTRS('U',nrdof,1,GRID(Igrid)%ptch(i)%zAp,zxp,nrdof,info)
#else
            call DPPTRS('U',nrdof,1,GRID(Igrid)%ptch(i)%zAp,zxp,nrdof,info)
#endif
         ! endif
!        ...construct the global vector
            do ip = 1,nrdof
!           ...pick up the global dof number
               ig = GRID(Igrid)%ptch(i)%lcon(ip)
!           ...accumulate for the global solution vector
               zaux(ig) = zaux(ig) + zxp(ip)
            enddo
         deallocate(zxp)
!     ...end of loop through the patches
      enddo
!$omp end do
!$omp end parallel
!
      zaux = w * zaux
!
      z = z + zaux
!  ...compute global residual
      call global_stiff_multiply(Igrid,zaux,r_new,n)
      r_new = r_old - r_new
      r_old = r_new
!
      ! rnorm = dznrm2(n,r_old,1)

   enddo
!
!
   end subroutine block_jacobi_smoother





