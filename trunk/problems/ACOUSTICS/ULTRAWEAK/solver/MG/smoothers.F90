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
   subroutine block_jacobi_smoother(r,n, maxit, z)
!
   use data_structure3D, only: ZERO, ZONE
   use patch_info
   use pcg_info,         only: local_stiff_multiply
   use macro_grid_info,  only: NRELES_COARSE, DLOC
   implicit none
!
   integer, intent(in)     :: n, maxit
#if C_MODE
   complex*16, intent(in)  :: r(n)
   complex*16, intent(out) :: z(n) 
   complex*16              :: r_new(n), r_old(n), zaux(n)
   complex*16, allocatable, save :: zxp(:) 
#else
   real*8,     intent(in)  :: r(n) 
   real*8,     intent(out) :: z(n) 
   real*8                  :: r_new(n), r_old(n), zaux(n)
   real*8, allocatable,  save :: zxp(:) 
#endif 
!$omp threadprivate(zxp)   
!..locals
   integer :: k, i, nrdof, ip,ig, info, kk, iel, ndof
   real*8  :: w, rnorm
   real*8, external :: dznrm2
!
!--------------------------------------------------------   
!..damping factor 
   w = 0.2d0
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
!$omp private(i,nrdof,ip,ig)
!$omp do reduction(+:zaux) schedule(dynamic)
      do i = 1,NRPATCH
!     ...get the size of the block
         nrdof = PTCH(i)%nrdof
!     ...allocate solution vector associated with the patch      
         allocate(zxp(nrdof))
!     ...loop through the patch dof
         do ip = 1,nrdof
!        ...pick up the global dof number
            ig = PTCH(i)%lcon(ip)
            zxp(ip) = r_old(ig)
!        ...end of loop through patch dof
         enddo
!     ...back substitution 
#if C_MODE
         call ZPPTRS('U',nrdof,1,PTCH(i)%zAp,zxp,nrdof,info)
#else
         call DPPTRS('U',nrdof,1,PTCH(i)%zAp,zxp,nrdof,info)
#endif 
!     ...construct the global vector
         do ip = 1,nrdof
!        ...pick up the global dof number
            ig = PTCH(i)%lcon(ip)
!        ...accumulate for the global solution vector        
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
      call local_stiff_multiply(zaux,r_new)
      r_new = r_old - r_new
      r_old = r_new
!
      rnorm = (dznrm2(n,r_old,1))
!
   enddo
!
!
   end subroutine block_jacobi_smoother





