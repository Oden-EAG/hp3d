!----------------------------------------------------------------------
!
!   module name        - pcg_info
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2018
!
!   purpose            - holds all required info for PCG solver
!
!----------------------------------------------------------------------
!
   module pcg_info
!
   integer :: CG_ITER, MAX_ITER
   real*8  :: TOL
!
#if C_MODE   
      complex*16, allocatable :: XSOL(:), RHS(:)
#else
      real*8,     allocatable :: XSOL(:), RHS(:) 
#endif

!..local data needed by the solver
   type sarray
#if C_MODE   
      complex*16, allocatable :: z(:),r(:)
#else
      real*8,     allocatable :: z(:),r(:) 
#endif
   endtype sarray

   type(sarray),  allocatable :: LOC(:)




   contains



   subroutine allocate_local_vectors
!      
   use macro_grid_info, only: NRELES_COARSE, DLOC
!
   implicit none
!
   integer :: iel, ndof
!
!---------------------------------------------------------
!
   allocate(LOC(NRELES_COARSE))

   do iel = 1, NRELES_COARSE   
      ndof = DLOC(iel)%ndof
!      
      allocate(LOC(iel)%r(ndof))
   enddo   

   end subroutine allocate_local_vectors



   subroutine deallocate_local_vectors
!      
   use macro_grid_info, only: NRELES_COARSE
!
   implicit none
!
   integer :: iel

   do iel = 1, NRELES_COARSE   
      deallocate(LOC(iel)%r)
   enddo   

   deallocate(LOC)

   end subroutine deallocate_local_vectors
!   
!-----------------------------------------------------------------------
!
!   routine name       - compute_global_residual
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2018
!
!   purpose            - 
!
!   arguments :
!
!                      - 
!-----------------------------------------------------------------------
! 
   subroutine local_stiff_multiply(p,Ap)
!
   use data_structure3D, only: ZERO, ZONE
   use macro_grid_info,  only: NRELES_COARSE, DLOC, NRDOF_MACRO
!
   implicit none
!
#if C_MODE
   complex*16, intent(in)  :: p(NRDOF_MACRO)
   complex*16, intent(out) :: Ap(NRDOF_MACRO)
   complex*16, allocatable, save :: p_loc(:), Ap_loc(:)
#else
   real*8,     intent(in)  :: p(NRDOF_MACRO)
   real*8,     intent(out) :: Ap(NRDOF_MACRO)
   real*8,     allocatable, save :: p_loc(:), Ap_loc(:)
#endif   
!$omp threadprivate(p_loc, Ap_loc)   
!..locals
   integer :: iel, ndof, k, i
! 
!--------------------------------------------------------------------------
!   
!..initialize residual
   Ap = ZERO
!
!$omp parallel default(shared)   &
!$omp private(iel,ndof,k,i)
!$omp do reduction(+:Ap) schedule(dynamic)
!..loop through the elements
   do iel = 1, NRELES_COARSE
!  ...number of dof
      ndof = DLOC(iel)%ndof
      allocate(p_loc(ndof), Ap_loc(ndof))
!  ...compute element solution
      do k = 1, ndof
         i  = DLOC(iel)%lcon(k)
         p_loc(k) = p(i)
      enddo
!      
!  ...compute the local residual vectors
#if C_MODE
      call ZHPMV('U',ndof,ZONE,DLOC(iel)%zstiff,p_loc,1,ZERO,Ap_loc,1)
#else
      call DSPMV('U',ndof,ZONE,DLOC(iel)%zstiff,p_loc,1,ZERO,Ap_loc,1)
#endif   
!      
!  ...assemble global residual vector
      do k = 1, ndof
         i = DLOC(iel)%lcon(k)
         Ap(i) = Ap(i) + Ap_loc(k)
      enddo
!
      deallocate(p_loc,Ap_loc)
!      
!..end of loop through the elements
   enddo
!$omp end do
!$omp end parallel

   end subroutine local_stiff_multiply
!   
!-----------------------------------------------------------------------
!
!   routine name       - compute_global_residual
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2018
!
!   purpose            - 
!
!   arguments :
!
!                      - 
!-----------------------------------------------------------------------
! 
   subroutine compute_global_residual(Xsol,Rsd)
!
   use data_structure3D, only: ZERO, ZONE
   use macro_grid_info,  only: NRELES_COARSE, DLOC, NRDOF_MACRO
!
   implicit none
!
#if C_MODE
   complex*16, intent(in)  :: Xsol(NRDOF_MACRO)
   complex*16, intent(out) :: Rsd(NRDOF_MACRO)
   complex*16, allocatable, save :: z_loc(:), r_loc(:)
#else
   real*8,     intent(in)  :: Xsol(NRDOF_MACRO)
   real*8,     intent(out) :: Rsd(NRDOF_MACRO)
   real*8,     allocatable, save :: z_loc(:), r_loc(:)
#endif   
!$omp threadprivate(z_loc,r_loc)   
!..locals
   integer :: iel, ndof, k, i
! 
!--------------------------------------------------------------------------
!   
!..initialize residual
   Rsd = ZERO
!
!..loop through the elements
!$omp parallel default(shared)   &
!$omp private(iel,ndof,k,i)
!$omp do reduction(+:Rsd) schedule(dynamic)
   do iel = 1, NRELES_COARSE
!  ...number of dof
      ndof = DLOC(iel)%ndof
      allocate(z_loc(ndof),r_loc(ndof))
!  ...compute element solution
      do k = 1, ndof
         i  = DLOC(iel)%lcon(k)
         z_loc(k) = Xsol(i)
      enddo
!      
!  ...compute the local residual vectors
#if C_MODE
      call ZHPMV('U',ndof,ZONE,DLOC(iel)%zstiff,z_loc,1,ZERO,r_loc,1)
#else
      call DSPMV('U',ndof,ZONE,DLOC(iel)%zstiff,z_loc,1,ZERO,r_loc,1)
#endif   
      r_loc = DLOC(iel)%zbload - r_loc
!      
!  ...assemble global residual vector
      do k = 1, ndof
         i = DLOC(iel)%lcon(k)
         Rsd(i) = Rsd(i) + r_loc(k)
      enddo
!
      deallocate(z_loc,r_loc)
!..end of loop through the elements
   enddo
!$omp end do
!$omp end parallel
!
   end subroutine compute_global_residual
!   
!-----------------------------------------------------------------------
!
!   routine name       - compute_global_residual
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2018
!
!   purpose            - 
!
!   arguments :
!
!                      - 
!-----------------------------------------------------------------------
! 
   subroutine compute_local_residuals(Xsol)
!
   use data_structure3D, only: ZERO, ZONE
   use macro_grid_info,  only: NRELES_COARSE, DLOC, NRDOF_MACRO
!
   implicit none
!
#if C_MODE
   complex*16, intent(in)  :: Xsol(NRDOF_MACRO)
   complex*16, allocatable, save :: z_loc(:), r_loc(:)
#else
   real*8,     intent(in)  :: Xsol(NRDOF_MACRO)
   real*8,     allocatable, save :: z_loc(:), r_loc(:)
#endif   
!$omp threadprivate(z_loc,r_loc)   
!..locals
   integer :: iel, ndof, k, i
! 
!--------------------------------------------------------------------------
!   
!..loop through the elements
!$omp parallel default(shared)   &
!$omp private(iel,ndof,k,i)
!$omp do schedule(dynamic)
   do iel = 1, NRELES_COARSE
!  ...number of dof
      ndof = DLOC(iel)%ndof
      allocate(z_loc(ndof))
!  ...compute element solution
      do k = 1, ndof
         i  = DLOC(iel)%lcon(k)
         z_loc(k) = Xsol(i)
      enddo
!      
!  ...compute the local residual vectors
#if C_MODE
      call ZHPMV('U',ndof,ZONE,DLOC(iel)%zstiff,z_loc,1,ZERO,LOC(iel)%r,1)
#else
      call DSPMV('U',ndof,ZONE,DLOC(iel)%zstiff,z_loc,1,ZERO,LOC(iel)%r,1)
#endif   
      LOC(iel)%r = DLOC(iel)%zbload - LOC(iel)%r
!      
      deallocate(z_loc)
!..end of loop through the elements
   enddo
!$omp end do
!$omp end parallel
!
!
   end subroutine compute_local_residuals
!
!
   end module pcg_info
