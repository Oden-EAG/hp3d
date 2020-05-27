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
#include "implicit_none.h"

   module pcg_info

   use derived_types
!
!   
   integer :: CG_ITER, MAX_ITER, MAX_SM_ITER
   real*8  :: TOL, THETA
!
   VTYPE, allocatable :: XSOL(:), RHS(:)
!

!..corrections for each level
   type correction
      VTYPE, allocatable     :: y1(:), y2(:), y3(:)
   endtype correction
   type(correction), allocatable  :: G(:)  

   contains
!   
!-----------------------------------------------------------------------
!
!   routine name       - global_stiff_multiply
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
   subroutine global_stiff_multiply(Igrid,p,Ap,n)
!
   use mg_data_structure, only: GRID
   use data_structure3D,  only: ZERO, ZONE
!
   implicit none
!
   integer, intent(in)  :: Igrid,n
   VTYPE,   intent(in)  :: p(n)
   VTYPE,   intent(out) :: Ap(n)
   VTYPE, allocatable   :: p_loc(:), Ap_loc(:)
! $omp threadprivate(p_loc, Ap_loc)   
!..locals
   integer :: iel, ndof, k, i
! 
!--------------------------------------------------------------------------
!   
!..initialize residual
   Ap = ZERO
!
!$omp parallel default(shared)   &
!$omp private(iel,ndof,k,i,p_loc,Ap_loc)
!$omp do reduction(+:Ap) schedule(dynamic)
!..loop through the elements
   do iel = 1, GRID(Igrid)%nreles
!  ...number of dof
      ndof = GRID(Igrid)%dloc(iel)%ndof
      allocate(p_loc(ndof), Ap_loc(ndof))
!  ...compute element solution
      do k = 1, ndof
         i  = GRID(Igrid)%dloc(iel)%lcon(k)
         p_loc(k) = p(i)
      enddo
!      
!  ...compute the local residual vectors
#if C_MODE
      call ZHPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,p_loc,1,ZERO,Ap_loc,1)
#else
      call DSPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,p_loc,1,ZERO,Ap_loc,1)
#endif   
!      
!  ...assemble global residual vector
      do k = 1, ndof
         i = GRID(Igrid)%dloc(iel)%lcon(k)
         Ap(i) = Ap(i) + Ap_loc(k)
      enddo
!
      deallocate(p_loc,Ap_loc)
!      
!..end of loop through the elements
   enddo
!$omp end do
!$omp end parallel

   end subroutine global_stiff_multiply
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
   subroutine compute_global_residual(Igrid,Xsol,Rsd,n)
!
   use data_structure3D,  only: ZERO, ZONE
   use mg_data_structure, only: GRID
!
   implicit none
!
   integer, intent(in)    :: Igrid,n
   VTYPE,   intent(in)    :: Xsol(n)
   VTYPE,   intent(out)   :: Rsd(n)
   VTYPE,   allocatable   :: z_loc(:), r_loc(:)
! $omp threadprivate(z_loc,r_loc)   
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
!$omp private(iel,ndof,k,i,z_loc,r_loc)
!$omp do reduction(+:Rsd) schedule(dynamic)
   do iel = 1, GRID(Igrid)%nreles
!  ...number of dof
      ndof = GRID(Igrid)%dloc(iel)%ndof
      allocate(z_loc(ndof),r_loc(ndof))
!  ...compute element solution
      do k = 1, ndof
         i  = GRID(Igrid)%dloc(iel)%lcon(k)
         z_loc(k) = Xsol(i)
      enddo
!      
!  ...compute the local residual vectors
#if C_MODE
      call ZHPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,z_loc,1,ZERO,r_loc,1)
#else
      call DSPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,z_loc,1,ZERO,r_loc,1)
#endif   
      GRID(Igrid)%loc(iel)%r = GRID(Igrid)%dloc(iel)%zbload - r_loc
!      
!  ...assemble global residual vector
      do k = 1, ndof
         i = GRID(Igrid)%dloc(iel)%lcon(k)
         ! Rsd(i) = Rsd(i) + r_loc(k)
         Rsd(i) = Rsd(i) + GRID(Igrid)%loc(iel)%r(k)
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
   subroutine update_global_residual(Igrid,Rsd,n)
!
   use data_structure3D,  only: ZERO, ZONE
   use mg_data_structure, only: GRID
!
   implicit none
!
   integer, intent(in)    :: Igrid,n
   VTYPE,   intent(out)   :: Rsd(n)
   VTYPE,   allocatable   :: z_loc(:), r_loc(:)
! $omp threadprivate(z_loc,r_loc)   
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
!$omp private(iel,ndof,k,i,z_loc,r_loc)
!$omp do reduction(+:Rsd) schedule(dynamic)
   do iel = 1, GRID(Igrid)%nreles
!  ...number of dof
      ndof = GRID(Igrid)%dloc(iel)%ndof
      allocate(z_loc(ndof),r_loc(ndof)) ; z_loc = zero
!  ...compute global solution
      do k = 1, ndof
         i  = GRID(Igrid)%dloc(iel)%lcon(k)
         G(Igrid)%y2(i) = GRID(Igrid)%loc(iel)%z(k)
      enddo
!      
!  ...compute the local residual vectors
#if C_MODE
      call ZHPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,GRID(Igrid)%loc(iel)%z,1,ZERO,r_loc,1)
#else
      call DSPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,GRID(Igrid)%loc(iel)%z,1,ZERO,r_loc,1)
#endif   
      GRID(Igrid)%loc(iel)%r = GRID(Igrid)%loc(iel)%r - r_loc
!      
!  ...assemble global residual vector
      do k = 1, ndof
         i = GRID(Igrid)%dloc(iel)%lcon(k)
         Rsd(i) = Rsd(i) + GRID(Igrid)%loc(iel)%r(k)
      enddo
!
      deallocate(z_loc,r_loc)
!..end of loop through the elements
   enddo
!$omp end do
!$omp end parallel
!
   end subroutine update_global_residual




!   
!-----------------------------------------------------------------------
!
!   routine name       - compute_local_residuals
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
   subroutine compute_local_residuals(Igrid,Xsol,n)
!
   use mg_data_structure, only: GRID
   use data_structure3D, only: ZERO, ZONE
!
   implicit none
!
   integer, intent(in)  :: Igrid,n
   VTYPE,   intent(in)  :: Xsol(n)
   VTYPE, allocatable   :: z_loc(:), r_loc(:)
! $omp threadprivate(z_loc,r_loc)   
!..locals
   integer :: iel, ndof, k, i
! 
!--------------------------------------------------------------------------
!   
!..loop through the elements
!$omp parallel default(shared)   &
!$omp private(iel,ndof,k,i,z_loc,r_loc)
!$omp do schedule(dynamic)
   do iel = 1, GRID(Igrid)%nreles
!  ...number of dof
      ndof = GRID(Igrid)%dloc(iel)%ndof
      allocate(z_loc(ndof))
!  ...compute element solution
      do k = 1, ndof
         i  = GRID(Igrid)%dloc(iel)%lcon(k)
         z_loc(k) = Xsol(i)
      enddo
!      
!  ...compute the local residual vectors
#if C_MODE
      call ZHPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,z_loc,1,ZERO,GRID(Igrid)%loc(iel)%r,1)
#else
      call DSPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,z_loc,1,ZERO,GRID(Igrid)%loc(iel)%r,1)
#endif   
      GRID(Igrid)%loc(iel)%r = GRID(Igrid)%dloc(iel)%zbload - GRID(Igrid)%loc(iel)%r
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
!-----------------------------------------------------------------------
!
!   routine name       - update_local_residuals
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
   subroutine update_local_residuals(Igrid,z,n)
!
   use mg_data_structure, only: GRID
   use data_structure3D, only: ZERO, ZONE
!
   implicit none
!
   integer, intent(in)  :: Igrid,n
   VTYPE,   intent(in)  :: z(n)
   VTYPE, allocatable   :: z_loc(:), r_loc(:)
! $omp threadprivate(z_loc,r_loc)   
!..locals
   integer :: iel, ndof, k, i
! 
!--------------------------------------------------------------------------
!   
!..loop through the elements
!$omp parallel default(shared)   &
!$omp private(iel,ndof,k,i,z_loc,r_loc)
!$omp do schedule(dynamic)
   do iel = 1, GRID(Igrid)%nreles
!  ...number of dof
      ndof = GRID(Igrid)%dloc(iel)%ndof
      allocate(z_loc(ndof),r_loc(ndof))
!  ...compute element solution
      do k = 1, ndof
         i  = GRID(Igrid)%dloc(iel)%lcon(k)
         z_loc(k) = z(i)
      enddo
!      
!  ...compute the local residual vectors
#if C_MODE
      call ZHPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,z_loc,1,ZERO,r_loc,1)
#else
      call DSPMV('U',ndof,ZONE,GRID(Igrid)%dloc(iel)%zstiff,z_loc,1,ZERO,r_loc,1)
#endif   
      GRID(Igrid)%loc(iel)%r = GRID(Igrid)%loc(iel)%r - r_loc
!      
      deallocate(z_loc,r_loc)
!..end of loop through the elements
   enddo
!$omp end do
!$omp end parallel
!
!
   end subroutine update_local_residuals


!
!..compute number of dof for a coarse grid element
   subroutine compute_ndof(Igrid,Iel,Ndof)
!      
   use data_structure3D
   use mg_data_structure, only: GRID
!
   implicit none
!   
   integer, intent(in)  :: Igrid, Iel
   integer, intent(out) :: Ndof
   integer :: nrnodm, nod, i, j, k
   integer :: index(NRINDEX)
!
!-----------------------------------------------------------
!
! NEED TO MODIFY IN ORDER TO SUPPORT MULTIPLE VARIABLES 
! OF THE SAME ENERGY SPACE

   nrnodm = GRID(igrid)%constr(iel)%nrnodm 
   Ndof = 0
!..count coarse grid dof excluding the Dirichlet ones
   do i = 1, nrnodm
      nod = GRID(igrid)%constr(iel)%nodm(i)
      call get_index(nod,index)
      do k=1,NRINDEX
         select case(index(k))
         case(2)
!     ...update counter
            do j=1,GRID(igrid)%constr(iel)%ndofmH(i)
               Ndof = Ndof+1
            enddo
         case(4)
!     ...update counter
            do j=1,GRID(igrid)%constr(iel)%ndofmE(i)
               Ndof = Ndof+1
            enddo                  
         case(6)
!     ...update counter
            do j=1,GRID(igrid)%constr(iel)%ndofmV(i)
               Ndof = Ndof+1
            enddo                  
         end select
      enddo      
   enddo


   end subroutine compute_ndof







   subroutine allocate_local_vectors(Nrgrids)
!      
   use mg_data_structure, only: GRID
   use macro_grid_info,   only: NRDOF_MACRO
   use parameters, only : ZERO
!
   implicit none
!
   integer, intent(in) :: Nrgrids
   integer :: iel, ndof,ndof_c, nreles, igrid, ndof_glob
!
!---------------------------------------------------------
!
   allocate(G(nrgrids))

   do igrid = 1, Nrgrids
      nreles = GRID(Igrid)%nreles
      allocate(GRID(Igrid)%loc(nreles))

      do iel = 1, nreles
         ndof = GRID(Igrid)%dloc(iel)%ndof
         call compute_ndof(igrid,iel,ndof_c)
         GRID(igrid)%DLOC(iel)%ndof_c = ndof_c
         ndof_glob = NRDOF_MACRO(igrid)
!         
         allocate(GRID(Igrid)%loc(iel)%r(ndof))
         allocate(GRID(Igrid)%loc(iel)%z(ndof))
         allocate(GRID(Igrid)%loc(iel)%r_c(ndof_c))
         allocate(GRID(Igrid)%loc(iel)%z_c(ndof_c))

         GRID(Igrid)%loc(iel)%r = ZERO
         GRID(Igrid)%loc(iel)%r_c = ZERO
         GRID(Igrid)%loc(iel)%z = ZERO
         GRID(Igrid)%loc(iel)%z_c = ZERO
!
      enddo   
      allocate(G(igrid)%y1(ndof_glob))
      allocate(G(igrid)%y2(ndof_glob))
      allocate(G(igrid)%y3(ndof_glob))
   enddo   
!
   end subroutine allocate_local_vectors



   subroutine deallocate_local_vectors(Nrgrids)
!      
   use mg_data_structure, only: GRID
!
   implicit none
!
   integer, intent(in) :: Nrgrids
   integer :: iel, igrid
!
!------------------------------------------------------
!
   do igrid = 1, nrgrids
!
      do iel = 1, GRID(Igrid)%nreles
         deallocate(GRID(Igrid)%loc(iel)%r)
         deallocate(GRID(Igrid)%loc(iel)%z)
         deallocate(GRID(Igrid)%loc(iel)%r_c)
         deallocate(GRID(Igrid)%loc(iel)%z_c)
      enddo   
      deallocate(G(igrid)%y1)
      deallocate(G(igrid)%y2)
      deallocate(G(igrid)%y3)
      deallocate(GRID(igrid)%loc)
   enddo   
!
   deallocate(G)
!
   end subroutine deallocate_local_vectors
!



   subroutine deallocate_macro_system(Igrid)
!
   use mg_data_structure, only: GRID
!   
   implicit none
   integer :: Igrid, iel
!   
   do iel = 1, GRID(Igrid)%nreles
      deallocate(GRID(Igrid)%dloc(iel)%lcon)
      deallocate(GRID(Igrid)%dloc(iel)%zstiff)
      deallocate(GRID(Igrid)%dloc(iel)%zbload)
   enddo
   deallocate(GRID(Igrid)%dloc)

   end subroutine deallocate_macro_system
!
   end module pcg_info
