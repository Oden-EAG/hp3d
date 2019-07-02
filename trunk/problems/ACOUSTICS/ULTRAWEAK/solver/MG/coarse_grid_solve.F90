!   
!-----------------------------------------------------------------------
!
!   routine name       - coarse_grid_solve
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2018
!
!   purpose            - perform the coarse grid solve of 
!                        the multigrid solver
!   arguments :
!
!              out 
!                 z    - the correction to the solution vector
!
!-----------------------------------------------------------------------
! 
   subroutine coarse_grid_solve(z)
!
   use data_structure3D, only: ZERO
   use macro_grid_info,  only: NRDOF_MACRO, NRELES_COARSE, DLOC, NRDOF_COARSE
   use pcg_info,         only: LOC
   use m_assembly,       only: CLOC
   use mumps,            only: MUMPS_PAR
!
   implicit none
!   
!-----------------------------------------------------------------------
!
#if C_MODE
   complex*16, intent(out) :: z(NRDOF_MACRO) 
   complex*16              :: rhs_coarse(NRDOF_COARSE)
   complex*16, allocatable, save :: r_coarse(:), z_coarse(:), z_macro(:)
#else
   real*8,     intent(out) :: z(NRDOF_MACRO) 
   real*8                  :: rhs_coarse(NRDOF_COARSE)
   real*8,     allocatable, save :: r_coarse(:), z_coarse(:), z_macro(:)
#endif 
!$omp threadprivate(r_coarse,z_coarse,z_macro)   
   integer                 :: iel, k, i
   integer                 :: ndof_coarse, ndof_macro
!
!-----------------------------------------------------------------------
!
   z =ZERO ; rhs_coarse = ZERO
!..transfer the element residuals to the coarse grid
!
!$omp parallel default(shared)   &
!$omp private(iel,ndof_macro, ndof_coarse,k,i)
!$omp do reduction(+:rhs_coarse) schedule(dynamic)
   do iel = 1, NRELES_COARSE
      ndof_macro  = DLOC(iel)%ndof ! macro dof
      ndof_coarse = CLOC(iel)%nHc + CLOC(iel)%nEc + CLOC(iel)%nVc
      allocate(r_coarse(ndof_coarse))
!
      call macro2coarse(iel,LOC(iel)%r,ndof_macro,r_coarse,ndof_coarse)
!
!  ...assemble global coarse residual
      do k=1,ndof_coarse
!     ...global coarse dof is:
         i = CLOC(iel)%con(k)
!     ...assemble global load vector
         rhs_coarse(i) =  rhs_coarse(i) + r_coarse(k)
      enddo
      deallocate(r_coarse)
!
   enddo
!$omp end do
!$omp end parallel
!
   mumps_par%RHS = rhs_coarse
!
!..solve the coarse grid system
!
   mumps_par%JOB = 3
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
!
   rhs_coarse = mumps_par%RHS
!
!..restrict the element, prolong to the macro grid and assemble
!$omp parallel default(shared)   &
!$omp private(iel,ndof_macro, ndof_coarse,k,i)
!$omp do schedule(dynamic)
   do iel=1,NRELES_COARSE
!
      ndof_macro = DLOC(iel)%ndof
      ndof_coarse = CLOC(iel)%nHc + CLOC(iel)%nEc + CLOC(iel)%nVc
!
      allocate(z_coarse(ndof_coarse)) 
      allocate(z_macro(ndof_macro)) 
!
      do k=1,ndof_coarse
         i = CLOC(iel)%con(k)
         z_coarse(k) = rhs_coarse(i)
      enddo
!      
!  ...apply the prolongation operator
      call coarse2macro(iel,z_coarse,ndof_coarse,z_macro,ndof_macro)
!
      deallocate(z_coarse)
!  ...assemble the correction
      do k=1,ndof_macro
         i = DLOC(iel)%lcon(k)
         z(i) = z_macro(k)
      enddo
!
      deallocate(z_macro)
!  ...end of loop through the elements
   enddo
!$omp end do
!$omp end parallel
!
! 
   end subroutine coarse_grid_solve