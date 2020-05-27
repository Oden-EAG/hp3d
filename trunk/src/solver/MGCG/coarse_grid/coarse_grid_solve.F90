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
#include "implicit_none.h"
!
   ! subroutine coarse_grid_solve(Igrid,z,m)
   subroutine coarse_grid_solve(Igrid)
!
   use data_structure3D,  only: ZERO
   use mg_data_structure, only: GRID, COARSE_SOLVER, MUMPS_SOLVER
   use macro_grid_info,   only: NRDOF_COARSE
   use stc,               only: CLOC
   use mumps,             only: MUMPS_PAR
   use pardiso_data
!
   implicit none
!
!-----------------------------------------------------------------------
!
   integer, intent(in) :: Igrid
   ! VTYPE,   intent(out):: z(m)
   VTYPE               :: rhs_coarse(NRDOF_COARSE)
   VTYPE, allocatable, save :: r_coarse(:), z_coarse(:), z_macro(:)
!$omp threadprivate(r_coarse,z_coarse,z_macro)
   integer             :: iel, k, i !,m
   integer             :: ndof_coarse, ndof_macro
!
!-----------------------------------------------------------------------
!
   rhs_coarse = ZERO !; z = ZERO
!..transfer the element residuals to the coarse grid
!
!$omp parallel default(shared)   &
!$omp private(iel,ndof_macro, ndof_coarse,k,i)
!$omp do reduction(+:rhs_coarse) schedule(dynamic)
   do iel = 1, GRID(Igrid)%nreles
      ndof_macro  = GRID(Igrid)%dloc(iel)%ndof ! macro dof
      ndof_coarse = CLOC(iel)%ni
      allocate(r_coarse(ndof_coarse))
!
      call macro2coarse(Igrid,iel,GRID(Igrid)%loc(iel)%r,ndof_macro,r_coarse,ndof_coarse)
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
   if (COARSE_SOLVER .eq. MUMPS_SOLVER) then
      mumps_par%RHS = rhs_coarse
!
      mumps_par%JOB = 3
#if C_MODE
      call zmumps(mumps_par)
#else
      call dmumps(mumps_par)
#endif
!
      rhs_coarse = mumps_par%RHS
   else
!
      PRDS_RHS = rhs_coarse
      PRDS_PHASE = 33 ! only solve
      PRDS_IPARM(8) = 1 ! max numbers of iterative refinement steps

      call pardiso(PRDS_PT,PRDS_MAXFCT,PRDS_MNUM,PRDS_MTYPE,PRDS_PHASE,   &
                   PRDS_N,PRDS_A(1:PRDS_NZ),PRDS_IA(1:NRDOF_COARSE+1),    &
                   PRDS_JA(1:PRDS_NZ),PRDS_PERM,PRDS_NRHS,PRDS_IPARM,     &
                   PRDS_MSGLVL,PRDS_RHS,PRDS_XSOL,PRDS_ERROR)

      rhs_coarse = PRDS_XSOL
   endif
!
!..restrict the element, prolong to the macro grid and assemble
!$omp parallel default(shared)   &
!$omp private(iel,ndof_macro, ndof_coarse,k,i)
!$omp do schedule(dynamic)
   do iel=1,GRID(Igrid)%nreles
!
      ndof_macro = GRID(Igrid)%dloc(iel)%ndof
      ndof_coarse = CLOC(iel)%ni
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
      call coarse2macro(Igrid,iel,z_coarse,ndof_coarse,z_macro,ndof_macro)
!
      deallocate(z_coarse)
!  ...assemble the correction
!
      GRID(Igrid)%loc(iel)%z = z_macro
!
      deallocate(z_macro)
!  ...end of loop through the elements
   enddo
!$omp end do
!$omp end parallel
!
!
   end subroutine coarse_grid_solve
