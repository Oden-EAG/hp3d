
! -----------------------------------------------------------------------
!
!    routine name       - pcg_solve
!
! -----------------------------------------------------------------------
!
!    latest revision    - May 2018
!
!    purpose            - routine interfaces with the CG solver and 
!                         after the solution stores the dof in the 
!                         data structure
!    Arguments
!    in  :        Igrid - the mesh index
!
!-----------------------------------------------------------------------
!
#include "implicit_none.h"
!
   subroutine pcg_solve(Igrid)

   use mg_data_structure,  ONLY: GRID, NRGRIDS
   use pcg_info,           ONLY: XSOL
   use patch_info,         ONLY: deallocate_patches
   use parameters,         ONLY: ZERO
   use macro_grid_info,    ONLY: A_MACRO, ZSOL_C, NRDOF_MACRO
!   
!-----------------------------------------------------------------------
!
   implicit none
!
   integer, intent(in) :: Igrid
   integer :: iel,ndof,k,i
!   
   VTYPE, allocatable, save :: zsol_macro(:)
!$omp threadprivate(zsol_macro)   
!-----------------------------------------------------------------------
!
   allocate(XSOL(NRDOF_MACRO(Igrid))) ; XSOL = ZERO
!
   call prolong_solution(Igrid)
!
   do iel = 1, GRID(Igrid)%nreles
      ndof = GRID(Igrid)%dloc(iel)%ndof
!
      do k = 1,ndof
         i = GRID(Igrid)%dloc(iel)%lcon(k)
         XSOL(i) = ZSOL_C(iel)%macro(k)
      enddo
      deallocate(ZSOL_C(iel)%macro)
!..end of loop through macro elements
   enddo
   deallocate(ZSOL_C)

   if (igrid .eq. 1) XSOL = ZERO
!
   call mg_pcg_solve(Igrid)
! 
!..store the solution into the data structure
!..loop through macro elements
!$omp parallel default(shared)    &
!$omp private(iel,ndof,k,i)
!$omp do schedule(guided)
   do iel = 1, GRID(Igrid)%nreles
      ndof = GRID(Igrid)%dloc(iel)%ndof
!
      allocate(zsol_macro(ndof)) 
!
      do k=1,ndof
         i = GRID(Igrid)%dloc(iel)%lcon(k)
         zsol_macro(k) = XSOL(i)
      enddo
!
      call solout_macro(Igrid, iel,zsol_macro,ndof)
      deallocate(A_MACRO(iel)%GLOC)
      deallocate(zsol_macro) 
!..end of loop through macro elements
   enddo
!$omp end do   
!$omp end parallel
!
   deallocate(XSOL,A_MACRO)
!   
   if (Igrid .lt. NRGRIDS) then
      call compute_sol_dof(Igrid+1)
   endif   
!
!
   end subroutine pcg_solve
