!-----------------------------------------------------------------------
!
!    routine name       - store_elem_sol
!
!-----------------------------------------------------------------------
!
!    latest revision    - Mar 2018
!
!    purpose            - store element modified solution vector,
!                         order and geometry dof
!
!----------------------------------------------------------------------
!
#include "implicit_none.h"
!
   subroutine store_elem_sol(Iel, Zelem, Ndof)
!
   use macro_grid_info,   only: ZSOL_C
!
   implicit none
!
   integer, intent(in) :: Iel, Ndof
   VTYPE,   intent(in) :: Zelem(Ndof)
   integer :: mdle

   ZSOL_C(iel)%ndof_coarse = Ndof
   allocate(ZSOL_C(iel)%coarse(Ndof))
   ZSOL_C(iel)%coarse = Zelem
!
!
   end subroutine store_elem_sol
!
!
!-----------------------------------------------------------------------
!
!    routine name       - prolong_solution
!
!-----------------------------------------------------------------------
!
!    latest revision    - Feb 2018
!
!    purpose            - driver for the prolongation of the coarse solution
!                         to the macro grid. Used for debugging purposes.
!
!----------------------------------------------------------------------
!
   subroutine prolong_solution(Igrid)
!
!
   use mg_data_structure, only: GRID
   use macro_grid_info,   only: ZSOL_C


   implicit none
!
   integer :: Igrid
   integer :: iel, ndof_coarse, ndof_macro, i
!
!----------------------------------------------------------------------
!

!$omp parallel default(shared)  &
!$omp private(iel, ndof_coarse, ndof_macro)
!$omp do schedule(guided)
!
   do iel = 1, GRID(Igrid)%nreles
      ndof_coarse = ZSOL_C(iel)%ndof_coarse
      ndof_macro  = GRID(Igrid)%macro_elem(iel)%ndof_macro
!
      allocate(ZSOL_C(iel)%macro(ndof_macro))
!
      call coarse2macro(Igrid,iel, ZSOL_C(iel)%coarse, ndof_coarse,     &
                                   ZSOL_C(iel)%macro,  ndof_macro )
!
      deallocate(ZSOL_C(iel)%coarse)
!
   enddo
!$omp end do
!$omp end parallel
!
!
   end subroutine prolong_solution
!
