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
   subroutine store_elem_sol(Iel, Zelem, Ndof)
! 
   use data_structure3D
   use macro_grid_info
!
   IMPLICIT NONE
!
   integer, intent(in) :: Iel, Ndof
#if C_MODE
   complex*16, intent(in) :: Zelem(Ndof)
#else
   real*8,     intent(in) :: Zelem(Ndof)
#endif   

   ZSOL_C(iel)%ndof_coarse = Ndof
   allocate(ZSOL_C(iel)%coarse(Ndof))
   ZSOL_C(iel)%coarse = Zelem

   call find_order(MDLE_MACRO(iel), ZSOL_C(iel)%norder)
   call find_orient(MDLE_MACRO(iel),ZSOL_C(iel)%nedge_orient,ZSOL_C(iel)%nface_orient)
   call nodcor(MDLE_MACRO(iel),ZSOL_C(iel)%xnod)

!..determine element dof
   allocate(ZSOL_C(iel)%zdofH(MAXEQNH,MAXbrickH))
   allocate(ZSOL_C(iel)%zdofE(MAXEQNE,MAXbrickE))
   allocate(ZSOL_C(iel)%zdofV(MAXEQNV,MAXbrickV))
   allocate(ZSOL_C(iel)%zdofQ(MAXEQNQ,MAXbrickQ))
!
   call get_sol_dof(Iel, ZSOL_C(iel)%zdofH, ZSOL_C(iel)%zdofE, ZSOL_C(iel)%zdofV) 
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
   subroutine prolong_solution
! 
   use macro_grid_info
! 
   IMPLICIT NONE
!
   integer :: iel, ndof_coarse, ndof_macro, i
!
!----------------------------------------------------------------------
!

!$omp parallel default(shared)  &
!$omp private(iel, ndof_coarse, ndof_macro)
!$omp do schedule(guided)
   do iel = 1, NRELES_COARSE
      ndof_coarse = ZSOL_C(iel)%ndof_coarse
      ndof_macro  = MACRO_ELEM(iel)%ndof_macro
!
      allocate(ZSOL_C(iel)%macro(ndof_macro))
!
      call coarse2macro(iel, ZSOL_C(iel)%coarse, ndof_coarse,     &
                             ZSOL_C(iel)%macro,  ndof_macro)
!
!                          
   enddo   
!$omp end do   
!$omp end parallel
!
!
   end subroutine prolong_solution
!
