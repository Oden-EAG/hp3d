!-----------------------------------------------------------------------
! 
!    routine name       - debug_driver
! 
!-----------------------------------------------------------------------
! 
!    latest revision    - Jan 2018
! 
!    purpose            - debugging driver for prolongation module          
! 
!    arguments          - none
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 
   subroutine debug_driver(nrule)
!
   use data_structure3D, only: NRELES 
   use mg_data_structure
   use common_prob_data
   use macro_grid_info,  only: NRELES_COARSE, MDLE_MACRO, MACRO_ELEM,ZSOL_C
   use mumps
   use constraints_info
   use patch_info
   use m_assembly,       only: CLOC
!
   IMPLICIT NONE
! 
   integer, intent (in) :: nrule 
   integer              :: nreflag, nstop, iel, mdle, i
   real*8               :: factor
!   
   integer              :: iprint, idec
   real*8               :: tm(10), start, finish, omp_get_wtime


!-----------------------------------------------------------------------
!
!
   iprint = 0
!  
   call mg_init
!
   nreflag = 1; factor = 0.1d0
   ! call refine_DPG(IUNIFORM,1,factor, nstop)
   ! call mumps_interface('H')   
   ! call refine_DPG(IUNIFORM,1,factor, nstop)
!..uniform refinement
   do i = 1,1
      ! call refine_DPG(IADAPTIVE,3,factor, nstop)
   enddo 
!
!..save mdle nodes of the coarse grid
!  
   NRELES_COARSE = NRELES
! 
   allocate(MACRO_ELEM(NRELES_COARSE))
   allocate(MDLE_MACRO(NRELES_COARSE))
   mdle = 0
   do iel = 1, NRELES_COARSE
      call nelcon(mdle,mdle)
      MDLE_MACRO(iel) = mdle
      NODES_MG(mdle)%iel = iel
   enddo
!
!..mark current mesh
   call mark_masters    
!
! 
   write(*,*) 'Solve coarse grid system...'
   call coarse_solve_mumps('H')
!
!..THIS WILL BE SUBSTITUTED BY THE TRUE MAX GENERATION
   do i = 1,6
   MAXGEN_PR = i
      ! do i = 1,MAXGEN_PR
      call refine_DPG(IADAPTIVE,nrule,factor, nstop)

      ! enddo     
      ! write(*,*) 'Solve coarse grid system...'
      ! call coarse_solve_mumpsMG('H')
   !
   !   
      write(*,*) 'Create patch data structure...'
      start = omp_get_wtime()
      call patch_nodes
      finish = omp_get_wtime()

      write(*,1000) finish-start
   1000 format(' Create patch data structure:t=    ', es13.4)    

   !
   !
      call allocate_constr
      write(*,*) 'Store prolongation coefficients...'
      start = omp_get_wtime()
      call store_constraints
      finish = omp_get_wtime()

      write(*,1001) finish-start
   1001 format(' Store prolongation coefficients:t=', es13.4)    
   !
   !

      write(*,*) 'Assemble macro-grid system...'
      call macro_assembly

   !
      ! write(*,*) 'Calculate prolongation error...'
      ! call prolong_error
   ! 
      call deallocate_constr
      ! call deallocate_coarse_grid

   enddo   
   !
   !
      call deallocate_coarse_grid
   !
      call mg_finalize
      call refine_DPG(0,nrule,factor, nstop)

      deallocate(MDLE_MACRO)
      deallocate(MACRO_ELEM)

      ! call graphb

   ! call mumps_interface('H')   
!
!
   end subroutine debug_driver


   subroutine deallocate_coarse_grid
   use mumps
   use macro_grid_info,  only: NRELES_COARSE, MDLE_MACRO, MACRO_ELEM,ZSOL_C
   use patch_info
   use m_assembly,       only: CLOC

   integer :: iel
!
   call mumps_destroy
   do iel = 1, NRELES_COARSE
      deallocate(CLOC(iel)%con)
      deallocate(ZSOL_C(iel)%coarse)
   enddo   
   deallocate(ZSOL_C)
   deallocate(CGRID_VERTICES)
   deallocate(PTCH)
   deallocate(CLOC)
   

   end subroutine deallocate_coarse_grid