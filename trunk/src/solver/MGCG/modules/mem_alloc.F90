!----------------------------------------------------------------------
!
!   module name        - mem_alloc
!
!----------------------------------------------------------------------
!
!   latest revision    - MAR 18
!
!   purpose            - modules handles memory allocation for mg solver
!
!----------------------------------------------------------------------
!
   module mem_alloc
!  
   use patch_info,        only: deallocate_patches
   use constraints_info,  only: deallocate_constr
   use pcg_info,          only: deallocate_macro_system
   use mg_data_structure, only: NRGRIDS
!
   contains
!
   subroutine mg_dealloc
!
   implicit none

   integer :: igrid

   write(*,*) 'mg_dealloc: deallocating patches, constraints and macro system'
   do igrid = 1,NRGRIDS
      call deallocate_constr(igrid)
      call deallocate_patches(Igrid)   
      call deallocate_macro_system(Igrid)
   enddo   
! 
!..release memory for the coarse grid
   write(*,*) 'mg_dealloc: deallocating coarse grid'
   call deallocate_coarse_grid

   end subroutine mg_dealloc
!
!
   subroutine deallocate_coarse_grid
!   
   use mumps
   use pardiso_data
   use mg_data_structure
   use patch_info
   use stc,       only: CLOC
!
   implicit none

   integer :: iel
!


   select case(COARSE_SOLVER)

   case(MUMPS_SOLVER)
      call mumps_destroy
   case(PARDISO_SOLVER)
      call finalize_pardiso
   end select   

!
   do iel = 1, GRID(1)%nreles
      deallocate(CLOC(iel)%con)
   enddo   
!
   deallocate(CGRID_VERTICES)
   deallocate(GRID(1)%ptch)
   deallocate(CLOC)
!   
   end subroutine deallocate_coarse_grid
!
!
   end module mem_alloc