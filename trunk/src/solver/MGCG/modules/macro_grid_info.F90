!----------------------------------------------------------------------
!
!   module name        - macro_grid_info
!
!----------------------------------------------------------------------
!
!   latest revision    - FEB 18
!
!   purpose            - holds all required info for the macro-grid
!
!----------------------------------------------------------------------
!
   module macro_grid_info
!    
   use derived_types  
!
   integer  :: NRELES_COARSE
   integer  :: NRDOF_COARSE  
   integer, allocatable :: NRDOF_MACRO(:)
!
   type (sol_struct),allocatable :: ZSOL_C(:)
!
   type(store_mat), allocatable :: A_MACRO(:)
!
!
   end module macro_grid_info