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
   use data_structure3D
   use mg_data_structure, only: MAXGEN_PR
   use parameters
!
   integer  :: NRELES_COARSE
   integer  :: NRDOF_MACRO, NRDOF_COARSE         
!
   integer, allocatable :: MDLE_MACRO(:)  
!

   type super_vector
      integer              :: ndof_macro, nrnod_macro
      integer, allocatable :: nod_macro(:)  , ndofH_macro(:)
      integer, allocatable :: ndofE_macro(:), ndofV_macro(:)
   end type super_vector
!..macro element nod list and number of nodes
   type(super_vector), allocatable :: MACRO_ELEM(:)
!
!..workspace for prolong solution
   
   type sol_struct
      integer    :: ndof_coarse
      integer    :: nedge_orient(12),nface_orient(6), norder(19)
      real*8     :: xnod(NDIMEN,MAXbrickH)

#if C_MODE
      complex*16, allocatable :: zdofH(:,:),zdofE(:,:)
      complex*16, allocatable :: zdofV(:,:),zdofQ(:,:)
      complex*16, allocatable :: coarse(:)
      complex*16, allocatable :: macro(:)  
#else 
      real*8,     allocatable :: zdofH(:,:),zdofE(:,:)
      real*8,     allocatable :: zdofV(:,:),zdofQ(:,:)
      real*8,     allocatable :: coarse(:)
      real*8,     allocatable :: macro(:)  
#endif
   endtype sol_struct

   type (sol_struct),allocatable :: ZSOL_C(:)
!
!
!..workspace for celem_macro
!..store matrices and connectivities for the macro element
   type con_info
!  ...local connectivity map
      integer, allocatable :: lcon(:)
!  ...size
      integer :: ndof
!  ...mdle number
      integer :: mdle      
!  ...element matrices
#if C_MODE
      complex*16, allocatable :: zstiff(:), zbload(:)
#else 
      real*8,     allocatable :: zstiff(:), zbload(:)
#endif
!    
   end type con_info
!
   type(con_info), allocatable :: DLOC(:)


   type store_mat
      type(con_info), allocatable :: GLOC(:)

!  ...local matrix (Abb)^-1 * Abi and vector (Abb)^-1 * bi
#if C_MODE
      complex*16, allocatable :: array(:,:), vect(:)
#else 
      real*8,     allocatable :: array(:,:), vect(:)
#endif
!  ...size
      integer :: ni, nb
!  ...number of mdle nodes within a macro_element      
      integer :: nrmdle   
      end type store_mat
!
      type(store_mat), allocatable :: A_MACRO(:)


   end module macro_grid_info