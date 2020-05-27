!----------------------------------------------------------------------
!
!   module name        - derived_types
!
!----------------------------------------------------------------------
!
!   latest revision    - Sept 18
!
!   purpose            - modules defined user derived types used
!                        in the multigrid solver
!
!----------------------------------------------------------------------
!
#include "implicit_none.h"

   module derived_types
!
   use parameters
!
!---------------------------------------------------------------------
!   
!..node_mg data structure
   type node_mg
!  ...order of a nod on the coarse grid
      integer, allocatable :: orderC(:)
!  ...master flag      
      integer, allocatable :: master(:)
!
!  ...natural ordering of elements
      integer, allocatable :: iel(:)
!  ...number of dof supported by the node
      integer, allocatable :: nod_ndof(:)
!  ...visitation flag      
      integer :: visit      
!            
   end type node_mg
!
!---------------------------------------------------------------------
!      
   type super_vector
      integer              :: ndof_macro, nrnod_macro
      integer, allocatable :: nod_macro(:)  , ndofH_macro(:)
      integer, allocatable :: ndofE_macro(:), ndofV_macro(:)
   end type super_vector
!
!---------------------------------------------------------------------
! 
!..workspace for prolong solution
   type sol_struct
      integer    :: ndof_coarse
      integer    :: nedge_orient(12),nface_orient(6), norder(19)
      real*8     :: xnod(NDIMEN,MAXbrickH)

      VTYPE, allocatable :: zdofH(:,:),zdofE(:,:)
      VTYPE, allocatable :: zdofV(:,:),zdofQ(:,:)
      VTYPE, allocatable :: coarse(:)
      VTYPE, allocatable :: macro(:)  
   endtype sol_struct
!
!---------------------------------------------------------------------
! 
!..workspace for celem_macro
!..store matrices and connectivities for the macro element
   type con_info
!  ...local connectivity map
      integer, allocatable :: lcon(:)
!  ...size
      integer :: ndof, ndof_c
!  ...mdle number
      integer :: mdle      
!
      integer :: iel      
!  ...element matrices
      VTYPE, allocatable :: zstiff(:), zbload(:)

!  ...Schur complement corresponding to element interior dof
!    
   end type con_info
!
!---------------------------------------------------------------------
! 
  type con_info2
!  ...local connectivity map
      integer, allocatable :: lcon(:)
!  ...size
      integer :: mdle     
      integer :: ndof
!
      integer :: iel

!  ...Schur complement corresponding to element interior dof
      VTYPE, allocatable :: array(:,:), vect(:)
      integer :: ni
!    
   end type con_info2


   type store_mat
!
      type(con_info2), allocatable :: GLOC(:)

!  ...local matrix (Abb)^-1 * Abi and vector (Abb)^-1 * bi
      VTYPE, allocatable :: array(:,:), vect(:)
!  ...size
      integer :: ni, nb
!  ...number of mdle nodes within a macro_element      
      integer :: nrmdle   
   end type store_mat
!
!---------------------------------------------------------------------
! 
   integer, parameter :: MAX_NREDGES = 150
   integer, parameter :: MAX_NRFACES = 150
   integer, parameter :: MAX_PATCH_MDLE = 100
!   
   type patch1
!   
!  ...mdle number
      integer :: mdle
!  ...number of local dof for each element in the patch
      integer :: ndof
!  ...connectivity arrays for each mdle node in the patch
      integer, allocatable :: lcon(:)
   end type patch1
!
!..patch information   
   type patch
!  ...number of coarse grid edges and faces in the patch
      integer :: nredges
      integer :: nrfaces
      integer :: nsz    
!  ...edge and face lists (these are coarse grid nodes)
      integer :: nedgel(MAX_NREDGES), nfacel(MAX_NRFACES)    
!  ...the list of nodes for each patch
      integer, allocatable :: nodl(:)
!
!  ...number of mdle nodes contributing to the patch
      integer :: nrmdle
!  ...list of mdle nodes contributing to the patch
      type(patch1) :: mdlel(MAX_PATCH_MDLE)
!  ...total number of dof in the patch      
      integer :: nrdof
!  ...connectivity map form patch to global
      integer, allocatable :: lcon(:)
!  ...Cholesky decomposition of the assembled block in Lapack packed form (dense matrix)
      VTYPE,  allocatable :: zAp(:) 
!  
   end type patch
!
!---------------------------------------------------------------------
! 
   type constraint
!  ...modified coarse grid element
      integer              :: nrnodm, nodm(MAXNODM), nrdof
      integer              :: ndofmH(MAXNODM), ndofmE(MAXNODM), ndofmV(MAXNODM)
!  ...macro-element
      integer              :: nrnod_macro, nrdofH_macro, nrdofE_macro, nrdofV_macro
      integer, allocatable :: nod_macro(:)
      integer, allocatable :: ndofH_macro(:),  ndofE_macro(:),  ndofV_macro(:)
      integer, allocatable :: nrconH_macro(:), nrconE_macro(:), nrconV_macro(:)
      integer, allocatable :: nacH_macro(:,:), nacE_macro(:,:), nacV_macro(:,:)
      real*8,  allocatable :: constrH_macro(:,:),constrE_macro(:,:), constrV_macro(:,:)

   end type constraint
!
!---------------------------------------------------------------------
! 
!..local data needed by the solver
   type sarray
      VTYPE, allocatable :: z(:), z_c(:) ,r(:), r_c(:)
   endtype sarray
!
!
!..required structure for multigrid
   type schur_compl
      type(con_info), allocatable :: GLOC(:)
!                 
!  ...local matrix A22 in sparse form 
      integer, allocatable :: IA(:), JA(:)
      VTYPE,   allocatable :: SA(:)
!  ...local vector  A22^(-1) * A21
      VTYPE,   allocatable :: A21(:,:), r2(:)
! ....size
      integer n1, n2, inz
!  ...number of mdle nodes within a macro_element      
      integer nrmdle
   end type schur_compl
!
!      
   type ssarray
!   
!  ...Number of macro elements
      integer :: nreles
!      
!  ...Number of patches  
      integer :: nrpatch
!      
!  ...list of macro mdle nodes   
      integer, allocatable :: mdlel(:)
!
!  ...macro element nod list and number of nodes
      type(super_vector), allocatable :: macro_elem(:)
!  
!  ...matrices and connectivities for the macro element
      type(con_info), allocatable :: dloc(:)
!
!  ...constraints coefficients for the prolongation
      type (constraint), allocatable :: constr(:)
!
!  ...patch structure required for smoothing
      type(patch), allocatable :: ptch(:)
!
!  ...local arrays for element residual and solution vectors during solution step
      type(sarray), allocatable :: loc(:)
!
!  ...structure holding the Schur complement extension operator
      type(schur_compl), allocatable :: sch(:)
!
!
   endtype ssarray
!
!
   end module derived_types
