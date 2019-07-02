c----------------------------------------------------------------------
c
c   module name        - assembly
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - module defines element arrays
c
c
c-----------------------------------------------------------------------
c
      module assembly
      save
c     
c  ...number of right-hand sides (load vectors)
      integer :: NR_RHS
c
c  ...maximum number of element local dof for a physics attribute
      integer, allocatable :: MAXDOFS(:)
c
c  ...maximum number of modified element dof for all variables
      integer :: MAXDOFM
c
c  ...maximum number of modified element dof for all variables after 
c     compresssion
      integer :: MAXDOFC
c
c
c
c**********************************************************************      
c  ...PG, Nov 10: this type is not used ANYWHERE in the code!!!
c                 consider removing it!!!      
      type super_vector
c
c  ...number of rows
      integer :: nrow
#if C_MODE
      complex*16, allocatable :: vect(:)
#else
      double precision, allocatable :: vect(:)
#endif
      endtype super_vector
c**********************************************************************      
c
c      
c
      type super_array
c
c  ...number of columns
      integer :: ncol
c
c  ...number of rows
      integer :: nrow
#if C_MODE
      complex*16, allocatable :: array(:,:)
#else
      double precision, allocatable :: array(:,:)
#endif
      endtype super_array
c
c  ...element local load vectors
      type(super_array), allocatable  :: BLOC(:)
c
c  ...element local stiffness matrices
      type(super_array), allocatable  :: ALOC(:,:)
c
c  ...half-way modified element stifness matrix
      type(super_array), allocatable  :: AAUX(:)
c
c  ...modified element load vector and stifness matrix
#if C_MODE
      complex*16, allocatable :: ZBMOD(:,:)
      complex*16, allocatable :: ZAMOD(:,:)
#else
      double precision, allocatable :: ZBMOD(:,:)
      double precision, allocatable :: ZAMOD(:,:)
#endif
c
c  ...dof extraction vector
      integer, allocatable :: NEXTRACT(:)
c
c  ...Dirichlet data flags
      integer, allocatable :: IDBC(:)
c
c  ...Dirichlet data
#if C_MODE
      complex*16, allocatable :: ZDOFD(:)
#else
      double precision, allocatable :: ZDOFD(:)
#endif
c
c      
      end module assembly
