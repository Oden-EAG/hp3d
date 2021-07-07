module PROJ
!      
      save
!     
!     number of right hand sides      
      integer, parameter :: NUMBER_OF_RHS=1
!
!     0 - interactive mode 
!     1 - test multiple projections
!     2 - test orientations        
!     3 - test constrained approximation
!     4 - test refinements for deadlock
      integer :: ITEST
!
!     0 - frontal solver
!     1 - MUMPS
      integer :: ISOLVER
!
!     isotropic orders of approximation for prism, hexa, tet, pyramid     
      integer :: IORDER_PRIS
      integer :: IORDER_BRIC
      integer :: IORDER_TETR
      integer :: IORDER_PYRA 
!
!     0 - lower  order modes for H1, H(curl), H(div), L2
!     1 - higher order modes for H1     
!     2 - higher order modes for H(curl)
!     3 - higher order modes for H(div)
      integer :: NEDELEC
!
!     0     - set all components of exact solution
!     1,2,3 - set only one component of exact solution (remaining ones are null)     
      integer :: ICOMP
!
!     exponents of monomials representing the exact solution
      integer :: NPX,NPY,NPZ
!
!     'pris' - solve projections on master prism
!     'hexa' - solve projections on master hexa
!     'tetr' - solve projections on master tet
!     'pyra' - solve projections on master pyramid
!     'hybr' - solve projections on hybrid mesh (prism + hexa + tet + pyramid)
      character(len=4) :: PROBLEM
!
!     percentage (0,1) of elements to refine when checking for deadlock
      real*8 :: PERC
!
!     number of refinement levels to perform when checking for deadlock
      integer :: NLEVEL
!
!     case tag
      integer :: ITAG = 0
!      
!     tags file
      character(len=128) :: FILE_TAGS
!
!
endmodule PROJ
