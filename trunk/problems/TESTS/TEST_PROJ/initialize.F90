subroutine initialize
!      
      use PROJ , only : NUMBER_OF_RHS
      use environment
      use control
      use parameters
      use gmp
      use frsolmod
      use data_structure3D
      use refinements
!!!      use uhm2
!
      implicit none
      character(len=1024) :: argv
!      
!-----------------------------------------------------------------------      
!     
!  ...UHM SOLVER........................................................
!!      UHM_VERBOSE            = .FALSE.
!!      UHM_SOLVER_REPORT_SHOW = .FALSE.
!    
!!      call get_command(argv)
!!      call uhm_initialize(argv)
!    
!!      call uhm_option_begin
!!      call uhm_option_end
!      
!!      call uhm_time_begin(UHM_TIME_LEVEL_FRONT)
!      call uhm_direct_lu_piv_initialize( &
!           UHM_DOUBLE, NUMBER_OF_RHS, 256, UHM_SOLVER_PTR)
!           
!  ...INPUT FILES.......................................................
      call open_history_file(trim(FILE_HISTORY))
      call init_refinements( trim(FILE_REFINE) )
!      
      call read_control(     trim(FILE_CONTROL))    
!
!  ...GMP DATA STRUCTURE................................................
!
!                               NDIM //   MANDIM //    
      call set_gmp_parameters(     3 ,         3 ,               &
!                              MAXSU //    MAXNP //    MAXNC //  
                                  30 ,     10000 ,     10000 ,   &
!                              MAXTR //    MAXRE //    MAXBT //  
                               10000 ,     10000 ,     10000 ,   &
!                              MAXHE //    MAXTE //    MAXPY //  
                               10000 ,     10000 ,     10000    )
!
!  ...HP3D DATA STRUCTURE...............................................
!=======================================================================      
! DS PARAMETERS ( aka "The product of a sick mind") :                  |
!                                                                      |
!    NRCOMS  - # of copies of the DS, i.e., the bloody # of rhs        |  
!                                                                      |
!    MAXEQNH - maximum # of equations posed in H1                      |
!      Ex: acoustic & elasticity. (3 + 1) * # of rhs                   |
!                                                                      |
!    MAXEQHE,V,Q - same as above                                       |
!                                                                      |
!    MAXNRHS - maximum number of rhs's (yet again...)                  |
!                                                                      |
!  Constraints on parameters to avoid catastrophe:                     |
!                                                                      |
!    MAX...      >= 1                                                  |
!    MAXNRHS     >= NRCOMS                                             |
!    MAXEQNH     >= NRCOMS * # of equations posed in H1                |
!    MAXEQNE,V,Q >= NRCOMS * # of equations posed in H(curl),H(div),L2 |
!=======================================================================
! 
!                                 NRCOMS //       MAXNRHS //
      call set_parameters( NUMBER_OF_RHS ,  NUMBER_OF_RHS ,             &
!                                MAXEQNH //       MAXEQNE // MAXEQNV // MAXEQNQ //
                                       1 ,              1 ,        1 ,        1    )
!     
      call read_geometry(trim(FILE_GEOM))
      call print_GMP_parameters
!     
!     frontal solver workspace
      call set_frsol_workspace(10000000)
!      
!     initialize parameters for X11 graphics      
!                              MXIGTR // MXIGSTR // MXRGTRZ 
      call set_x11_workspace( 1000000 ,  1000000 ,  100000  )
!      
!     generate mesh
      call hp3gen(trim(FILE_PHYS))
!      
!     mesh refinement filters
      call set_mesh_refinement
!
!      
endsubroutine initialize
