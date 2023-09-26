!----------------------------------------------------------------------------
!> @brief   initialize problem dependent environments, solvers, graphics and
!!          create initial mesh.
!> @date    Feb 2023
!----------------------------------------------------------------------------
subroutine initialize
   use environment
   use data_structure3D
   use refinements
   use control
   use GMP
   use common_prob_data
   use frsolmod
   use upscale
   use paraview
!
   implicit none
!
!--------------------------------------------------------------------------
!  G E O M E T R Y    S E T T I N G
   integer, parameter :: NDIM_PROB    = 3    ! Dimension
   integer, parameter :: MANDIM_PROB  = 3    ! Manifold dimension
                                             ! EXPECTED MAX NUMBER OF:
   integer, parameter :: MAXSU_PROB   = 1    ! surfaces
   integer, parameter :: MAXNP_PROB   = 500  ! points
   integer, parameter :: MAXNC_PROB   = 500  ! curves
   integer, parameter :: MAXTR_PROB   = 1    ! triangles
   integer, parameter :: MAXRE_PROB   = 100  ! rectangles
   integer, parameter :: MAXBT_PROB   = 1    ! prisms
   integer, parameter :: MAXHE_PROB   = 100  ! hexas
   integer, parameter :: MAXTE_PROB   = 1    ! tetras
   integer, parameter :: MAXPY_PROB   = 1    ! pyramids
!--------------------------------------------------------------------------
!  E Q U A T I O N    S E T T I N G
   integer, parameter :: NRCOMS_PROB  = 1     ! number of solution copies
   integer, parameter :: NRRHS_PROB   = 1     ! number of rhs
!--------------------------------------------------------------------------
   integer :: tmp
!--------------------------------------------------------------------------
!  output file open for the history of refinements (not for MPI)
   !call open_history_file(trim(FILE_HISTORY))
!
!  initialize refinements arrays
   call init_refinements(trim(FILE_REFINE))
!
!  generate constraint arrays - needed???
   call init_cnstr
!
!  read control file
   call read_control(trim(FILE_CONTROL))
!
!  read physics file quietly first to automatically setup equation settings
   call read_input(trim(FILE_PHYS))
!
!  setup equation settings based on physics data and common problem data
   call set_parameters(NRCOMS_PROB,NRRHS_PROB)
!
!  setup geometry settings
   call set_gmp_parameters( &
         NDIM_PROB, MANDIM_PROB, MAXSU_PROB,             &
         MAXNP_PROB, MAXNC_PROB, MAXTR_PROB, MAXRE_PROB, &
         MAXBT_PROB, MAXHE_PROB, MAXTE_PROB, MAXPY_PROB)
!
!  read geometry
   call read_geometry(trim(FILE_GEOM))
!
!  generate mesh and read physics file
   tmp = INTEGRATION
   call hp3gen(trim(FILE_PHYS))
   INTEGRATION = tmp
!
end subroutine initialize
