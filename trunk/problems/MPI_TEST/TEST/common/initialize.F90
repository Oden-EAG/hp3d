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
  use gmp
  use common_prob_data
  use testvars
  use frsolmod
  use upscale
  use paraview
  use zoltan_wrapper
  use mpi_param, only: NUM_PROCS
  !
  implicit none
  !
  !--------------------------------------------------------------------------
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB    = 3    ! Dimension
  integer, parameter :: MANDIM_PROB  = 3    ! Manifold dimension
                                            ! EXPECTED MAX NUMBER OF:
  integer, parameter :: MAXSU_PROB   = 5    ! surfaces
  integer, parameter :: MAXNP_PROB   = 512  ! points
  integer, parameter :: MAXNC_PROB   = 1344 ! curves
  integer, parameter :: MAXTR_PROB   = 1    ! triangles
  integer, parameter :: MAXRE_PROB   = 1176 ! rectangles
  integer, parameter :: MAXBT_PROB   = 1    ! prisms
  integer, parameter :: MAXHE_PROB   = 342  ! hexas
  integer, parameter :: MAXTE_PROB   = 1    ! tetras
  integer, parameter :: MAXPY_PROB   = 1    ! pyramids
  !--------------------------------------------------------------------------
  ! E Q U A T I O N    S E T T I N G
  integer, parameter :: NRCOMS_PROB  = 1     ! number of component
  integer :: MAXNRHS_PROB                    ! MAX nr rhs
  integer :: MAXEQNH_PROB                    ! MAX H1 var
  integer :: MAXEQNE_PROB                    ! MAX Hcurl var
  integer :: MAXEQNV_PROB                    ! MAX Hdiv var
  integer :: MAXEQNQ_PROB                    ! MAX L2 var
  !--------------------------------------------------------------------------
  integer :: iflag,i,INTEGRATION_tmp
  character(len=1024) :: argv
  !--------------------------------------------------------------------------
  ! output file open for the history of refinements (not for MPI)
  !call open_history_file(trim(FILE_HISTORY))
  !
  ! initialize refinements arrays
  call init_refinements(trim(FILE_REFINE))
  !
  ! generate constraint arrays - needed???
  call init_cnstr
  !
  ! read control file
  call read_control(trim(FILE_CONTROL))
  !
  ! the following is a method to automatically setup the equation settings
  ! by silently reading the physics file and using the data there
  ! alternatively one could set this manually if desired (depends on problem)
  !
  ! read physics file quietly first to automatically setup equation settings
  call read_input(trim(FILE_PHYS))
  !
  ! setup equation settings based on physics data and common problem data
  MAXNRHS_PROB = NR_RHS_PROB !from common_prob_data
  MAXEQNH_PROB = max(1,NRHVAR) !from physics - after quietly reading physics
  MAXEQNE_PROB = max(1,NREVAR) !from physics - after quietly reading physics
  MAXEQNV_PROB = max(1,NRVVAR) !from physics - after quietly reading physics
  MAXEQNQ_PROB = max(1,NRQVAR) !from physics - after quietly reading physics
  call set_parameters(NRCOMS_PROB,MAXNRHS_PROB,  &
      MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
  !
  ! setup geometry settings
  call set_gmp_parameters( &
       NDIM_PROB, MANDIM_PROB, MAXSU_PROB, &
       MAXNP_PROB, MAXNC_PROB, MAXTR_PROB, MAXRE_PROB, &
       MAXBT_PROB, MAXHE_PROB, MAXTE_PROB, MAXPY_PROB)
  !
  ! read geometry
  call read_geometry(trim(FILE_GEOM))
  !
  ! optionally, choose partitioning algorithm for the initial mesh
  !             if using MPI parallelism
  if (NUM_PROCS > 1) call zoltan_w_set_lb(ZOLTAN_LB_DEFAULT)
  !
  ! generate mesh and read physics file
  ! keep integration flag value
  INTEGRATION_tmp = INTEGRATION
  call hp3gen(trim(FILE_PHYS))
  INTEGRATION = INTEGRATION_tmp
  !
end subroutine initialize
