!----------------------------------------------------------------------------
!> Purpose : initialize problem dependent environments, solvers, graphics and
!!           create initial mesh.
!! @date Jun 15
!----------------------------------------------------------------------------
subroutine initialize
  use environment
  use common_prob_data
  use data_structure3D
  use refinements
  use control
  use gmp
  use geometry_transformations
  use testvars
  use frsolmod
  use upscale
  use paraview
  implicit none
  !--------------------------------------------------------------------------
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB    = 3    ! Dimension
  integer, parameter :: MANDIM_PROB  = 3    ! Manifold dimension
                                            ! EXPECTED MAX NUMBER OF:
  integer, parameter :: MAXSU_PROB   = 5    ! surfaces
  integer, parameter :: MAXNP_PROB   = 100  ! points
  integer, parameter :: MAXNC_PROB   = 120  ! curves
  integer, parameter :: MAXTR_PROB   = 105  ! triangles
  integer, parameter :: MAXRE_PROB   = 100  ! rectangles
  integer, parameter :: MAXBT_PROB   = 1    ! prisms
  integer, parameter :: MAXHE_PROB   = 100  ! hexas
  integer, parameter :: MAXTE_PROB   = 100  ! tetras
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
  ! output file open for the history of refinements
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
  QUIET_MODE = .TRUE.
  call read_input(trim(FILE_PHYS))
  QUIET_MODE = .TRUE.
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
  ! apply transformations (if any)
  if (COORD_TRANS) then
     write(*,*) 'Applying geometrical transformation'
     call geometry_transformation_GMPpoints
!     call dumpout_GMP
  end if
  !
  ! generate mesh and read physics file
  if ((INPUT_FILE.eq.RECONSTRUCT_).and.(IDOMAIN_SMOOTHE.gt.0)) then
     call smooth_domain(IDOMAIN_SMOOTHE)
     ! generate mesh and read physics file
     call hp3gen(trim(FILE_PHYS))
     write(*,*) 'Mesh checking begin'
     call check_left_oriented(iflag)
     if (iflag.ne.1) then
        write(*,*) 'ALERT :: Mesh is not right oriented'
        call pause
     end if
     write(*,*) 'Mesh checking begin'
     ! call update_gdof
  elseif (INPUT_FILE.eq.DUMPIN_) then
!     write(*,*)'dumping in physics...'
!     call dumpin_physics_from_file('files/dumpPHYS' )
!     write(*,*) 'dumping in GMP'
!     call dumpin_GMP('files/dumpGMP')
     write(*,*)'dumping in hp3d...'
     call dumpin_hp3d('files/dumpc3Dhp')
  else
     ! generate mesh and read physics file
     ! keep integration flag value
     INTEGRATION_tmp = INTEGRATION
     call hp3gen(trim(FILE_PHYS))
     INTEGRATION = INTEGRATION_tmp
  end if
  !
  ! read test variables
  call read_testvars(trim(FILE_TESTVARS))
  !
  !SOLVERS
  ! frontal solver: initialize workspace
  call set_frsol_workspace (100000)
  !
  ! UHM solver: initialize
  UHM_VERBOSE            = .FALSE.
  UHM_SOLVER_REPORT_SHOW = .FALSE.
  ! Fortran Mystery or Our code mystery...
  call get_command(argv)
  call uhm_initialize(argv)
  !
  call uhm_option_begin
  call uhm_option_end
  !
  call uhm_time_begin(UHM_TIME_LEVEL_FRONT)
  call uhm_direct_lu_piv_initialize(UHM_DOUBLE, NR_RHS_PROB, 256, UHM_SOLVER_PTR)
  !
  !GRAPHICS
  ! X11 graphics optional: initialize parameters
  call set_x11_workspace(40000000,60000000,200000000)
  !
  ! vis file - for vtk file for paraview
  call load_vis(TETR_VIS, trim(FILE_VIS)//'/tetra_'//trim(VLEVEL), 'tetr')
  call load_vis(PRIS_VIS, trim(FILE_VIS)//'/prism_'//trim(VLEVEL), 'pris')
  call load_vis(HEXA_VIS, trim(FILE_VIS)//'/hexa_'//trim(VLEVEL), 'hexa')
  !
  !
end subroutine initialize
