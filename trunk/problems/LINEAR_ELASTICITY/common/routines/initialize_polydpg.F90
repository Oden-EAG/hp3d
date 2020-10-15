!----------------------------------------------------------------------------
!> Purpose : initialize problem dependent environments, solvers, graphics and
!!           create initial mesh.
!! @date Jun 15
!----------------------------------------------------------------------------
subroutine initialize
  use environment
  use common_prob_data
  use data_structure3D_poly
  use control
  use geometry_polydpg
  use testvars
  use frsolmod
  use upscale
  use paraview
  implicit none
  !--------------------------------------------------------------------------
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB    = 3     ! Dimension
  integer, parameter :: MANDIM_PROB  = 3     ! Manifold dimension
                                             ! EXPECTED MAX NUMBER OF:
  integer, parameter :: MAXNV_PROB   = 14000 ! vertices
  integer, parameter :: MAXNE_PROB   = 40000 ! edges
  integer, parameter :: MAXNF_PROB   = 24000 ! faces
  integer, parameter :: MAXNC_PROB   = 9500  ! cells
  integer, parameter :: MAXF_NE_PROB   = MAX_NRVF !15  ! edges per face
  integer, parameter :: MAXC_NF_PROB   = MAX_NRFC !50 ! faces per cell
  !--------------------------------------------------------------------------
  ! E Q U A T I O N    S E T T I N G
  integer, parameter :: NRCOMS_PROB  = 1     ! number of component
  integer :: MAXNRHS_PROB                    ! MAX nr rhs
  integer :: MAXEQNH_PROB                    ! MAX H1 var
  integer :: MAXEQNE_PROB                    ! MAX Hcurl var
  integer :: MAXEQNV_PROB                    ! MAX Hdiv var
  integer :: MAXEQNU_PROB                    ! MAX Discon FaceH1 var
  integer :: MAXEQNF_PROB                    ! MAX FaceL2 var
  integer :: MAXEQNQ_PROB                    ! MAX L2 var
  !--------------------------------------------------------------------------
  integer :: iflag,i,INTEGRATION_tmp
  character(len=1024) :: argv
  !--------------------------------------------------------------------------
  !
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
  call read_input_poly(trim(FILE_PHYS))
  QUIET_MODE = .FALSE.
  ! setup equation settings based on physics data and common problem data
  MAXNRHS_PROB = NR_RHS_PROB !from common_prob_data
  MAXEQNH_PROB = max(1,NRHVAR) !from physics - after quietly reading physics
  MAXEQNE_PROB = max(1,NREVAR) !from physics - after quietly reading physics
  MAXEQNV_PROB = max(1,NRVVAR) !from physics - after quietly reading physics
  MAXEQNU_PROB = max(1,NRUVAR)
  MAXEQNF_PROB = max(1,NRFVAR) !from physics - after quietly reading physics
  MAXEQNQ_PROB = max(1,NRQVAR) !from physics - after quietly reading physics
  call set_parameters_poly(NRCOMS_PROB,MAXNRHS_PROB,  &
      MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNU_PROB,MAXEQNF_PROB,MAXEQNQ_PROB)
  !
  ! setup geometry settings
  call set_geom_poly_parameters(NDIM_PROB    , &
                                      MANDIM_PROB    , &
                                      MAXNV_PROB     , &
                                      MAXNE_PROB     , &
                                      MAXNF_PROB     , &
                                      MAXNC_PROB     , &
                                      MAXF_NE_PROB  , &
                                      MAXC_NF_PROB)
  !
  ! read geometry
  call input_polydpg(trim(FILE_GEOM))
  !
  ! 
  ! generate mesh and read physics file
  ! keep integration flag value
  INTEGRATION_tmp = INTEGRATION
  call hp3gen_poly(trim(FILE_PHYS),IP)
  INTEGRATION = INTEGRATION_tmp

  !
  ! read test variables
  call read_testvars(trim(FILE_TESTVARS))
  !
  !SOLVERS
  ! ! frontal solver: initialize workspace
  ! call set_frsol_workspace (100000)
  ! !
  ! !GRAPHICS
  ! ! X11 graphics optional: initialize parameters
  ! call set_x11_workspace(40000000,60000000,200000000)
  !
  ! vis file - for vtk file for paraview
  call load_vis(TETR_VIS, trim(FILE_VIS)//'/tetra_'//trim(VLEVEL), 'tetr')
  call load_vis(PRIS_VIS, trim(FILE_VIS)//'/prism_'//trim(VLEVEL), 'pris')
  call load_vis(HEXA_VIS, trim(FILE_VIS)//'/hexa_'//trim(VLEVEL), 'hexa')

end subroutine initialize
