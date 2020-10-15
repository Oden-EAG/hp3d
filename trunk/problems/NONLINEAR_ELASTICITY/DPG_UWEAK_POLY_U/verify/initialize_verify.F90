subroutine initialize_verify
  use environment
  use common_prob_data
  use data_structure3D_poly
  use control
  use geometry_polydpg
  implicit none
  !--------------------------------------------------------------------------
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB    = 3     ! Dimension
  integer, parameter :: MANDIM_PROB  = 3     ! Manifold dimension
                                             ! EXPECTED MAX NUMBER OF:
  integer, parameter :: MAXNV_PROB   = 1000 ! vertices
  integer, parameter :: MAXNE_PROB   = 1000 ! edges
  integer, parameter :: MAXNF_PROB   = 1000 ! faces
  integer, parameter :: MAXNC_PROB   = 512  ! cells
  integer, parameter :: MAXF_NE_PROB   = 20  ! edges per face
  integer, parameter :: MAXC_NF_PROB   = 40  ! faces per cell
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
  MAXEQNU_PROB = 1
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
  call input_polydpg(trim('input/cube'))

  write(*,*) &
   'initialize_verify: default CONTROL, PHYSICS and CUBE GEOMETRY files have been read'
  !
  

  ! generate mesh and read physics file
  ! keep integration flag value
  INTEGRATION_tmp = INTEGRATION
  call hp3gen_poly(trim(FILE_PHYS),IP)
  INTEGRATION = INTEGRATION_tmp
  !

end subroutine initialize_verify