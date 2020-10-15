subroutine set_local_environment
!
  use data_structure3D
  use environment
!
  implicit none
!-------------------------------------------------------------------------
!   E Q U A T I O N    S E T T I N G
  integer, parameter :: NRCOMS_PROB  = 1     ! number of component
  integer, parameter :: MAXNRHS_PROB = 1     ! MAX nr rhs
  integer, parameter :: MAXEQNH_PROB = 3     ! MAX H1 var
  integer, parameter :: MAXEQNE_PROB = 1     ! MAX Hcurl var
  integer, parameter :: MAXEQNV_PROB = 3     ! MAX Hdiv var
  integer, parameter :: MAXEQNQ_PROB = 12    ! MAX L2 var
! integer, parameter :: MAXEQN_PROB = &
!        max(MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
  !-------------------------------------------------------------------------

  call get_option_string(  &
       '-file-control', 'Control file',  &
       '../files/control/control_DPG_UWEAK', FILE_CONTROL)
  call get_option_string(  &
       '-file-err', 'Error file',  &
       '../files/errorlogs/DPG_UWEAK/errorlog.txt', FILE_ERR)
  call get_option_string(  &
       '-file-history', 'History file',  &
       './history', FILE_HISTORY)
  call get_option_string(  &
       '-file-phys', 'Physics file',  &
       './physics', FILE_PHYS)

  ! set problem dependent parameters
  call set_parameters(NRCOMS_PROB,MAXNRHS_PROB,  &
        MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)

end subroutine set_local_environment
