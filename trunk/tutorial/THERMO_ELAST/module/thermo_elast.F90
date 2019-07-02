!> Purpose : module defines data for the thermoealsticity
!----------------------------------------------------------------------
module thermo_elast
  
  use error
  use data_structure3D
  save
  !
  ! problem dependent parameters
  !-------------------------------------------------------------------------
  ! E Q U A T I O N    S E T T I N G
  integer, parameter :: NRCOMS_PROB  = 1     ! number of component
  integer, parameter :: NR_RHS_PROB  = 1     ! rhs
  integer, parameter :: MAXNRHS_PROB = 1     ! MAX nr rhs
  integer, parameter :: MAXEQNH_PROB = 4     ! MAX H1 var
  integer, parameter :: MAXEQNE_PROB = 1     ! MAX Hcurl var
  integer, parameter :: MAXEQNV_PROB = 1     ! MAX Hdiv var
  integer, parameter :: MAXEQNQ_PROB = 1     ! MAX L2 var
  integer, parameter :: MAXEQN_PROB  = max(MAXEQNH_PROB,MAXEQNE_PROB,  &
                                           MAXEQNV_PROB,MAXEQNQ_PROB)
  !
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB    = 3      ! Dimension
  integer, parameter :: MANDIM_PROB  = 3      ! Manifold dimension
  integer, parameter :: MAXSU_PROB   = 100    ! surfaces
  integer, parameter :: MAXNP_PROB   = 10000  ! point
  integer, parameter :: MAXNC_PROB   = 10000  ! curve
  integer, parameter :: MAXTR_PROB   = 10000  ! triangle
  integer, parameter :: MAXRE_PROB   = 1000   ! rectangle
  integer, parameter :: MAXBT_PROB   = 1000   ! prism
  integer, parameter :: MAXHE_PROB   = 1000   ! hexa
  integer, parameter :: MAXTE_PROB   = 10000  ! tetra
  integer, parameter :: MAXPY_PROB   = 1000   ! pyramid    

  integer :: NPX, NPY, NPZ

  integer :: ICHOOSE_DISP = 1

contains

  subroutine disp_soldis(Nstream)
    implicit none
    integer, intent(in) :: Nstream

    write(Nstream,100)
    write(Nstream,310)

    select case (ICHOOSE_DISP)
    case(1); write(Nstream,201)
    case(2); write(Nstream,202)
    case(3); write(Nstream,203)
    case(4); write(Nstream,204)
    end select

    write(Nstream,310)

100 format('DISPLAY SETUP')
201 format('U X DIRECTION IS CHOSEN')
202 format('U Y DIRECTION IS CHOSEN')
203 format('U Z DIRECTION IS CHOSEN')
204 format('TEMP IS CHOSEN')
310 format('-----------------------')
  end subroutine disp_soldis



end module thermo_elast
