!-------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!-------------------------------------------------------------------------
!
module dome
  !
  use error
  use data_structure3D
  save
  !
  ! problem dependent parameters
  !-----------------------------------------------------------------------
  ! E Q U A T I O N    S E T T I N G
  integer, parameter :: NRCOMS_PROB  = 1   ! copies of variables
  integer, parameter :: NR_RHS_PROB  = 1   ! rhs
  integer, parameter :: MAXNRHS_PROB = 1   ! MAX nr rhs
  integer, parameter :: MAXEQNH_PROB = 3   ! MAX H1 var
  integer, parameter :: MAXEQNE_PROB = 1   ! MAX Hcurl var
  integer, parameter :: MAXEQNV_PROB = 1   ! MAX Hdiv var
  integer, parameter :: MAXEQNQ_PROB = 1   ! MAX L2 var
  integer, parameter :: MAXEQN_PROB = &
       max(MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB   = 3    ! Dimension
  integer, parameter :: MANDIM_PROB = 3    ! Manifold dimension
  integer, parameter :: MAXSU_PROB = 100   ! surfaces
  integer, parameter :: MAXNP_PROB = 1000  ! point
  integer, parameter :: MAXNC_PROB = 1000  ! curve
  integer, parameter :: MAXTR_PROB = 1000  ! triangle
  integer, parameter :: MAXRE_PROB = 1000  ! rectangle
  integer, parameter :: MAXBT_PROB = 1000  ! prism
  integer, parameter :: MAXHE_PROB = 1000  ! hexa
  integer, parameter :: MAXTE_PROB = 1000  ! tetra
  integer, parameter :: MAXPY_PROB = 1000  ! pyramid
  !
  ! display setting
  !-----------------------------------------------------------------------
  integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP, ICOMPLEX_DISP
  !
  !-----------------------------------------------------------------------
  ! material data
  real*8, parameter :: LAMBDA = 1.d0
  real*8, parameter :: MU     = 1.d0
  real*8, parameter :: RHO    = 1.d0
  real*8, parameter :: OMEGA  = 1.d0
  !
  !
  !  
contains
  !
  !
  !-----------------------------------------------------------------------
  subroutine disp_soldis(Nstream)
    implicit none
    integer, intent(in) :: Nstream

    write(Nstream,100)
    write(Nstream,310)
    select case (IEXACT_DISP)
    case(1); write(Nstream,101)
    case(2); write(Nstream,102)
    end select
    select case (ITANGENT_DISP)
    case(1); write(Nstream,103)
    case(2); write(Nstream,104)
    end select
    select case (ICOMPLEX_DISP)
    case(1); write(Nstream,201)
    case(2); write(Nstream,202)
    end select
    select case (ICHOOSE_DISP)
    case(1); write(Nstream,301)
    case(2); write(Nstream,302)
    case(3); write(Nstream,303)
    end select
    write(Nstream,310)

100 format('DISPLAY SETUP')
101 format('EXACT SOLUTION IS CHOSEN')
102 format('APPROX SOLUTION IS CHOSEN')
103 format('VECTOR COMPONENT IS CHOSEN')
104 format('TANGENT COMPONENT IS CHOSEN')
201 format('REAL VALUE')
202 format('IMAGINARY VALUE')
301 format('Ex IS CHOSEN')
302 format('Ey IS CHOSEN')
303 format('Ez IS CHOSEN')
304 format('SAR IS CHOSEN')
310 format('-----------------------')
  end subroutine disp_soldis
!
!
end module dome
