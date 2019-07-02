!---------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!---------------------------------------------------------------------------
module em
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
  integer, parameter :: MAXEQNH_PROB = 1     ! MAX H1 var
  integer, parameter :: MAXEQNE_PROB = 1     ! MAX Hcurl var
  integer, parameter :: MAXEQNV_PROB = 1     ! MAX Hdiv var
  integer, parameter :: MAXEQNQ_PROB = 1     ! MAX L2 var
  integer, parameter :: MAXEQN_PROB  = max(MAXEQNH_PROB,MAXEQNE_PROB,  &
                                           MAXEQNV_PROB,MAXEQNQ_PROB)
  !
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB    = 3     ! Dimension
  integer, parameter :: MANDIM_PROB  = 3     ! Manifold dimension
  integer, parameter :: MAXSU_PROB   = 100   ! surfaces
  integer, parameter :: MAXNP_PROB   = 1000  ! point
  integer, parameter :: MAXNC_PROB   = 1000  ! curve
  integer, parameter :: MAXTR_PROB   = 1000  ! triangle
  integer, parameter :: MAXRE_PROB   = 1000  ! rectangle
  integer, parameter :: MAXBT_PROB   = 1000  ! prism
  integer, parameter :: MAXHE_PROB   = 1000  ! hexa
  integer, parameter :: MAXTE_PROB   = 1000  ! tetra
  integer, parameter :: MAXPY_PROB   = 1000  ! pyramid
  !
  ! constants
  !-------------------------------------------------------------------------
  complex*16, parameter :: Z_0 = cmplx(0.d0, 0.d0)
  complex*16, parameter :: Z_1 = cmplx(1.d0, 0.d0)
  complex*16, parameter :: Z_i = cmplx(0.d0, 1.d0)
  !
  ! material data
  !-------------------------------------------------------------------------
  real*8, parameter :: OMEGA = 1.d0          ! frequency
  real*8, parameter :: EPSIL = 1.d0          ! permittivity
  real*8, parameter :: SIGMA = 1.d0          ! conductivity 
  real*8, parameter :: MU    = 1.d0          ! permeability
  !
  ! display setting
  !-------------------------------------------------------------------------
  integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP, ICOMPLEX_DISP
  ! 
  !
  !
contains
  !
  !
  !-------------------------------------------------------------------------
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
!!!  !----------------------------------------------------------
!!!  subroutine disp_bioem(Nstream)
!!!    implicit none
!!!    integer, intent(in) :: Nstream
!!!    write(Nstream,*) ' '
!!!    write(Nstream,1000)
!!!    write(Nstream,2000)
!!!    write(Nstream,3000) PHYSICS_TYPE
!!!    write(Nstream,3001) NR_RHS_PROB
!!!    write(Nstream,3002) IS_CUSTOMIZE, NR_SPHERE_PROB, THICKNESS_PROB
!!!    write(Nstream,3003) NR_DOMAIN_PROB
!!!    write(Nstream,3004) PERMITTIVITY(1:NR_DOMAIN_PROB)
!!!    write(Nstream,3005) PERMEABILITY(1:NR_DOMAIN_PROB)
!!!    write(Nstream,3006) CONDUCTIVITY(1:NR_DOMAIN_PROB)
!!!    write(Nstream,4000) PML_RMIN, PML_RMAX
!!!    write(Nstream,5000) OMEGA_PROB
!!!    write(Nstream,6000) PROB_TYPE
!!!    write(Nstream,6001)
!!!    write(Nstream,6002)
!!!    write(Nstream,2000)
!!!1000 format(' MODULE BIOEM ')
!!!2000 format('------------------------------')
!!!3000 format(' PROBLEM TYPE 1) EM 2) HEAT   ', i6 )
!!!3001 format(' NR_RHS                       ', i6 )
!!!3002 format(' CUSTOMIZE, SPHERE, THICKNESS ', l4, i3, e12.5 ) 
!!!3003 format(' NR_DOMAIN                    ', i6 )
!!!3004 format(' PERMITTIVITY                 ', 5(e12.5,2x) )
!!!3005 format(' PERMEABILITY                 ', 5(e12.5,2x) )
!!!3006 format(' CONDUCTIVITY                 ', 5(e12.5,2x) )
!!!4000 format(' PML REGION                   ', 2(e12.5,2x) )
!!!5000 format(' OMEGA_PROB                   ', e12.5 )
!!!6000 format(' INCIDENT WAVE                ', i6 )
!!!6001 format(' 1) PLANE ' )
!!!6002 format(' 2) PLANE SCATTER ' )
!!!  end subroutine disp_bioem

end module em
