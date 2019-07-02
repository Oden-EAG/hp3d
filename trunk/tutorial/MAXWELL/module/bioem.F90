!> Purpose : define all necessary problem dependent variables
!! rev@Feb.2011
!----------------------------------------------------------
module bioem
  use error
  use data_structure3D
  save
  ! problem dependent parameters
  !----------------------------------------------------------
  integer, parameter :: NRCOMS_PROB=1  ! number of component
  integer, parameter :: MAXNRHS_PROB=1 ! MAX nr rhs
  integer, parameter :: MAXEQNH_PROB=1 ! MAX H1 var
  integer, parameter :: MAXEQNE_PROB=1 ! MAX Hcurl var
  integer, parameter :: MAXEQNV_PROB=1 ! MAX Hdiv var
  integer, parameter :: MAXEQNQ_PROB=1 ! MAX L2 var
  integer, parameter :: MAXEQN_PROB = &
       max(MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
  !
  integer, parameter :: NDIM_PROB   = 3    ! Dimension
  integer, parameter :: MANDIM_PROB = 3    ! Manifold dimension
  !
  integer, parameter :: MAXSU_PROB = 100   ! surfaces
  integer, parameter :: MAXNP_PROB = 10000 ! point
  integer, parameter :: MAXNC_PROB = 10000 ! curve
  integer, parameter :: MAXTR_PROB = 10000 ! triangle
  integer, parameter :: MAXRE_PROB = 10000 ! rectangle
  integer, parameter :: MAXBT_PROB = 10000 ! prism
  integer, parameter :: MAXHE_PROB = 1     ! hexa
  integer, parameter :: MAXTE_PROB = 10000 ! tetra
  integer, parameter :: MAXPY_PROB = 1     ! pyramid
  
  ! display setting
  !----------------------------------------------------------
  integer :: IEXACT_DISP, ITANGENT_DISP, ICHOOSE_DISP, ICOMPLEX_DISP

  ! problem data
  !----------------------------------------------------------
  integer, parameter :: EM_PHYSICS = 1
  integer, parameter :: HEAT_PHYSICS = 2
  integer :: PHYSICS_TYPE

  integer :: NR_RHS_PROB = 1

  ! geometry customization
  logical :: IS_CUSTOMIZE
  integer :: NR_SPHERE_PROB
  real*8  :: THICKNESS_PROB

  ! constant
  !----------------------------------------------------------
  real*8,     parameter :: PROBLEM_TOL = 10e-10
  real*8,     parameter :: PI = 3.14159265358979323846
  real*8,     parameter :: PERMITTIVITY_0 = 8.8541878176e-12
  real*8,     parameter :: PERMEABILITY_0 = 4.0e-7*PI
  real*8,     parameter :: SPEED_LIGHT = &
       1.d0/sqrt(PERMITTIVITY_0*PERMEABILITY_0)

  complex*16, parameter :: Z_0 = cmplx(0.d0, 0.d0)
  complex*16, parameter :: Z_1 = cmplx(1.d0, 0.d0)
  complex*16, parameter :: Z_i = cmplx(0.d0, 1.d0)

  ! non-dimensionalization
  !----------------------------------------------------------
  real*8,     parameter :: SCALE_X = 1.d0
  real*8,     parameter :: SCALE_E = 1.d0
  
  ! material data
  !----------------------------------------------------------
  integer, parameter :: MAX_NR_DOMAIN_PROB = 5
  integer :: NR_DOMAIN_PROB
  real*8, dimension(MAX_NR_DOMAIN_PROB) :: &
       PERMITTIVITY, PERMEABILITY, CONDUCTIVITY
  complex*16, dimension(3, MAX_NR_DOMAIN_PROB) :: &
       CURRENT_IMPRESSED 

  ! PML data
  !----------------------------------------------------------
  real*8 :: PML_RMIN, PML_RMAX

  ! frequency
  !----------------------------------------------------------
  real*8 :: OMEGA_PROB

  ! incident wave
  !----------------------------------------------------------
  integer, parameter :: PROB_PLANE_INCIDENT = 1
  integer, parameter :: PROB_PLANE_SCATTER  = 2
  integer, parameter :: PROB_DIPOLE_SCATTER = 3
  integer :: PROB_TYPE

contains
  !----------------------------------------------------------
  subroutine set_disp(Kexact, Ktangent, Kcomplex, Kchoose)
    IEXACT_DISP   = Kexact
    ITANGENT_DISP = Ktangent
    ICOMPLEX_DISP = Kcomplex
    ICHOOSE_DISP  = Kchoose
  end subroutine set_disp
  !----------------------------------------------------------
  subroutine set_prob_type(K)
    implicit none
    integer, intent(in) :: K
    PROB_TYPE = K
  end subroutine set_prob_type
  !----------------------------------------------------------
  subroutine get_prob_type(K)
    implicit none
    integer, intent(out) :: K
    K = PROB_TYPE
  end subroutine get_prob_type
  !----------------------------------------------------------
  subroutine set_customize_geometry(Nr_sph, Thick)
    implicit none
    integer, intent(in) :: Nr_sph
    real*8,  intent(in) :: Thick
    
    NR_SPHERE_PROB = Nr_sph
    THICKNESS_PROB = Thick
    
  end subroutine set_customize_geometry
  !----------------------------------------------------------
  subroutine enable_customize_geometry
    implicit none
    IS_CUSTOMIZE = .TRUE.
  end subroutine enable_customize_geometry
  !----------------------------------------------------------
  subroutine disable_customize_geometry
    implicit none
    IS_CUSTOMIZE = .FALSE.
  end subroutine disable_customize_geometry
  !----------------------------------------------------------
  subroutine set_pml_thickness( Rmin, Rmax )
    implicit none
    real*8, intent(in) :: Rmin, Rmax
    PML_RMIN = Rmin
    PML_RMAX = Rmax
  end subroutine set_pml_thickness
  !----------------------------------------------------------
  subroutine get_pml_thickness( Rmin, Rmax )
    implicit none
    real*8, intent(out) :: Rmin, Rmax
    Rmin = PML_RMIN
    Rmax = PML_RMAX
  end subroutine get_pml_thickness
  !----------------------------------------------------------
  subroutine set_angular_velocity(Omega)
    implicit none
    real*8, intent(in) :: Omega
    OMEGA_PROB = Omega
  end subroutine set_angular_velocity
  !----------------------------------------------------------
  subroutine get_angular_velocity(Omega)
    implicit none
    real*8, intent(out) :: Omega
    Omega = OMEGA_PROB
  end subroutine get_angular_velocity
  !----------------------------------------------------------
  subroutine set_nr_domain(Nrdom)
    implicit none
    integer, intent(in) :: Nrdom
    NR_DOMAIN_PROB = Nrdom
  end subroutine set_nr_domain
  !----------------------------------------------------------
  subroutine set_em_current(I, Z)
    implicit none
    integer,                  intent(in) :: I  ! domain number
    complex*16, dimension(3), intent(in) :: Z  ! current
    CURRENT_IMPRESSED(1:3,I) = Z(1:3)
  end subroutine set_em_current
  !----------------------------------------------------------
  subroutine get_em_current(I, Z)
    implicit none
    integer,                  intent(in)  :: I  ! domain number
    complex*16, dimension(3), intent(out) :: Z  ! current
    Z(1:3) = CURRENT_IMPRESSED(1:3, I)
  end subroutine get_em_current
  !----------------------------------------------------------
  subroutine set_em_material(I, E, U, S)
    implicit none
    integer, intent(in) :: I        ! domain number
    real*8,  intent(in) :: E, U, S  ! material constant
    PERMITTIVITY(I) = E
    PERMEABILITY(I) = U
    CONDUCTIVITY(I) = S
  end subroutine set_em_material
  !----------------------------------------------------------
  subroutine get_em_material(I, E, U, S)
    implicit none
    integer, intent(in)  :: I        ! domain number
    real*8,  intent(out) :: E, U, S  ! material constant
    E = PERMITTIVITY(I) 
    U = PERMEABILITY(I)
    S = CONDUCTIVITY(I)
  end subroutine get_em_material
  !----------------------------------------------------------
  subroutine solve_problem
    implicit none
    call solve1(NR_RHS_PROB)
  end subroutine solve_problem
  !----------------------------------------------------------
  subroutine set_physics_type(Itype)
    implicit none
    integer, intent(in) :: Itype
    !
    select case (Itype)
    case (EM_PHYSICS, HEAT_PHYSICS)
       PHYSICS_TYPE = Itype
    case default
       call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
    end select
    !
  end subroutine set_physics_type
  !----------------------------------------------------------
  subroutine get_physics_type(Itype)
    implicit none
    integer, intent(out) :: Itype
    Itype = PHYSICS_TYPE
  end subroutine get_physics_type
  !----------------------------------------------------------
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
  !----------------------------------------------------------
  subroutine disp_bioem(Nstream)
    implicit none
    integer, intent(in) :: Nstream
    write(Nstream,*) ' '
    write(Nstream,1000)
    write(Nstream,2000)
    write(Nstream,3000) PHYSICS_TYPE
    write(Nstream,3001) NR_RHS_PROB
    write(Nstream,3002) IS_CUSTOMIZE, NR_SPHERE_PROB, THICKNESS_PROB
    write(Nstream,3003) NR_DOMAIN_PROB
    write(Nstream,3004) PERMITTIVITY(1:NR_DOMAIN_PROB)
    write(Nstream,3005) PERMEABILITY(1:NR_DOMAIN_PROB)
    write(Nstream,3006) CONDUCTIVITY(1:NR_DOMAIN_PROB)
    write(Nstream,4000) PML_RMIN, PML_RMAX
    write(Nstream,5000) OMEGA_PROB
    write(Nstream,6000) PROB_TYPE
    write(Nstream,6001)
    write(Nstream,6002)
    write(Nstream,2000)
1000 format(' MODULE BIOEM ')
2000 format('------------------------------')
3000 format(' PROBLEM TYPE 1) EM 2) HEAT   ', i6 )
3001 format(' NR_RHS                       ', i6 )
3002 format(' CUSTOMIZE, SPHERE, THICKNESS ', l4, i3, e12.5 ) 
3003 format(' NR_DOMAIN                    ', i6 )
3004 format(' PERMITTIVITY                 ', 5(e12.5,2x) )
3005 format(' PERMEABILITY                 ', 5(e12.5,2x) )
3006 format(' CONDUCTIVITY                 ', 5(e12.5,2x) )
4000 format(' PML REGION                   ', 2(e12.5,2x) )
5000 format(' OMEGA_PROB                   ', e12.5 )
6000 format(' INCIDENT WAVE                ', i6 )
6001 format(' 1) PLANE ' )
6002 format(' 2) PLANE SCATTER ' )
  end subroutine disp_bioem

end module bioem
