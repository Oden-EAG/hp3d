!---------------------------------------------------------------------------
!> Purpose : define all necessary problem dependent variables
!---------------------------------------------------------------------------
!
module proj
!
  use error
  use data_structure3D
  use environement
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
  integer, parameter :: MAXEQN_PROB  = &
       max(MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
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

  integer, parameter :: DISP_EXACT_PROB  = 1
  integer, parameter :: DISP_APPROX_PROB = 2
  integer :: IEXACT_DISP, ICHOOSE_DISP

  integer, parameter :: IERROR_L2 = 1
  integer, parameter :: IERROR_H1 = 2
  integer :: IERROR_PROB

  ! I N P U T   A R G U M E N T S
  character(len=128) :: &
       FILE_CONTROL, FILE_HISTORY, FILE_REFINE, FILE_VIS, FILE_GEOM, FILE_PHYS,&
       VLEVEL, PREFIX

  ! projection function polynomial
  integer :: NPX, NPY, NPZ
  !
contains
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
    select case (ICHOOSE_DISP)
    case(1); write(Nstream,301)
    case(2); write(Nstream,302)
    case(3); write(Nstream,303)
    end select

    write(Nstream,310)

100 format('DISPLAY SETUP')
101 format('EXACT SOLUTION IS CHOSEN')
102 format('APPROX SOLUTION IS CHOSEN')
301 format('u IS CHOSEN')
302 format('u_x IS CHOSEN')
303 format('u_y IS CHOSEN')
304 format('u_z IS CHOSEN')
310 format('-----------------------')

  end subroutine disp_soldis
!
!
end module proj
