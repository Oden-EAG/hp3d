!----------------------------------------------------------------------------
!> Purpose : initialize problem dependent environments, solvers, graphics and
!!           create initial mesh.
!!
!! @date Jun 15
!----------------------------------------------------------------------------
subroutine initialize
   use environment
   use common_prob_data_UW
   use data_structure3D
   use refinements
   use control
   use gmp
   use geometry_transformations
   use testvars
   use frsolmod
   use upscale
   use paraview
   use mpi_param
   use salt_data
!
   implicit none
!
!--------------------------------------------------------------------------
! Geometry Parameters
!--------------------------------------------------------------------------
   integer, parameter :: NDIM_PROB    = 3      ! Dimension
   integer, parameter :: MANDIM_PROB  = 3      ! Manifold dimension
                                               ! EXPECTED MAX NUMBER OF:
   integer, parameter :: MAXSU_PROB   = 200    ! surfaces
   integer, parameter :: MAXNP_PROB   = 1024   ! points
   integer, parameter :: MAXNC_PROB   = 2688   ! curves
   integer, parameter :: MAXTR_PROB   = 200    ! triangles
   integer, parameter :: MAXRE_PROB   = 2352   ! rectangles
   integer, parameter :: MAXBT_PROB   = 200    ! prisms
   integer, parameter :: MAXHE_PROB   = 684    ! hexas
   integer, parameter :: MAXTE_PROB   = 200    ! tetras
   integer, parameter :: MAXPY_PROB   = 200    ! pyramids
!
!--------------------------------------------------------------------------
! Equation parameters
!--------------------------------------------------------------------------
   integer, parameter :: NRCOMS_PROB  = 1     ! number of component
   integer :: MAXNRHS_PROB                    ! MAX nr rhs
   integer :: MAXEQNH_PROB                    ! MAX H1 var
   integer :: MAXEQNE_PROB                    ! MAX Hcurl var
   integer :: MAXEQNV_PROB                    ! MAX Hdiv var
   integer :: MAXEQNQ_PROB                    ! MAX L2 var
!
!..misc
   integer :: iflag,i,INTEGRATION_tmp
   character(len=1024) :: argv
!--------------------------------------------------------------------------
!
!..read history
   call open_history_file(trim(FILE_HISTORY))
!
!..intialize refinements
   call init_refinements(trim(FILE_REFINE))
!
!..generate constraint arrays
   call init_cnstr
!
!..read control file
   call read_control(trim(FILE_CONTROL))
!
   if (MAXNODS_USER .gt. 0) MAXNODS = MAXNODS_USER
   if (.not. QUIET_MODE) then
      write(*,*) 'User specified MAXNODS value:'
      write(*,9999) MAXNODS
      9999 format(' MAXNODS = ',i12)
   endif
!
!..read physics file quietly first to automatically setup equation settings
   QUIET_MODE = .true.
   if (RANK .eq. ROOT) QUIET_MODE = .false.
   call read_input(trim(FILE_PHYS))
!
!..set equation sizes
   MAXNRHS_PROB = NR_RHS_PROB !from common_prob_data
   MAXEQNH_PROB = max(0,NRHVAR) !from physics - after quietly reading physics
   MAXEQNE_PROB = max(0,NREVAR) !from physics - after quietly reading physics
   MAXEQNV_PROB = max(0,NRVVAR) !from physics - after quietly reading physics
   MAXEQNQ_PROB = max(0,NRQVAR) !from physics - after quietly reading physics
   call set_parameters(NRCOMS_PROB,MAXNRHS_PROB,                  &
            MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
!
!..set geometry sizes
   call set_gmp_parameters( &
            NDIM_PROB, MANDIM_PROB, MAXSU_PROB, &
            MAXNP_PROB, MAXNC_PROB, MAXTR_PROB, MAXRE_PROB, &
            MAXBT_PROB, MAXHE_PROB, MAXTE_PROB, MAXPY_PROB)
!
!..read geometry
   call read_geometry(trim(FILE_GEOM))
!
!..generate mesh
   INTEGRATION_tmp = INTEGRATION
   call hp3gen(trim(FILE_PHYS))
   INTEGRATION = INTEGRATION_tmp
! 
!..get command line arguments
   call get_command(argv)
!
   call load_vis(TETR_VIS, trim(FILE_VIS)//'/tetra_'//trim(VLEVEL), TETR)
   call load_vis(PRIS_VIS, trim(FILE_VIS)//'/prism_'//trim(VLEVEL), PRIS)
   call load_vis(HEXA_VIS, trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),  BRIC)
!
end subroutine initialize
