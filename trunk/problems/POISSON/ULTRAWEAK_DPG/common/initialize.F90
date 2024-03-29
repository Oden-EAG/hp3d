!----------------------------------------------------------------------------
!> @brief   initialize problem dependent environments, solvers, graphics and
!!          create initial mesh.
!> @date May 2023
!----------------------------------------------------------------------------
subroutine initialize
   use environment
   use data_structure3D
   use refinements
   use control
   use gmp
   use common_prob_data
   use frsolmod
   use upscale
   use paraview
!
   implicit none
!
!--------------------------------------------------------------------------
!  G E O M E T R Y    S E T T I N G
   integer, parameter :: NDIM_PROB    = 3    ! Dimension
   integer, parameter :: MANDIM_PROB  = 3    ! Manifold dimension
                                             ! EXPECTED MAX NUMBER OF:
   integer, parameter :: MAXSU_PROB   = 1    ! surfaces
   integer, parameter :: MAXNP_PROB   = 32768! points
   integer, parameter :: MAXNC_PROB   = 49152! curves
   integer, parameter :: MAXTR_PROB   = 1    ! triangles
   integer, parameter :: MAXRE_PROB   = 24576! rectangles
   integer, parameter :: MAXBT_PROB   = 1    ! prisms
   integer, parameter :: MAXHE_PROB   = 4096 ! hexas
   integer, parameter :: MAXTE_PROB   = 1    ! tetras
   integer, parameter :: MAXPY_PROB   = 1    ! pyramids
!--------------------------------------------------------------------------
!  E Q U A T I O N    S E T T I N G
   integer, parameter :: NRCOMS_PROB  = 1     ! number of solution copies
   integer, parameter :: NRRHS_PROB   = 1     ! number of rhs
!--------------------------------------------------------------------------
   integer :: INTEGRATION_tmp
!--------------------------------------------------------------------------
!..output file open for the history of refinements (not for MPI)
!  call open_history_file(trim(FILE_HISTORY))
!
!..initialize refinements arrays
   call init_refinements(trim(FILE_REFINE))
!
!..generate constraint arrays
   call init_cnstr
!
!..read control file
   call read_control(trim(FILE_CONTROL))
!
!..Overwrite MAXNODS if specified by user input via argument list
   if (MAXNODS_USER .gt. 0) MAXNODS = MAXNODS_USER
   if (.not. QUIET_MODE) then
      !write(*,*) 'User specified MAXNODS value:'
      !write(*,9999) MAXNODS
9999 format(' MAXNODS = ',i12)
   endif
!
!..setup equation settings based on physics data and common problem data
   call set_parameters(NRCOMS_PROB,NRRHS_PROB)
!
!..setup geometry settings
   call set_gmp_parameters( &
       NDIM_PROB, MANDIM_PROB, MAXSU_PROB, &
       MAXNP_PROB, MAXNC_PROB, MAXTR_PROB, MAXRE_PROB, &
       MAXBT_PROB, MAXHE_PROB, MAXTE_PROB, MAXPY_PROB)
!
!..read geometry
   call read_geometry(trim(FILE_GEOM))
!
!..generate mesh and read physics file
!..keep integration flag value
   INTEGRATION_tmp = INTEGRATION
   call hp3gen(trim(FILE_PHYS))
   INTEGRATION = INTEGRATION_tmp
!
!..frontal solver: initialize workspace
   call set_frsol_workspace(1000000)
!
!..X11 graphics optional: initialize parameters
   call set_x11_workspace(40000000,60000000,200000000)
!
end subroutine initialize
