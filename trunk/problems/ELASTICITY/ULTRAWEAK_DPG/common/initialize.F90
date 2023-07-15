!----------------------------------------------------------------------------
!> @brief   Initialize problem dependent environments and create mesh
!!
!> @date    July 2023
!----------------------------------------------------------------------------
   subroutine initialize
!
      use environment
      use common_prob_data
      use data_structure3D
      use refinements
      use control
      use gmp
      use upscale
      use paraview
      use assembly,           only: NR_RHS
!
      implicit none
!
!--------------------------------------------------------------------------
!     G E O M E T R Y    S E T T I N G
!--------------------------------------------------------------------------
      integer, parameter :: NDIM_PROB    = 3    ! Dimension
      integer, parameter :: MANDIM_PROB  = 3    ! Manifold dimension
                                                ! EXPECTED MAX NUMBER OF:
      integer, parameter :: MAXSU_PROB   = 5    ! surfaces
      integer, parameter :: MAXNP_PROB   = 100  ! points
      integer, parameter :: MAXNC_PROB   = 120  ! curves
      integer, parameter :: MAXTR_PROB   = 105  ! triangles
      integer, parameter :: MAXRE_PROB   = 100  ! rectangles
      integer, parameter :: MAXBT_PROB   = 1    ! prismsm
      integer, parameter :: MAXHE_PROB   = 100  ! hexas
      integer, parameter :: MAXTE_PROB   = 100  ! tetras
      integer, parameter :: MAXPY_PROB   = 1    ! pyramids
!
!--------------------------------------------------------------------------
! E Q U A T I O N    S E T T I N G
!--------------------------------------------------------------------------
      integer, parameter :: NRCOMS_PROB  = 1     ! number of component
      integer :: MAXNRHS_PROB                    ! MAX nr rhs
      integer :: MAXEQNH_PROB                    ! MAX H1 var
      integer :: MAXEQNE_PROB                    ! MAX Hcurl var
      integer :: MAXEQNV_PROB                    ! MAX Hdiv var
      integer :: MAXEQNQ_PROB                    ! MAX L2 var
!
      integer :: iflag,i,INTEGRATION_tmp
      character(len=1024) :: argv
!--------------------------------------------------------------------------
!
!  ...initialize refinements arrays
      call init_refinements(trim(FILE_REFINE))
!
!  ...generate constraint arrays - needed???
      call init_cnstr
!
!  ...read control file
      call read_control(trim(FILE_CONTROL))
!
!  ...read physics file quietly first to automatically setup equation settings
      QUIET_MODE = .TRUE.
      call read_input(trim(FILE_PHYS))
      QUIET_MODE = .FALSE.
!
      if (MAXNODS_USER .gt. 0) MAXNODS = MAXNODS_USER
!
!  ...setup equation settings based on physics data and common problem data
      MAXNRHS_PROB = NR_RHS         !from common_prob_data
      MAXEQNH_PROB = max(1,NRHVAR)  !from physics - after quietly reading physics
      MAXEQNE_PROB = max(1,NREVAR)  !from physics - after quietly reading physics
      MAXEQNV_PROB = max(1,NRVVAR)  !from physics - after quietly reading physics
      MAXEQNQ_PROB = max(1,NRQVAR)  !from physics - after quietly reading physics
      call set_parameters(NRCOMS_PROB,MAXNRHS_PROB,  &
                          MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
!
!  ...setup geometry settings
      call set_gmp_parameters(NDIM_PROB,  MANDIM_PROB, MAXSU_PROB,             &
                              MAXNP_PROB, MAXNC_PROB,  MAXTR_PROB, MAXRE_PROB, &
                              MAXBT_PROB, MAXHE_PROB,  MAXTE_PROB, MAXPY_PROB)
!
!  ...read geometry
      call read_geometry(trim(FILE_GEOM))
!
!  ...generate mesh and read physics file
      INTEGRATION_tmp = INTEGRATION
      call hp3gen(trim(FILE_PHYS))
      INTEGRATION = INTEGRATION_tmp
!
!  ...initialize paraview data
      call load_vis(TETR_VIS, trim(FILE_VIS)//'/tetra_'//trim(VLEVEL), TETR)
      call load_vis(PRIS_VIS, trim(FILE_VIS)//'/prism_'//trim(VLEVEL), PRIS)
      call load_vis(HEXA_VIS, trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),  BRIC)
!
!  ...if using PETSc
#ifdef __PETSC_USE__
      call PetscInitialize("",iflag)
      write(*,*) 'PETSC initialized'
#endif
!
   end subroutine initialize
