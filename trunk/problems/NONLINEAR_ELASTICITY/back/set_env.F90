!> Purpose : initialize problem dependent environments and create initial mesh
subroutine initialize
  use control
  use frsolmod
  use gmp
  use data_structure3D
  use refinements
  use upscale
  use common_prob_data
  use geometry_transformations, only : geometry_transformation_GMPpoints
  use environment
  use paraview
  use uhm2
  implicit none
  !-------------------------------------------------------------------------
  ! G E O M E T R Y    S E T T I N G
  integer, parameter :: NDIM_PROB    = 3     ! Dimension
  integer, parameter :: MANDIM_PROB  = 3     ! Manifold dimension

                                            ! NUMBER OF:
  integer, parameter :: MAXSU_PROB   = 1    ! surfaces
  integer, parameter :: MAXNP_PROB   = 100  ! points
  integer, parameter :: MAXNC_PROB   = 100  ! curves
  integer, parameter :: MAXTR_PROB   = 100  ! triangles
  integer, parameter :: MAXRE_PROB   = 100  ! rectangles
  integer, parameter :: MAXBT_PROB   = 100  ! prisms
  integer, parameter :: MAXHE_PROB   = 100  ! hexas
  integer, parameter :: MAXTE_PROB   = 100  ! tetras
  integer, parameter :: MAXPY_PROB   = 100  ! pyramids
  !-------------------------------------------------------------------------
  integer :: iflag
  character(len=1024) :: argv
  !-------------------------------------------------------------------------
  ! read control and geometry
  call read_control(trim(FILE_CONTROL))
  call set_gmp_parameters( &
       NDIM_PROB, MANDIM_PROB, MAXSU_PROB, &
       MAXNP_PROB, MAXNC_PROB, MAXTR_PROB, MAXRE_PROB, &
       MAXBT_PROB, MAXHE_PROB, MAXTE_PROB, MAXPY_PROB)
  !
  ! output file open for the history of refinements
  call open_history_file(trim(FILE_HISTORY))
  !
  ! internal workspace set up
  call set_frsol_workspace (100000)
  call set_x11_workspace(40000000,60000000,200000000)
  !
  ! generate quadrature data
  !call const
  !
  ! initialize refinements array
  call init_refinements(trim(FILE_REFINE))
  !
  ! vis file
  call load_vis(TETR_VIS, trim(FILE_VIS)//'/tetra_'//trim(VLEVEL), 'tetr')
  call load_vis(PRIS_VIS, trim(FILE_VIS)//'/prism_'//trim(VLEVEL), 'pris')
  call load_vis(HEXA_VIS, trim(FILE_VIS)//'/hexa_'//trim(VLEVEL), 'hexa')
  !
  ! generate constraint arrays
  call init_cnstr
  !
  ! read geometry

  call read_geometry(trim(FILE_GEOM))
  !
  ! apply transformations (if any)
  if (COORD_TRANS) then
     write(*,*) 'Applying geometrical transformation'
     call geometry_transformation_GMPpoints
!     call dumpout_GMP
  end if
  !
  if ((INPUT_FILE.eq.RECONSTRUCT_).and.(IDOMAIN_SMOOTHE.gt.0)) then
     call smooth_domain(IDOMAIN_SMOOTHE)

     call hp3gen(trim(FILE_PHYS))
     write(*,*) 'Mesh checking begin'
     call check_left_oriented(iflag)
     if (iflag.ne.1) then
        write(*,*) 'ALERT :: Mesh is not right oriented'
        call pause
     end if
     write(*,*) 'Mesh checking begin'
     ! call update_gdof
  else
     ! generate mesh
     call hp3gen(trim(FILE_PHYS))
  end if

  ! initialize UHM solver
  UHM_VERBOSE            = .FALSE.
  UHM_SOLVER_REPORT_SHOW = .FALSE.

  ! Fortran Mystery or Our code mystery...
  call get_command(argv)
  call uhm_initialize(argv)

  call uhm_option_begin
  call uhm_option_end

  call uhm_time_begin(UHM_TIME_LEVEL_FRONT)
  call uhm_direct_lu_piv_initialize(UHM_DOUBLE, NR_RHS_PROB, 256, UHM_SOLVER_PTR)

#ifdef __PETSC_USE__
  call PetscInitialize("",iflag)
  write(*,*) 'PETSC initialized'
#endif

end subroutine initialize

!> Purpose : finalize the program
subroutine finalize
  use data_structure3D
  implicit none
  call close_history_file
end subroutine finalize

subroutine set_environment
  use environment
  use common_prob_data
  use paraview
  use parametersDPG, only: NORD_ADD
  implicit none
  call get_option_string(  &
       '-file-refinement', 'Refinement lookup file location',  &
       '../../../files/ref', FILE_REFINE)
  call get_option_string(  &
       '-file-vis-upscale', 'Visualization upscale file location',  &
       '../../../files/vis', FILE_VIS)
  call get_option_string(  &
       '-file-geometry', 'Geometry file',  &
       '../files/hexa/hexa', FILE_GEOM)
  call get_option_string(  &
       '-vis-level', 'Visualization upscale level (0-3)',  &
       '1', VLEVEL)
  call get_option_int(  &
       '-p', 'order of approximation (p>0 : uniform order)',  &
       1, IP)
  call get_option_int(  &
       '-dp', 'enrichment order     (FOR DPG ONLY)',  &
       1, NORD_ADD)
  call get_option_int(  &
       '-bc', 'Boundary Conditions 1) Dirichlet, 2) Neumann, 8) Mixed',  &
       BC_DIRICHLET, IBC_PROB)
  call get_option_int(  &
       '-exact', 'Manufactured solution 1) - 5) (see common_prob_data.F90)',  &
     IEXACT_TRILINEAR, IEXACT_PROB)
  call get_option_int(  &
       '-error', '1) L2, 2) H1, 3) Energy, 4) Strain Energy',  &
       IERROR_H1, IERROR_PROB)
  call get_option_int(  &
       '-norm', '1) Adjont Graph, 2) H1, 3) H1_symm     (FOR DPG ONLY)',  &
       H1, TEST_NORM)
  ! call get_option_int( &
  !      '-solver', 'Solver: 1) frontal, 2) MUMPS, 3) UHM', &
  !      1,isolver)
  ! call get_option_int( &
  !      '-adapt-mode', '0) iso, 1) aniso', &
  !      0, imode)

  ISEL_PARAVIEW(1:10) = 0
  ISEL_PARAVIEW(5:10) = 1

end subroutine set_environment


!  ...fake subroutines for linking

subroutine display_domains


end subroutine display_domains
