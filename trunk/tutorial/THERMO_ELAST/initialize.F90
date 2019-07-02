!> Purpose : initialize problem dependent environments and create initial mesh
!! @param[in] Args Three necessary input arguments from main
subroutine initialize(Args)
  use control
  use frsolmod
  use gmp
  use data_structure3D
  use refinements

  use thermo_elast
  character(len=*), dimension(3) :: Args

  ! read control and geometry
  call read_control        (Args(1))
  call set_gmp_parameters( NDIM_PROB, MANDIM_PROB, MAXSU_PROB, &
       MAXNP_PROB, MAXNC_PROB, MAXTR_PROB, MAXRE_PROB, &
       MAXBT_PROB, MAXHE_PROB, MAXTE_PROB, MAXPY_PROB)
  ! set problem dependent parameters
  ! - input arguments
  ! (NRCOMS, MAXNRHS,     MAXEQNH, MAXEQNE, MAXEQNV, MAXEQNQ)
  call set_parameters(NRCOMS_PROB,MAXNRHS_PROB,   &
                      MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
  !
  ! output file open for the history of refinements
  call open_history_file   ('files/history')
  !     
  ! frontal solver workspace setting
  call set_frsol_workspace (1000000)
  !      
  ! generate quadrature data
  call const
  !
  ! initialize refinements array
  call init_refinements    ('../../files/ref')
  !
  ! generate constraint arrays
  call init_cnstr
  !
  ! read geometry
  call read_geometry       (Args(2))
  !     
  ! generate mesh
  call hp3gen              (Args(3))

  call set_mesh_ref_filter

end subroutine initialize

!> Purpose : finalize the program
subroutine finalize
  use data_structure3D
  ! close hisotry file
  call close_history_file
end subroutine finalize

