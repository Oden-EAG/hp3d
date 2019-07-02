!> Purpose : initialize problem dependent environments and create initial mesh
subroutine initialize
  use control
  use frsolmod
  use gmp
  use data_structure3D
  use refinements
  use upscale

  use environment
  use lapl

  ! read control and geometry
  call read_control(trim(FILE_CONTROL))
  call set_gmp_parameters( &
       NDIM_PROB, MANDIM_PROB, MAXSU_PROB, &
       MAXNP_PROB, MAXNC_PROB, MAXTR_PROB, MAXRE_PROB, &
       MAXBT_PROB, MAXHE_PROB, MAXTE_PROB, MAXPY_PROB)
  ! set problem dependent parameters
  call set_parameters(NRCOMS_PROB,MAXNRHS_PROB,   &
                      MAXEQNH_PROB,MAXEQNE_PROB,MAXEQNV_PROB,MAXEQNQ_PROB)
  !
  ! output file open for the history of refinements
  call open_history_file(trim(FILE_HISTORY))
  !     
  ! internal workspace set up
  call set_frsol_workspace (10000)
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
  !
  ! generate constraint arrays
  call init_cnstr
  !
  ! read geometry
  call read_geometry(trim(FILE_GEOM))
  !     
  ! generate mesh
  call hp3gen(trim(FILE_PHYS))

end subroutine initialize

!> Purpose : finalize the program
subroutine finalize
  use data_structure3D
  ! close hisotry file
  call close_history_file
end subroutine finalize

