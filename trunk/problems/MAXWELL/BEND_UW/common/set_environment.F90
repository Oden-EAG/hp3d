!----------------------------------------------------------------------
!
!     routine name      - set_environment_maxwell
!
!----------------------------------------------------------------------
!
!     latest revision:  - Sep 2023
!
!> @brief         - define options for main file specific to the problem.
!                         These can be consulted with the -help option when running
!                         the executable. It is especially important to define
!                         the global environment variables in module/environment.
!                         The other options are problem specific.
!
!     arguments:        - none
!
!----------------------------------------------------------------------
!
subroutine set_environment_maxwell
!
   use environment
   use commonParam
   use paraview
   use parametersDPG
!
   implicit none
!
#if HP3D_USE_OPENMP
   integer :: nthreads
#endif
!
! =============================
! ======= Control file ========
! =============================
   call get_option_string( '-file_control', 'Control file', './control/control', FILE_CONTROL)
!
! =============================
! ========= GEOMETRIES ========
! =============================
   call get_option_string( '-file_geometry', 'Geometry file', './geometries/bent_square_waveguide_cart', FILE_GEOM )
!
! =============================
! ========== PHYSICS ==========
! =============================
   call get_option_string( '-file_phys', 'Physics file', './input/physics', FILE_PHYS)
!
!..Misc (refinement, history error files)
   call get_option_string( '-file_refinement', 'Refinement files location', '../../../files/ref', FILE_REFINE )
   call get_option_string( '-file_history'   , 'History file'             , './history'         , FILE_HISTORY)
   call get_option_string( '-file_err'       , 'Error file'               , './dump_err'        , FILE_ERR    )
!
! =============================
! ====== PROBLEM PARAMS =======
! =============================
!..read in problem dependent parameters
!..option label // explanation // default value // parameter
!
!..polynomial orders
   call get_option_int('-p'  , 'IP'      , 2, IP)
   call get_option_int('-dp' , 'NORD_ADD', 1, NORD_ADD    )
   call get_option_int('-npx', 'NPX'     , 5, NPX         )
   call get_option_int('-npy', 'NPY'     , 5, NPY         )
   call get_option_int('-npz', 'NPZ'     , 5, NPZ         )
!
!  ...... ICOMP, ISOL
   call get_option_int('-comp', 'ICOMP_EXACT', 1, ICOMP_EXACT  )
   call get_option_int('-isol', 'ISOL'       ,10, ISOL         )
!
   call get_option_int('-imax', 'IMAX', 3, IMAX)
   call get_option_int('-job' , 'JOB' , 0, JOB )
!
   call get_option_int('-maxnods','MAXNODS_USER', 0, MAXNODS_USER)
!
!..ALPHA (scaling coefficient in UW test norm)
   call get_option_real('-alpha'  , 'ALPHA_NORM', 1.d0, ALPHA_NORM)
!
!..MU, EPSILON
   call get_option_real('-mu'     , 'MU'     , 1.d0, MU     )
   call get_option_real('-epsilon', 'EPSILON', 1.d0, EPSILON)
!
!..Set frequency OMEGA and impedance constant GAMMA
   call get_option_real('-omega' , 'OMEGA', 2.d0*PI, OMEGA)
   call get_option_real('-gamma' , 'GAMMA', 1.0d0  , GAMMA)
!
!..Set envelope wavenumber WAVENUM_K
   call get_option_real('-k' , 'WAVENUM_K', 16.d0, WAVENUM_K)
!
!..Set bending radius RBEND
   call get_option_real('-rbend' , 'RBEND', 10.d0, RBEND)
!
!..IBCFLAG: 0 (dirichlet)
!           2 (impedance via penalty method)
!           3 (impedance via elimination)
   call get_option_int ( '-ibc', 'IBCFLAG', 0, IBCFLAG )
!
! =============================
! ========= PARAVIEW ==========
! =============================
!
!     -- Paraview Interface --
! Variables relevant to src/modules/paraview
! option label // explanation // default value // parameter
   call get_option_string('-prefix'          ,'Prefix paraview file'               ,'bend'              , PREFIX  )
   call get_option_string('-file_vis_upscale','Visualization upscale file location','../../../files/vis', FILE_VIS)
   call get_option_string('-vis_level'       ,'Visualization upscale level (0-3)'  ,'2'                 , VLEVEL  )
!
   call get_option_bool('-paraview_ho' , 'Enable higher order element output', .false., SECOND_ORDER_VIS)
   call get_option_bool('-paraview_vtu', 'Enable VTU output format'          , .false., VIS_VTU         )
!
!..I/O
   call get_option_string('-dir_output','Paraview root directory','../outputs/',OUTPUT_DIR)
   PARAVIEW_DIR = trim(OUTPUT_DIR)//'paraview/'
!
!..Paraview MISC
   call get_option_bool('-paraview_geom', 'Dump geom at every Paraview call', .true., PARAVIEW_DUMP_GEOM)
   call get_option_bool('-paraview_attr', 'Dump solution to Paraview'       , .true., PARAVIEW_DUMP_ATTR)
!
#if HP3D_USE_OPENMP
!..number of OpenMP threads
   call get_option_int( '-nthreads', 'Number of OpenMP threads', 1, nthreads)
   call omp_set_num_threads(nthreads)
#endif
!
end subroutine set_environment_maxwell
!
