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
   call get_option_string( '-file_geometry', 'Geometry file', './geometries/torus_part', FILE_GEOM )
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
!..Base amplitudes of electrid and magnetic fields E_AML, H_AMPL
   call get_option_real('-e0' , 'E_AMPL', 1.d0, E_AMPL )
   call get_option_real('-h0' , 'H_AMPL', 1.d0, H_AMPL )
!
!..Set frequency OMEGA and impedance constant GAMMA
   call get_option_real('-omega' , 'OMEGA', 2.d0*PI, OMEGA)
   call get_option_real('-gamma' , 'GAMMA', 1.0d0  , GAMMA)
!
!..Set envelope wavenumber ENVELOPEK
   call get_option_real('-k' , 'ENVELOPEK', 0.996d0*2.d0*PI, ENVELOPEK)
!
!..Set bending radius RBEND
   call get_option_real('-rbend' , 'RBEND', 25.d0, RBEND)
!
!..Set upper and lower bounds for coordinate THETA
   call get_option_real('-thup' , 'THUP', 0.5235987755982989d0, THUP)
   call get_option_real('-thlo' , 'THLO', 0.d0, THLO)
!
!..Set upper and lower bounds for coordinate R
   call get_option_real('-rup' , 'RUP', 1310.d0, RUP)
   call get_option_real('-rlo' , 'RLO', 1290.d0, RLO)
!
!..Set upper and lower bounds for coordinate X
   call get_option_real('-xup' , 'XUP', 0.d0, XUP)
   call get_option_real('-xlo' , 'XLO', 1.d0, XLO)
!
!..Set upper and lower bounds for coordinate RHO
   call get_option_real('-rhoup' , 'RHOUP', 0.d0, RHOUP)
   call get_option_real('-rholo' , 'RHOLO',10.d0, RHOLO)
!
!..Set PML proportion for upper and lower bounds of coordinate TH (≤0.5)
   call get_option_real('-pmlthup' ,'PMLTHUP', 0.d0, PMLTHUP)
   call get_option_real('-pmlthlo' ,'PMLTHLO', 0.d0, PMLTHLO)
!
!..Set PML proportion for upper and lower bounds of coordinate R (≤0.5)
   call get_option_real('-pmlrup' , 'PMLRUP', 0.d0, PMLRUP)
   call get_option_real('-pmlrlo' , 'PMLRLO', 0.d0, PMLRLO)
!
!..Set PML proportion for upper and lower bounds of coordinate X (≤0.5)
   call get_option_real('-pmlxup' , 'PMLXUP', 0.d0, PMLXUP)
   call get_option_real('-pmlxlo' , 'PMLXLO', 0.d0, PMLXLO)
!
!..Set PML proportion for upper and lower bounds of coordinate RHO (≤0.5)
   call get_option_real('-pmlrhoup','PMLRHOUP', 0.d0, PMLRHOUP)
   call get_option_real('-pmlrholo','PMLRHOLO', 0.d0, PMLRHOLO)
!
!..Set PML flag in toroidal geometry (radial direction of bent fiber's cross section)
   call get_option_int('-tordmn' , 'TOROIDAL_DMN', 0 , TOROIDAL_DMN)
!
!..Set flag for slab waveguide geometry
   call get_option_int('-slab' , 'SLAB_GUIDE', 0 , SLAB_GUIDE)
!
!..Set core's radius RCORE, cladding's radius RCLAD, coating's radius RCOAT (halfwidths in slab guide)
   call get_option_real('-rcore' , 'RCORE', 0.5d0, RCORE)
   call get_option_real('-rclad' , 'RCLAD', 5.0d0, RCLAD)
   call get_option_real('-rcoat' , 'RCOAT', 10.d0, RCOAT)
!
!..Set real part of refractive indices of core, cladding, coating: REFRCORE, REFRCLAD, REFRCOAT, REFRAIR
   call get_option_real('-ncore' , 'REFRCORE', 1.4512d0, REFRCORE)
   call get_option_real('-nclad' , 'REFRCLAD', 1.45d0, REFRCLAD)
   call get_option_real('-ncoat' , 'REFRCOAT', 1.38d0, REFRCOAT)
   call get_option_real('-nair' , 'REFRAIR', 1.00026897d0, REFRAIR)
!
!..Set COATING attenuation in ratio per unit length
   call get_option_real('-attncoat' , 'ATTNCOAT', log(10.d0)*0.003d0*25.4d-6, ATTNCOAT)
!
!..Set ELASTOOPTIC flag
   call get_option_int('-elast' , 'ELASTOOPTIC', 0 , ELASTOOPTIC)
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
   call get_option_bool('-paraview_ho' , 'Enable higher order element output', .true., SECOND_ORDER_VIS)
   call get_option_bool('-paraview_vtu', 'Enable VTU output format'          , .true., VIS_VTU         )
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

! write(*,*) 'set_environment_maxwell: OMEGA = ',OMEGA
! write(*,*) 'set_environment_maxwell: MU    = ',MU
! write(*,*) 'set_environment_maxwell: ENVELOPEK = ',ENVELOPEK
!
end subroutine set_environment_maxwell
!
