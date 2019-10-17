!----------------------------------------------------------------------
!
!     routine name      - set_environment_maxwell
!
!----------------------------------------------------------------------
!
!     latest revision:  - Apr 2019
!
!     purpose:          - define options for main file specific to the problem.
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
   use laserParam
   use paraview
   use parametersDPG
!
   implicit none
!
   integer :: nthreads
!
! =============================
! ======= Control file ========
! =============================
!  read in HP3D input files location
!  (if option is not present, the default value is used)
!  option label // explanation // default value // parameter
   call get_option_string( '-file_control', 'Control file', '../COMMON_FILES/control', FILE_CONTROL)
!
! =============================
! ========= GEOMETRIES ========
! =============================
!
!..Cube waveguide (set GEOM_NO = 1)
   !call get_option_string( '-file_geometry', 'Geometry file', '../GEOMETRIES/cubes/cube', FILE_GEOM )
!
!..Rectangular waveguide (set GEOM_NO = 1)
   !call get_option_string( '-file_geometry', 'Geometry file', '../GEOMETRIES/waveguide/rect', FILE_GEOM )
!
!..Fiber core - 5 hexas (set GEOM_NO = 4)
   !call get_option_string( '-file_geometry', 'Geometry file', '../GEOMETRIES/fiber/fhcor_curv', FILE_GEOM   )
!
!..Full fiber: core (5 hexas or 4 prisms+4 hexas) and cladding (4 or 8 hexas) (set GEOM_NO = 5)
   call get_option_string( '-file_geometry', 'Geometry file', '../GEOMETRIES/fiber/fiber_prism', FILE_GEOM )
!
! =============================
! ========== PHYSICS ==========
! =============================
!
!..physics files -- for signal/pump combined Maxwell + heat
   call get_option_string( '-file_phys', 'Physics file', '../UW_COUPLED/input/physLaserUW', FILE_PHYS)
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
   call get_option_int('-px' , 'ORDER_APPROX_X', 5, ORDER_APPROX_X)
   call get_option_int('-py' , 'ORDER_APPROX_Y', 5, ORDER_APPROX_Y)
   call get_option_int('-pz' , 'ORDER_APPROX_Z', 5, ORDER_APPROX_Z)
   call get_option_int('-dp' , 'NORD_ADD'    , 1, NORD_ADD    )
   call get_option_int('-npx', 'NPX'         , 0, NPX         )
   call get_option_int('-npy', 'NPY'         , 0, NPY         )
   call get_option_int('-npz', 'NPZ'         , 0, NPZ         )
!
!  ...... ICOMP, ISOL, GEOM_NO and NO_PROBLEM
   call get_option_int('-comp', 'ICOMP_EXACT', 1, ICOMP_EXACT  )
   call get_option_int('-isol', 'ISOL'       , 1, ISOL         )
   call get_option_int('-geom', 'GEOM_NO'    , 5, GEOM_NO      )
   call get_option_int('-prob', 'NO_PROBLEM' , 3, NO_PROBLEM   )
!
   call get_option_int('-imax', 'IMAX', 3, IMAX)
   call get_option_int('-job' , 'JOB' , 0, JOB )
!
   call get_option_int('-maxnods','MAXNODS_USER',0 ,MAXNODS_USER)
!
   call get_option_int('-iproduct' , 'INNER_PRODUCT', 1, INNER_PRODUCT)
!
!..MU, EPSILON, SIGMA
   call get_option_real('-mu'     , 'MU'     , 1.d0, MU     )
   call get_option_real('-epsilon', 'EPSILON', 1.d0, EPSILON)
   call get_option_comp('-sigma'  , 'SIGMA'  , ZERO, SIGMA  )
!
!..PML VARIABLES
!..exp_coeff is only used for manufactured solution tests
   call get_option_int ('-usepml' , 'USE_PML'  , 1      , USE_PML  )
   call get_option_real('-expA'   , 'EXP_COEFF', 0.2d0  , EXP_COEFF)
   call get_option_real('-pmlfrac', 'PML_FRAC' , 0.25d0 , PML_FRAC )
   call get_option_real('-zl'     , 'ZL'       , 0.40d0 , ZL       )
!
!..Set frequency OMEGA and impedance constant GAMMA
!..the propagation constant determining number of wavelengths depends on both OMEGA and GAMMA
   !call get_option_real('-omega' , 'OMEGA', 4.5d0  , OMEGA)
   !call get_option_real('-gamma' , 'GAMMA', 1.0d0  , GAMMA)
!..OMEGA, GAMMA for FIBER
   call get_option_real('-omega' , 'OMEGA', OMEGA_SIGNAL, OMEGA) ! LP01 LMA
   !call get_option_real('-omega' , 'OMEGA', 25.7d0, OMEGA) ! LP01 single-mode
   !call get_option_real('-omega' , 'OMEGA', 40.0d0, OMEGA) ! LP11 multi-mode
   !call get_option_real('-omega' , 'OMEGA', 8.1d0*PI, OMEGA)
   call get_option_real('-gamma' , 'GAMMA', 1.0d0, GAMMA)
!
!..RECTANGULAR WAVEGUIDE
!..OMEGA, GAMMA for rectangular waveguide,
!  gamma is related to propagation constant, and is the impedance constant for TE10
!..TE10 (0<x<1)
   !call get_option_real('-omega' , 'OMEGA', (sqrt(5.d0)/2.d0)*PI          , OMEGA)
   !call get_option_real('-gamma' , 'GAMMA', sqrt(1.d0-(PI**2)/(OMEGA**2)) , GAMMA)
!..TE20 (0<x<1)
   !call get_option_real('-omega' , 'OMEGA', sqrt(5.d0)*PI                        , OMEGA)
   !call get_option_real('-gamma' , 'GAMMA', sqrt(1.d0-((2.d0*PI)**2)/(OMEGA**2)) , GAMMA)
!
!..IBCFLAG: 0 (dirichlet)
!..         3 (impedance)
   call get_option_int ( '-ibc', 'IBCFLAG', 0, IBCFLAG )
!
! =============================
! ========= PARAVIEW ==========
! =============================
!
!     -- Paraview Interface --
! Variables relevant to src/modules/paraview
! option label // explanation // default value // parameter
   !call get_option_string('-prefix'          ,'Prefix paraview file'               ,'laserUW'           , PREFIX  )
   !call get_option_string('-file_vis_upscale','Visualization upscale file location','../../../files/vis', FILE_VIS)
   !call get_option_string('-vis_level'       ,'Visualization upscale level (0-3)'  ,'3'                 , VLEVEL  )
!
!..I/O
   !call get_option_string('-dir_output','Paraview root directory','../outputs/',OUTPUT_DIR)
   !PARAVIEW_DIR = trim(OUTPUT_DIR)//'paraview/'
!
!..Paraview MISC
   !call get_option_bool('-paraview_geom', 'Dump geom at every Paraview call', .FALSE., PARAVIEW_DUMP_GEOM)
   !call get_option_bool('-paraview_attr', 'Dump solution to Paraview'       , .TRUE., PARAVIEW_DUMP_ATTR)
!
!..FOR DUMPING OUT POWER
!   call get_option_int('-dumppower', 'DUMP_POWER', 0, DUMP_POWER)
!
end subroutine set_environment_maxwell
!
!
!----------------------------------------------------------------------
!
!     routine name      - set_environment_laser
!
!----------------------------------------------------------------------
!
!     latest revision:  - Mar 2019
!
!     purpose:          - 
!
!     arguments:        - none
!
!----------------------------------------------------------------------
!
subroutine set_environment_laser
!
   use environment
   use commonParam
   use laserParam
   use paraview
   use parametersDPG
!
   implicit none
!
   integer :: nthreads
!
   call set_environment_maxwell
!
   call get_option_int ('-nlflag'         , 'NONLINEAR_FLAG' , 0       , NONLINEAR_FLAG )
   call get_option_int ('-heat'           , 'HEAT_FLAG'      , 0       , HEAT_FLAG      )
   call get_option_int ('-aniso_heat'     , 'ANISO_HEAT'     , 0       , ANISO_HEAT     )
   call get_option_int ('-aniso_ref_index', 'ANISO_REF_INDEX', 0       , ANISO_REF_INDEX)
   call get_option_int ('-copump'         , 'COPUMP'         , 1       , COPUMP         )
   call get_option_real('-raman'          , 'RAMAN_GAIN'     , 1.d-3   , RAMAN_GAIN     )
   call get_option_real('-gain'           , 'ACTIVE_GAIN'    , 1.d3    , ACTIVE_GAIN    )
!
!..time stepping: number of steps, and step size
   call get_option_int ('-nsteps'   , 'NSTEPS'  , 1       , NSTEPS  )
   call get_option_real('-dt'       , 'DELTA_T' , 0.1d0   , DELTA_T )
!
!..refractive indices and beta
   call get_option_real('-ref_core' , 'REF_INDEX_CORE', 1.4512d0, REF_INDEX_CORE)
   call get_option_real('-ref_clad' , 'REF_INDEX_CLAD', 1.4500d0, REF_INDEX_CLAD)
!
!..numerical aperture and V-number of the fiber, and
!..core-cladding index of refraction average
   NA = sqrt(REF_INDEX_CORE**2-REF_INDEX_CLAD**2)
   VNUM = OMEGA*R_CORE*NA
!
!..number of time steps for Heat equation NSTEPS = T_MAX/DELTA_T
   T_MAX = NSTEPS*DELTA_T
!
!..short fiber scaling parameter
   ALPHA_Z = ZL/FIBER_LENGTH
   !ALPHA_Z = 1.0d-6
!
!..anisotropic refractive index
   CORE_NX = REF_INDEX_CORE
   CORE_NY = REF_INDEX_CORE
   CORE_NZ = REF_INDEX_CORE
!
   CLAD_NX = REF_INDEX_CLAD
   CLAD_NY = REF_INDEX_CLAD
   CLAD_NZ = REF_INDEX_CLAD
!
   if (ANISO_REF_INDEX .eq. 1) then
      CORE_NY = 1.6510d0
      CLAD_NY = 1.6500d0
   endif
!
   CORE_N = 0; CLAD_N = 0
   CORE_N(1,1) = CORE_NX; CLAD_N(1,1) = CLAD_NX
   CORE_N(2,2) = CORE_NY; CLAD_N(2,2) = CLAD_NY
   CORE_N(3,3) = CORE_NZ; CLAD_N(3,3) = CLAD_NZ
!
!..number of OpenMP threads
   call get_option_int( '-nthreads', 'Number of OpenMP threads', 1, nthreads)
   call omp_set_num_threads(nthreads)
!
end subroutine set_environment_laser
