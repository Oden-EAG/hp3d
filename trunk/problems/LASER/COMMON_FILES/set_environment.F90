!----------------------------------------------------------------------
!
!     routine name      - set_environment_maxwell
!
!----------------------------------------------------------------------
!
!     latest revision:  - Jul 18
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
   use CommonParam
   use paraview
   use parametersDPG
!
   implicit none
!
! =============================
! ======= Control file ========
! =============================
!  read in HP3D input files location
!  (if option is not present, the default value is used)
!  option label // explanation // default value // parameter
   call get_option_string( '-file-control', 'Control file', '../COMMON_FILES/control', FILE_CONTROL)
!
! =============================
! ========= GEOMETRIES ========
! =============================
!
!..geometry files -- cube and rectangular waveguides
   call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/cubes/cube'             , FILE_GEOM )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/cubes/cube_waveguide8'  , FILE_GEOM )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/cubes/cube_waveguide4'  , FILE_GEOM )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/cubes/cube_waveguide16' , FILE_GEOM )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/cubes/cube_waveguide512', FILE_GEOM )
!
!..hexa geometry files -- fiber problem
!
!..only the core (set GEOM=4)
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/hexas/fhcor_curv'        , FILE_GEOM   )
!
!..full fiber: core and cladding (set GEOM=5)
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/hexas/fhcc_curv_04_thick', FILE_GEOM  )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/hexas/fhcc_curv_08_thick', FILE_GEOM   )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/hexas/fhcc_curv_16_thick', FILE_GEOM   )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/hexas/fhcc_curv_32_thick', FILE_GEOM   )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/hexas/fhcc_curv_4_thick' , FILE_GEOM   )
   !call get_option_string( '-file-geometry', 'Geometry file', '../GEOMETRIES/hexas/fhcc'              , FILE_GEOM  )
!
! =============================
! ========== PHYSICS ==========
! =============================
!
!..physics files -- for signal/pump combined Maxwell + heat
   call get_option_string( '-file-phys', 'Physics file', '../UW_COUPLED/input/physLaserUW', FILE_PHYS)
!
!
!..Misc (refinement, history error files)
   call get_option_string( '-file-refinement', 'Refinement files location', '../../../files/ref', FILE_REFINE )
   call get_option_string( '-file-history'   , 'History file'             , './history'         , FILE_HISTORY)
   call get_option_string( '-file-err'       , 'Error file'               , './dump_err'        , FILE_ERR    )
!
!
! =============================
! ====== PROBLEM PARAMS =======
! =============================
!..read in problem dependent parameters
!..option label // explanation // default value // parameter
!
!..polynomial orders
   call get_option_int('-p' , 'ORDER_APPROX', 5, ORDER_APPROX)
   call get_option_int('-dp', 'NORD_ADD'    , 1, NORD_ADD    )
   call get_option_int('-px', 'NPX'         , 0, NPX         )
   call get_option_int('-py', 'NPY'         , 0, NPY         )
   call get_option_int('-pz', 'NPZ'         , 0, NPZ         )
!
!  ...... ICOMP, ISOL, GEOM_NO and NO_PROBLEM
!  ...... GEOM_NO = 1 Unit Cube
!  ......         = 2 Prism Core
!  ......         = 3 Prism Fiber
!  ......         = 4 Hexa Core
!  ......         = 5 Hexa Fiber
!  ......         = 6 Unit Prism
   call get_option_int('-comp'     , 'ICOMP_EXACT'  , 1, ICOMP_EXACT  )
   call get_option_int('-isol'     , 'ISOL'         , 1, ISOL         )
   call get_option_int('-geom'     , 'GEOM_NO'      , 5, GEOM_NO      )
   call get_option_int('-prob'     , 'NO_PROBLEM'   , 3, NO_PROBLEM   )
   call get_option_int('-iproduct' , 'INNER_PRODUCT', 1, INNER_PRODUCT)
!
!..MU, EPSILON, SIGMA
   call get_option_real('-mu'     , 'MU'     , 1.d0, MU     )
   call get_option_real('-epsilon', 'EPSILON', 1.d0, EPSILON)
   call get_option_comp('-sigma'  , 'SIGMA'  , ZERO, SIGMA  )
!
!..PML VARIABLES
   call get_option_int ('-usepml' , 'USE_PML'  , 1      , USE_PML  )
   call get_option_real('-expA'   , 'EXP_COEFF', 0.2d0  , EXP_COEFF)
   call get_option_real('-pmlfrac', 'PML_FRAC' , 0.05d0 , PML_FRAC )
   call get_option_real('-zl'     , 'ZL'       , 0.40d0 , ZL       )
   !call get_option_real('-zl'     , 'ZL'       , 0.80d0 , ZL       )
   !call get_option_real('-zl'     , 'ZL'       , 1.6d0  , ZL       )
   !call get_option_real('-zl'     , 'ZL'       , 3.2d0  , ZL       )
   !call get_option_real('-zl'     , 'ZL'       , 4.0d0  , ZL       )
   !call get_option_real('-zl'     , 'ZL'       , 3.0d0  , ZL       )
!
!..OMEGA, GAMMA for fiber
!
!TODO stefan [check GAMMA parameter]
!
   call get_option_real( '-omega' , 'OMEGA', 2.d0*PI                     , OMEGA)
   !call get_option_real( '-omega' , 'OMEGA', 30.d0*PI                     , OMEGA)
   !call get_option_real( '-gamma' , 'GAMMA', sqrt(1.d0-(PI**2)/(OMEGA**2)), GAMMA)
!..OMEGA, GAMMA for mfd solutions
   !call get_option_real( '-omega' , 'OMEGA', 1.001d0                      , OMEGA)
   !call get_option_real( '-gamma' , 'GAMMA', 1.001d0                      , GAMMA)
!..OMEGA, GAMMA for rectangular waveguide
   !call get_option_real( '-omega' , 'OMEGA', sqrt(6.d0)*PI                , OMEGA)
   !call get_option_real( '-gamma' , 'GAMMA', sqrt(1.d0-(PI**2)/(OMEGA**2)), GAMMA)
!
   call get_option_real( '-waist' , 'BEAM_WAIST'     , 0.65d0             , BEAM_WAIST  )
   call get_option_int ( '-ibc'   , 'IBCFLAG'        , 0                  , IBCFLAG     )
   call get_option_real( '-iratio', 'INTENSITY_RATIO', sqrt(25.d0/1000.d0), INTENSITY_RATIO)
!
! =============================
! ========= PARAVIEW ==========
! =============================
!
!     -- Parview Interface --
! Variables relevant to src/modules/paraview
! option label // explanation // default value // parameter
   call get_option_string( '-prefix'         ,'Prefix paraview file'               ,'laserUW'           , PREFIX  )
   call get_option_string('-file-vis-upscale','Visualization upscale file location','../../../files/vis', FILE_VIS)
   call get_option_string('-vis-level'       ,'Visualization upscale level (0-3)'  ,'3'                 , VLEVEL  )
!
! ... Paraview -- L=04, Nonlin Problem
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL04/NL04gain000001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL04/NL04gain00001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL04/NL04gain0001/' ,PARAVIEW_DIR      )


! ... Paraview -- L=08, Nonlin Problem
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL08/NL08gain000001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL08/NL08gain00001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL08/NL08gain0001/' ,PARAVIEW_DIR      )


! ... Paraview -- L=1.6, Nonlin Problem
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL16/NL16gain000001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL16/NL16gain00001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL16/NL16gain0001/' ,PARAVIEW_DIR      )

! ... Paraview -- L=3.2, Nonlin Problem
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL32/NL32gain000001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL32/NL32gain00001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL32/NL32gain0001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL32/NL32gain00025/' ,PARAVIEW_DIR      )


! ... Paraview -- L=4, Nonlin Problem
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL4/NL4gain000001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL4/NL4gain00001/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL4/NL4gain0001/' ,PARAVIEW_DIR      )
!

! ... Paraview -- L=4, Nonlin Problem, CounterPump
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/NL4/NL4gain0002_COUNTER/' ,PARAVIEW_DIR      )


! ... Paraview -- L=3.2, Linear Problem
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/Lin32/' ,PARAVIEW_DIR      )


! ... Paraview -- L=4, Linear Problem
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/THESIS/Lin4/' ,PARAVIEW_DIR      )


! ... Paraview -- TEST FOLDERS and waveguide folder
!     the root for the paraview path is the problem directory (UW coupled)
   !call get_option_string('-dir-paraview','Paraview root directory','../UW_COUPLED/paraview/',PARAVIEW_DIR)
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/LINTESTS/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/NLTESTS_COPUMP/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/NLTESTS_COUNTER/' ,PARAVIEW_DIR      )
   !call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'../COMMON_FILES/outputs/laser/laserCoupled/waveguide_512_thesis/' ,PARAVIEW_DIR      )
!
!..Paraview MISC
   call get_option_bool('-paraview-geom', 'Dump geom at every Paraview call', .TRUE., PARAVIEW_DUMP_GEOM)
   call get_option_bool('-paraview-attr', 'Dump solution to Paraview'       , .TRUE., PARAVIEW_DUMP_ATTR)
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
!     latest revision:  - Jul 18
!
!     purpose:          - TODO
!
!     arguments:        - none
!
!----------------------------------------------------------------------
!
subroutine set_environment_laser
!
   use environment
   use CommonParam
   use LaserParam
   use paraview
   use parametersDPG
!
   implicit none
!
   call set_environment_maxwell
!
   call get_option_int ('-nlflag'   , 'NONLINEAR_FLAG', 1       , NONLINEAR_FLAG)
   call get_option_int ('-lasermode', 'LASER_MODE'    , 0       , LASER_MODE    )
   call get_option_int ('-copump'   , 'COPUMP'        , 1       , COPUMP        )
   call get_option_real('-ntheta'   , 'NTHETA'        , 1.d0    , NTHETA        )
   call get_option_real('-raman'    , 'RAMAN_GAIN'    , 1.D-4   , RAMAN_GAIN    )
   call get_option_real('-lsscale'  , 'LS_SCALE'      , 0.5D0   , LS_SCALE      )
!
!..heat diffusion coeff (should be scaled appropriately for short fiber)
   call get_option_real('-kappa'    , 'kappa'         , 1.d0    , KAPPA         )
!
!..thermo optic coefficient (should be scaled appropriately for short fiber)
   call get_option_real('-alphaheat', 'ALPHA_HEAT'    , 1.d0    , ALPHA_HEAT    )
!
!..timestep delta t, and max time (tmax)
   call get_option_real('-deltat'   , 'deltaT'        , 1.0d0   , DELTAT        )
   call get_option_real('-tmax'     , 'Tmax'          , 1.0d0   , TMAX          )
   call get_option_real('-bessord'  , 'ORDER_BESSEL'  , 1.0d0   , ORDER_BESSEL  )
!
!..refractive indices and beta
   call get_option_real('-ref_core' , 'REF_INDEX_CORE', 1.4515d0, REF_INDEX_CORE)
   call get_option_real('-ref_clad' , 'REF_INDEX_CLAD', 1.450d0 , REF_INDEX_CLAD)
!
!..Socratis checks
!call get_option_real('-ref_core' , 'REF_INDEX_CORE', 2.0d0   , REF_INDEX_CORE)
!call get_option_real('-ref_clad' , 'REF_INDEX_CLAD', 1.0d0   , REF_INDEX_CLAD)
!
!..
   call get_option_real('-beta_prop', 'BETA_PROPAGATION', 1.45075d0, BETA_PROPAGATION)
!
!..Number of time steps for Heat equation NSTEPS = TMAX/DELTAT
   NSTEPS = int(TMAX/DELTAT)
!
!..Numerical aperture and V-number
   NA = dsqrt(REF_INDEX_CORE**2-REF_INDEX_CLAD**2)
   VNUM = OMEGA*R_CORE*NA
!
!..
   E_0 = INTENSITY_RATIO*E_INC
!
end subroutine set_environment_laser
