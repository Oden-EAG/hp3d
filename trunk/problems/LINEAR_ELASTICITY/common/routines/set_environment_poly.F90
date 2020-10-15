!----------------------------------------------------------------------------
!> Purpose : define options for main file specific to the problem. These can
!!           be consulted with the -help option when running the executable.
!!           It is especially important to define the global environment
!!           variables in module/environment.
!!           The other options are problem specific.
!! @date Jun 15
!----------------------------------------------------------------------------
subroutine set_environment_poly
  use environment
  use common_prob_data
  use paraview
  use parametersDPG, only: NORD_ADD, VARIABLE_DP
  ! use testvars
  implicit none
!
! Variables relevant to src/modules/environment
!                        option label     // explanation               // default value             // parameter
  call get_option_string('-file-control'    ,'Control file'              ,'../common/control/control' ,FILE_CONTROL)
  call get_option_string('-file-geometry'   ,'Geometry file'             ,'./input/cube_21' ,FILE_GEOM   )
  call get_option_string('-file-phys'       ,'Physics file'              ,'./input/physics'           ,FILE_PHYS   )
  ! call get_option_string('-file-history'    ,'History file'              ,'./input/history'           ,FILE_HISTORY)
  call get_option_string('-file-err'        ,'Error file'                ,'./output/errorlogs/log.txt',FILE_ERR    )
  ! call get_option_string('-file-refinement' ,'Refinement files location' ,'../../../files/ref'        ,FILE_REFINE )
  call get_option_string('-prefix'          ,'Prefix for paraview files' ,'cube21_'                     ,PREFIX      )
!
! Variables relevant to src/modules/paraview
!                        option label     // explanation                        // default value     // parameter
  call get_option_string('-file-vis-upscale','Visualization upscale file location','../../../files/vis',FILE_VIS          )
  call get_option_string('-vis-level'       ,'Visualization upscale level (0-3)'  ,'0'                 ,VLEVEL            )
  call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'./output/figures/' ,PARAVIEW_DIR      )
  call get_option_bool(  '-paraview-geom'   ,'Dump geom at every Paraview call'   ,.TRUE.              ,PARAVIEW_DUMP_GEOM)
  call get_option_bool(  '-paraview-attr'   ,'Dump solution to Paraview'          ,.TRUE.              ,PARAVIEW_DUMP_ATTR)
!
! LOCAL variables relevant to the problem - see ../module/common_prob_data
!                        option label     // explanation                                   // default       // parameter
  call get_option_int(   '-p'               ,'Uniform order p of initial mesh approximation'  ,4               ,IP         )
  call get_option_int(   '-bc'              ,'Bound. Cond.: 1)Dirichlet, 2)Neumann, 8)Mixed'  ,BC_DIRICHLET    ,IBC_PROB   )
  call get_option_int(   '-exact'           ,'Manufactured solution (integer: 1-5)'           ,2,IEXACT_PROB)
  ! call get_option_int(   '-error-attribute' ,'1)Displacement, 2)Stress, 3)Combined, 4)Custom' ,DISPLACEMENT    ,IERROR_ATTR)
  call get_option_int(   '-norm-trial'      ,'1)L2, 2)Natural, 3)Custom'                      ,IERROR_NATURAL  ,IERROR_PROB)
  call get_option_int(   '-norm-test'       ,'1) Adj Graph, 2)Natural, 3)Custom'              ,2  ,TEST_NORM  )
!Other typical LOCAL options - ignored for this problem
!  call get_option_int(   '-solver'          ,'Solver: 1)frontal, 2)MUMPS, 3)UHM'            ,1               ,isolver    )
!  call get_option_int(   '-adapt-mode'      ,'0)iso, 1)aniso'                               ,0               ,imode      )
!
! The variable dp defining enrichment in DPG problems - src/modules/parametersDPG
!                        option label // explanation                             // default // parameter
  call get_option_int(   '-dp'        ,'p-enrichment order for DPG'               ,4         ,NORD_ADD)
  call get_option_int(   '-var_dp'    ,'p-enrichment order varies over elements?' ,0         ,VARIABLE_DP)
!
! The DPG test variables (could be specified manually if desired) - src/modules/testvars
!                        option label     // explanation                        // default value     // parameter
  ! call get_option_string('-file-testvars','Test variables file for DPG'           ,'./input/testvars'  ,FILE_TESTVARS)
!
!
end subroutine set_environment_poly
