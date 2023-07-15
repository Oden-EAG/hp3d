!---------------------------------------------------------------------
!> @brief      Defines command line options; can be referenced with
!!             '-help' option
!!
!> @date       July 2023
!----------------------------------------------------------------------
   subroutine set_environment
!
      use environment
      use common_prob_data
      use paraview
      use parametersDPG, only: NORD_ADD
      use testvars
!
      implicit none
!
!  ...Variables relevant to src/modules/environment
!                            option label     // explanation               // default value             // parameter
      call get_option_string('-file-control'   ,'Control file'             ,'../common/control/control' ,FILE_CONTROL)
      call get_option_string('-file-geometry'  ,'Geometry file'            ,'../common/geometries/hexa' ,FILE_GEOM   )
      call get_option_string('-file-phys'      ,'Physics file'             ,'./input/physics'           ,FILE_PHYS   )
      call get_option_string('-file-history'   ,'History file'             ,'./input/history'           ,FILE_HISTORY)
      call get_option_string('-file-err'       ,'Error file'               ,'./output/errorlogs/log.txt',FILE_ERR    )
      call get_option_string('-file-refinement','Refinement files location','../../../files/ref'        ,FILE_REFINE )
      call get_option_string('-prefix'         ,'Prefix for paraview files','elasticityUW_'             ,PREFIX      )
!
!  ...Variables relevant to src/modules/paraview
!                            option label     // explanation                        // default value     // parameter
      call get_option_string('-file-vis-upscale','Visualization upscale file location','../../../files/vis',FILE_VIS          )
      call get_option_string('-vis-level'       ,'Visualization upscale level (0-3)'  ,'1'                 ,VLEVEL            )
      call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'./output/figures/' ,PARAVIEW_DIR      )
      call get_option_bool(  '-paraview-geom'   ,'Dump geom at every Paraview call'   ,.TRUE.              ,PARAVIEW_DUMP_GEOM)
      call get_option_bool(  '-paraview-attr'   ,'Dump solution to Paraview'          ,.TRUE.              ,PARAVIEW_DUMP_ATTR)
!
!  ...LOCAL variables relevant to the problem - see ../module/common_prob_data
!                         option label // explanation                                   // default       // parameter
      call get_option_int('-p'         , 'Uniform order p of initial mesh approximation',2               ,IP         )
      call get_option_int('-bc'        , 'Bound. Cond.: 1)Dirichlet, 2)Neumann, 8)Mixed',BC_DIRICHLET    ,IBC_PROB   )
      call get_option_int('-exact'     , 'Manufactured solution (integer: 1-5)'         ,IEXACT_TRILINEAR,IEXACT_PROB)
      call get_option_int('-norm-trial', '1)L2, 2)Natural, 3)Custom'                    ,IERROR_NATURAL  ,IERROR_PROB)
      call get_option_int('-norm-test' , '1) Adj Graph, 2)Natural, 3)Custom'            ,ADJOINT_GRAPH   ,TEST_NORM  )
      call get_option_int('-alpha'     , 'DPG test norm L2 weight'                      ,1.d0            ,ALPHA      )
!
!  ...The variable dp defining enrichment in DPG problems - src/modules/parametersDPG
      call get_option_int('-dp', 'p-enrichment order for DPG', 1, NORD_ADD)
!
   end subroutine set_environment
