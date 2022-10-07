!----------------------------------------------------------------------------
!> Purpose : define options for main file specific to the problem. These can
!!           be consulted with the -help option when running the executable.
!!           It is especially important to define the global environment
!!           variables in module/environment.
!!           The other options are problem specific.
!! @date Jun 15
!----------------------------------------------------------------------------
subroutine set_environment
  use environment
  use common_prob_data
  use paraview
  use parametersDPG, only: NORD_ADD
  use testvars
  implicit none
!
  integer :: nthreads
!
! Variables relevant to src/modules/environment
!                        option label       // explanation               // default value             // parameter
  call get_option_string('-file-control'    ,'Control file'              ,'./control/control_unknown' ,FILE_CONTROL)
  call get_option_string('-file-geometry'   ,'Geometry file'             ,'./geometries/hose'         ,FILE_GEOM   )
  call get_option_string('-file-phys'       ,'Physics file'              ,'./input/physics'           ,FILE_PHYS   )
  call get_option_string('-file-history'    ,'History file'              ,'./input/history'           ,FILE_HISTORY)
  call get_option_string('-file-err'        ,'Error file'                ,'./output/errorlogs/log.txt',FILE_ERR    )
  call get_option_string('-file-refinement' ,'Refinement files location' ,'../../files/ref'           ,FILE_REFINE )
!
! Variables relevant to src/modules/paraview
!                        option label       // explanation                        // default value     // parameter
  call get_option_string('-prefix'          ,'Prefix for paraview files'          ,'hose_'             ,PREFIX      )
  call get_option_string('-file-vis-upscale','Visualization upscale file location','../../files/vis'   ,FILE_VIS          )
  call get_option_string('-vis-level'       ,'Visualization upscale level (0-3)'  ,'2'                 ,VLEVEL            )
  call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'./output/figures/' ,PARAVIEW_DIR      )
  call get_option_bool(  '-paraview-geom'   ,'Dump geom at every Paraview call'   ,.TRUE.              ,PARAVIEW_DUMP_GEOM)
  call get_option_bool(  '-paraview-attr'   ,'Dump solution to Paraview'          ,.TRUE.              ,PARAVIEW_DUMP_ATTR)
!
! LOCAL variables relevant to the problem - see ../module/common_prob_data
!                        option label       // explanation                                   // default      // parameter
  call get_option_int(   '-p'               ,'Uniform order p of initial mesh approximation' ,2              ,IP         )
  call get_option_int(   '-bc'              ,'Bound. Cond.: 1)Dirichlet, 2)Neumann, 3)Mixed' ,BC_MIXED       ,IBC_PROB   )
  call get_option_int(   '-norm-trial'      ,'1)L2, 2)Natural, 3)Custom'                     ,IERROR_NATURAL ,IERROR_PROB)
  call get_option_int(   '-norm-test'       ,'1) Adj Graph, 2)Natural, 3)Custom'             ,MATHEMATICIANS ,TEST_NORM  )
!
!..number of OpenMP threads
   call get_option_int( '-nthreads', 'Number of OpenMP threads', 1, nthreads)
!  TODO: If OMP environment is active, call omp_set_num_threads
!   call omp_set_num_threads(nthreads)
!
! The variable dp defining enrichment in DPG problems - src/modules/parametersDPG
!                        option label     // explanation                // default  // parameter
  call get_option_int(   '-dp'            ,'p-enrichment order for DPG' ,1          ,NORD_ADD)
!
! The DPG test variables (could be specified manually if desired) - src/modules/testvars
!                        option label     // explanation                           // default value     // parameter
  call get_option_string('-file-testvars' ,'Test variables file for DPG'           ,'./input/testvars'  ,FILE_TESTVARS)
!
!
end subroutine set_environment
