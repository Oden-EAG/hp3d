!----------------------------------------------------------------------------
!> Purpose : define options for main file specific to the problem. These can
!!           be consulted with the -help option when running the executable.
!!           It is especially important to define the global environment
!!           variables in module/environment.
!!           The other options are problem specific.
!!
!!           Only for Galerkin (not DPG)
!! @date Jun 20
!----------------------------------------------------------------------------
subroutine set_environment
  use environment
  use common_prob_data
  use paraview
  implicit none
!
! Variables relevant to src/modules/environment
!                        option label     // explanation               // default value             // parameter
  call get_option_string('-file-control'    ,'Control file'              ,'../common/control/control' ,FILE_CONTROL)
  call get_option_string('-file-geometry'   ,'Geometry file'             ,'../common/geometries/hexa' ,FILE_GEOM   )
  call get_option_string('-file-phys'       ,'Physics file'              ,'./input/physics'           ,FILE_PHYS   )
  call get_option_string('-file-history'    ,'History file'              ,'./input/history'           ,FILE_HISTORY)
  call get_option_string('-file-err'        ,'Error file'                ,'./output/errorlogs/log.txt',FILE_ERR    )
  call get_option_string('-file-refinement' ,'Refinement files location' ,'../../../files/ref'        ,FILE_REFINE )
  call get_option_string('-prefix'          ,'Prefix for paraview files' ,'cube_'                     ,PREFIX      )
!
! Variables relevant to src/modules/paraview
!                        option label     // explanation                        // default value     // parameter
  call get_option_string('-file-vis-upscale','Visualization upscale file location','../../../files/vis',FILE_VIS          )
  call get_option_string('-vis-level'       ,'Visualization upscale level (0-3)'  ,'2'                 ,VLEVEL            )
  call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'./output/figures/' ,PARAVIEW_DIR      )
  call get_option_bool(  '-paraview-geom'   ,'Dump geom at every Paraview call'   ,.FALSE.              ,PARAVIEW_DUMP_GEOM)
  call get_option_bool(  '-paraview-attr'   ,'Dump solution to Paraview'          ,.TRUE.              ,PARAVIEW_DUMP_ATTR)
!
! LOCAL variables relevant to the problem - see ../module/common_prob_data
!                        option label     // explanation                                   // default       // parameter
  call get_option_int(   '-p'               ,'Uniform order p of initial mesh approximation'  ,1               ,IP         )
  call get_option_int(   '-bc'              ,'Bound. Cond.: 1)Dirichlet, 2)Neumann, 8)Mixed'  ,BC_DIRICHLET    ,IBC_PROB   )
  call get_option_int(   '-exact'           ,'Manufactured solution (integer: 1-5)'           ,IEXACT_POLYNOMIAL,IEXACT_PROB)

! MPI-RELATED
!
   call get_option_int(  '-imax'   ,'Number of refinements for job script',3,IMAX)
   call get_option_int(  '-job'    ,'Type of job submission'              ,0,JOB)
   call get_option_int(  '-maxnods','Maximum expected number of nodes'    ,0,MAXNODS_USER)
!
end subroutine set_environment