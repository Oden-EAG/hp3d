!
!----------------------------------------------------------------------
!
!     routine name      - set environment
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!> @brief         - define options for main file specific to the problem.
!                         These can be consulted with the -help option when running
!                         the executable. It is especially important to define
!                         the global environment variables in module/environment.
!                         The other options are problem specific.
!
!----------------------------------------------------------------------
!
subroutine set_environment
!
   use environment
   use common_prob_data
   use paraview
!
   implicit none
!
#if HP3D_USE_OPENMP
   integer :: nthreads
#endif
!
!..Variables relevant to src/modules/environment
!         option label // explanation // default value // parameter
   call get_option_string  &
        ('-file-control'   ,'Control file'    ,'./control/control'         ,FILE_CONTROL)

   call get_option_string  &
        ('-file-geometry'  ,'Geometry file'   ,'./geometries/hexa_orient0' ,FILE_GEOM)
!
   call get_option_string  &
        ('-file-phys'      ,'Physics file'    ,'./input/physics'           ,FILE_PHYS)
!
   call get_option_string  &
        ('-file-history'   ,'History file'    ,'./input/history'           ,FILE_HISTORY)
!
   call get_option_string  &
        ('-file-err'       ,'Error file'      ,'./output/errorlogs/log.txt',FILE_ERR)
!
   call get_option_string  &
        ('-file-refinement','Refinement files','../../../files/ref'        ,FILE_REFINE)
!
!..Variables relevant to this particular application
   call get_option_int     &
        ('-isol'   ,'Exact manufactured solution'         ,0,ISOL)
!
   call get_option_int     &
        ('-imax'   ,'Number of refinements for job script',3,IMAX)
!
   call get_option_int     &
        ('-p'      ,'Uniform order initial mesh'          ,3,IP)
!
   call get_option_int     &
        ('-job'    ,'Type of job submission'              ,0,JOB)
!
   call get_option_int     &
        ('-maxnods','Maximum expected number of nodes'    ,0,MAXNODS_USER)
!
! =============================
! ========= PARAVIEW ==========
! =============================
!
   call get_option_string  &
        ('-prefix'          ,'Prefix paraview file'             ,'maxw'              ,PREFIX  )
   call get_option_string  &
        ('-file_vis_upscale','Visualization upscale files'      ,'../../../files/vis',FILE_VIS)
   call get_option_string  &
        ('-vis_level'       ,'Visualization upscale level (0-3)','2'                 ,VLEVEL  )
!
   call get_option_bool    &
        ('-paraview_ho'    ,'Enable higher order element output',.false.             ,SECOND_ORDER_VIS)
   call get_option_bool    &
        ('-paraview_vtu'   ,'Enable VTU output format'          ,.false.             ,VIS_VTU)
!
   call get_option_string  &
        ('-dir_output'      ,'Paraview root directory'          ,'../outputs/'       ,OUTPUT_DIR)
   PARAVIEW_DIR = trim(OUTPUT_DIR)//'paraview/'
!
   call get_option_bool    &
        ('-paraview_geom'   ,'Dump geom at every Paraview call' ,.true.              ,PARAVIEW_DUMP_GEOM)
   call get_option_bool    &
        ('-paraview_attr'   ,'Dump solution to Paraview'        ,.true.              ,PARAVIEW_DUMP_ATTR)
!
#if HP3D_USE_OPENMP
!..number of OpenMP threads
   call get_option_int( '-nthreads', 'Number of OpenMP threads', 1, nthreads)
   call omp_set_num_threads(nthreads)
#endif
!
end subroutine set_environment

