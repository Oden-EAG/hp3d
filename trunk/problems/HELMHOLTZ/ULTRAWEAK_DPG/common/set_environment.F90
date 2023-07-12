!---------------------------------------------------------------------
!> @brief      Defines command line options; can be referenced with
!!             '-help' option
!!
!> @date       July 2023
!----------------------------------------------------------------------
subroutine set_environment
!    
   use environment
   use common_prob_data_UW
   use paraview
   use parametersDPG, only: NORD_ADD, MAXNORD_ADD
! 
   implicit none
!
   integer :: nthreads
!
!..Variables relevant to src/modules/environment
!                  option label     // explanation // default value // parameter
   call get_option_string  &
        ('-file-control' ,'Control file' ,'./control/control',FILE_CONTROL)
!
   call get_option_int     &
        ('-prob','0)Free space' ,PROB_FREESPACE, PROB_KIND)

   call get_option_string  &
        ('-file-geometry','Geometry file','./geometries/hexa',FILE_GEOM)
!
   call get_option_string  &
        ('-file-phys'    ,'Physics file' ,'./input/physics'  ,FILE_PHYS)
!
   call get_option_string  &
        ('-file-history' ,'History file' ,'./input/history'  ,FILE_HISTORY)
!
   call get_option_string  &
        ('-file-err'     ,'Error file'   ,'./output/errorlogs/log.txt',FILE_ERR)
!
   call get_option_string  &
        ('-file-refinement','Refinement files location','../../../files/ref',FILE_REFINE )
!
   call get_option_string  &
        ('-prefix' ,'Prefix for paraview files' ,'acoustics_' ,PREFIX)
!
   call get_option_int('-job' , 'JOB' , 0, JOB )
!
   call get_option_string  &
        ('-file-vis-upscale','Visualization upscale file','../../../files/vis',FILE_VIS)
!
   call get_option_string  &
        ('-vis-level'       ,'Visualization upscale level (0-3)','0',VLEVEL)
!        
   call get_option_string  &
        ('-dir-paraview'    ,'Paraview root directory','../output/paraview' ,PARAVIEW_DIR)
!
   call get_option_int('-maxnods','MAXNODS_USER',0 ,MAXNODS_USER)
!
   call get_option_bool    &
        ('-paraview-geom'   ,'Dump geom at every Paraview call',.TRUE.,PARAVIEW_DUMP_GEOM)
!
   call get_option_bool    &
        ('-paraview-attr'   ,'Dump solution to Paraview'       ,.TRUE.,PARAVIEW_DUMP_ATTR)
!
   call get_option_int     &
        ('-p'  ,'Uniform order initial mesh',2 ,IP)
!        
   call get_option_int     &
        ('-dp' ,'p-enrichment order for DPG',1 ,NORD_ADD)
!
   if (NORD_ADD .gt. MAXNORD_ADD) then
      write(*,*) 'set_enviroment: NORD_ADD greater then MAXNORD_ADD'
      write(*,*) 'set_enviroment: NORD_ADD, MAXNORD_ADD = ', NORD_ADD, MAXNORD_ADD
      stop 
   endif 
   call get_option_int     &
        ('-bc','1)Dirichlet, 2)Neumann, 3)Impedance' ,3,IBC_PROB)
!        
   call get_option_int     &
        ('-exact','Manufactured solution (integer: 1-5)',3,IEXACT_PROB)
!
!..acoustics specific
!..number of wavelengths
   call get_option_real( '-rnum', 'Number of wavelengths', 4.d0, RNUM)
!
!..angular frequency   
   OMEGA = 2.0d0*dacos(-1.d0)*RNUM
!   
   call get_option_real( '-alpha', 'Scaling parameter for the test norm' , 1.d0, ALPHA)
!
#if HP3D_USE_OPENMP
!..number of OpenMP threads
   call get_option_int( '-nthreads', 'Number of OpenMP threads', 1, nthreads)
   call omp_set_num_threads(nthreads)
#endif
!
end subroutine set_environment

