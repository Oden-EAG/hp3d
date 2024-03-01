!
!----------------------------------------------------------------------
!
!     routine name      - set environment
!
!----------------------------------------------------------------------
!
!     latest revision:  - Setp 2018
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
   use parametersDPG, only: NORD_ADD, MAXNORD_ADD
!
   implicit none
!
#if HP3D_USE_OPENMP
   integer :: nthreads
#endif
!
!..Variables relevant to src/modules/environment
!                  option label     // explanation // default value // parameter
   call get_option_string  &
        ('-file-control' ,'Control file' ,'./control/control',FILE_CONTROL)

   call get_option_string  &
        ('-file-geometry','Geometry file','./geometries/hexa_orient0',FILE_GEOM)
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
   call get_option_int     &
        ('-p'  ,'Uniform order initial mesh',3 ,IP)
!
   call get_option_int     &
        ('-dp' ,'p-enrichment order for DPG',1 ,NORD_ADD)
!
   if (NORD_ADD .gt. MAXNORD_ADD) then
      write(*,*) 'set_enviroment: NORD_ADD greater then MAXNORD_ADD'
      write(*,*) 'set_enviroment: NORD_ADD, MAXNORD_ADD = ', NORD_ADD, MAXNORD_ADD
      stop
   endif
!
#if HP3D_USE_OPENMP
!..number of OpenMP threads
   call get_option_int( '-nthreads', 'Number of OpenMP threads', 1, nthreads)
   call omp_set_num_threads(nthreads)
#endif
!
   IBC_PROB = BC_DIRICHLET
   !IP = 3
   !NORD_ADD = 1
!
end subroutine set_environment

