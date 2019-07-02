!
!----------------------------------------------------------------------
!                                                                     
!     routine name      - set environment
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - define options for main file specific to the problem.
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
  use testvars
!  
  implicit none
!
!..Variables relevant to src/modules/environment
!                  option label     // explanation // default value // parameter
   call get_option_string  &
        ('-file-control' ,'Control file' ,'./control/control',FILE_CONTROL)
!   

   call get_option_int     &
        ('-prob','0)Free space, 1)Cavity' ,PROB_FREESPACE, PROB_KIND)

   select case(PROB_KIND)
   case(PROB_FREESPACE) 
      call get_option_string  &
           ('-file-geometry','Geometry file','./geometries/hexa',FILE_GEOM)
           ! ('-file-geometry','Geometry file','./geometries/trapezoid',FILE_GEOM)
   case(PROB_CAVITY) 
      call get_option_string  &
           ('-file-geometry','Geometry file','./geometries/hexa26',FILE_GEOM)
   end select         
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
        ('-prefix' ,'Prefix for paraview files' ,'cube_' ,PREFIX)
!
!
   call get_option_string  &
        ('-file-vis-upscale','Visualization upscale file','../../../files/vis',FILE_VIS)
!
   call get_option_string  &
        ('-vis-level'       ,'Visualization upscale level (0-3)','2',VLEVEL)
!        
   call get_option_string  &
        ('-dir-paraview'    ,'Paraview root directory','./output/figures/' ,PARAVIEW_DIR)
!
   call get_option_bool    &
        ('-paraview-geom'   ,'Dump geom at every Paraview call',.TRUE.,PARAVIEW_DUMP_GEOM)
!
   call get_option_bool    &
        ('-paraview-attr'   ,'Dump solution to Paraview'       ,.TRUE.,PARAVIEW_DUMP_ATTR)
!
!
   call get_option_int     &
        ('-p'  ,'Uniform order initial mesh',1 ,IP)
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
        ('-bc','1)Dirichlet, 2)Neumann, 3)Impedance' ,BC_IMPEDANCE,IBC_PROB)
!        
   call get_option_int     &
        ('-exact','Manufactured solution (integer: 1-5)',IEXACT_GAUSS,IEXACT_PROB)
!        
   call get_option_int     &
        ('-norm-test','1) Adj Graph, 2)Mathematicians'  ,ADJOINT_GRAPH,TEST_NORM)
!
!..acoustics specific
!..number of wavelengths
   call get_option_real( '-rnum', 'Number of wavelengths' , 2.d0, RNUM)
!..angular frequency   
   OMEGA = 2.0d0*dacos(-1.d0)*RNUM
!   
   call get_option_real( '-alpha', 'Scaling parameter for the test norm' , 1.d0, ALPHA)
!
!
    call get_option_int     &
        ('-solver','1) Multigrid, 2) Smoother '  ,ISOLVER_MG,ISOLVER)
   end subroutine set_environment