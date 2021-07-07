program main
!
      use environment  
      use paraview    , only : VLEVEL,PARAVIEW_DUMP_ATTR,PARAVIEW_DUMP_GEOM
      use PROJ        , only : ITEST,PROBLEM,ISOLVER,IORDER_PRIS,IORDER_BRIC, &
                               IORDER_TETR,IORDER_PYRA,PERC,NLEVEL,FILE_TAGS
!     
      implicit none
      integer :: idec,mdle,kref,ip,iq,ir
!     
!--------------------------------------------------------------------------------------------      
!
!     display menu    
      idec=1
!
!  ...Initialize MPI environment
      call mpi_w_init
!      
!     initialize environment
      call begin_environment
!
!     -- HP3D --
      call get_option_string('-file-control','Control file','',FILE_CONTROL)
      call get_option_string('-file-geometry','Geometry file','',FILE_GEOM)
      call get_option_string('-file-phys','Physics file', '',FILE_PHYS)
      call get_option_string('-file-err','Error file','./files/dump_error',FILE_ERR)
      call get_option_string('-file-history','History file','./files/history',FILE_HISTORY)       
      call get_option_bool(  '-quiet-mode','Quiet mode',.FALSE.,QUIET_MODE  )
      call get_option_string('-prefix'          ,'Prefix to outputfile','proj_',PREFIX)
!
!     -- Parview Interface --
      call get_option_string('-vis-level'       ,'Visualization upscale level (0-4)','2',VLEVEL)
      call get_option_bool(  '-paraview-geom','Dump geometry at every Paraview invocation', &
                             .FALSE.,PARAVIEW_DUMP_GEOM  )
      call get_option_bool(  '-paraview-attr','Dump solution to Paraview', &
                             .FALSE.,PARAVIEW_DUMP_ATTR  )
!
!     -- TEST_PROJ --
      call get_option_string('-problem','Problem : pris, hexa, tetr, pyra, hybr','',PROBLEM)
      call get_option_string('-file-tags','Tags file','./files/dump_tags',FILE_TAGS)
      call get_option_int('-testing', &
        '0 - interactive ; 1 - projections ; 2 - orientations ; 3 - constr. approx. ; 4 - refinements',&
         0,ITEST)
      call get_option_int('-solver','Solver : 1 - frontal ; 2 - MUMPS',1,ISOLVER)
      call get_option_int('-order-pris','Prism order of approx.',1,IORDER_PRIS)
      call get_option_int('-order-bric','Hexa order of approx.',1,IORDER_BRIC)
      call get_option_int('-order-tetr','Tet order of approx.',1,IORDER_TETR)
      call get_option_int('-order-pyra','Pyramid order of approx.',1,IORDER_PYRA)
      call get_option_real('-percentage','Percentage (0,1) of elements to refine',0.d0,PERC)
      call get_option_int('-ref-levels','Levels of refinements',1,NLEVEL)
! 
      call end_environment
!
!     print header
IF (.NOT. QUIET_MODE) THEN
      write(*,*)'                     '
      write(*,*)'// -- TEST PROJ -- //'
      write(*,*)'                     '
ENDIF
!
!     initialize library
      call initialize
!
!
!     -- TESTING MODE --
      select case(ITEST)
!
!     1. multiple projections
      case(1) ; call test_1 ; idec=0
!        
!     2. multiple projections & orientations
      case(2) ; call test_2 ; idec=0
!
!     3. multiple projections & orientations & constrained approximation
      case(3) ; call test_3 ; idec=0
!
!     4. random refinements
      case(4) ; call test_4 ; idec=0
      endselect 
!
!
!     -- INTERACTIVE MODE --
!
!     suppress quiet mode
      QUIET_MODE = .FALSE.
!
!     display menu in infinite loop
      do while(idec /= 0)
!      
        write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
        write(*,*) 'Quit ...................................0'
        write(*,*) '                                         '
        write(*,*) 'Geometry graphics (X11) ................1'
        write(*,*) 'HP3D graphics (X11) ....................3'
        write(*,*) 'Paraview ...............................4'
        write(*,*) 'Print Data Structure arrays ............5'
        write(*,*) 'Dumpout Data Structure ................ 6'
        write(*,*) '                                         '
        write(*,*) ' -- Refinements  --                      '
        write(*,*) 'Interactive H-refinements .............31'
!!        write(*,*) 'Interactive P-enrichments .............32'
        write(*,*) 'Uniform     H-refinements .............33'
        write(*,*) 'Uniform     P-enrichments .............34'
        write(*,*) 'Clean up    H-refinements .............35'
        write(*,*) 'Set mesh order ........................36'
        write(*,*) 'Test 4 (random refinements) ...........37'
        write(*,*) '                                         '
        write(*,*) ' -- Projection problem --                '
        write(*,*) 'Test 1 (projection) ...................50'
        write(*,*) 'Test 2 (orientations) .................51'
        write(*,*) 'Test 3 (const. approx.) ...............52'
        write(*,*) '                                         '
        write(*,*) 'Activate dof .........................300'
!!!        write(*,*) 'Basis functions --> Paraview .........301'
        write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
        read( *,*) idec
!
        select case(idec)
!                
!       GMP x11 graphics
        case(1) ; call graphg
!                
!       hp3d x11 graphics        
        case(3) ; call graphb
!                
!       Paraview graphics
        case(4) ; call paraview_driver
!
!       print data structure        
        case(5) ; call result  
!                
!       dump out
        case(6) ; call dumpout
!
!       interactive h-refinements
        case(31)
!
!         print active elements, refinement kinds
          call display_act_elem
          call display_ref_kind
!
!         select element, refinement kind
          write(*,*) 'Select: mdle,kref ='
          read( *,*) mdle,kref
          write(*,7011) mdle,kref
 7011     format(' Refining mdle,kref = ',i4,2x,i3,' ...')       
!        
!         refine element
          call refine(mdle,kref)
          call close_mesh
          call update_gdof
          call update_ddof
          call verify_orient
          call verify_neig
!
!       interactive p-enrichments
        case(32)
!                
!         print active elements
          call display_act_elem
!
!         select element, refinement kind
          write(*,*) 'Select: mdle,ip,iq,ir ='
          read( *,*) mdle,ip,iq,ir
          write(*,7012) mdle,ip,iq,ir
 7012     format(' Refining mdle,ip,iq,ir = ',i4,2x,3(i2,2x),' ...')       
!        
!         refine element
!!!          call p_refine(mdle,ip,iq,ir)
          call enforce_min_rule 
          call update_gdof
          call update_ddof
          call verify_orient
          call verify_neig
!
!       uniform global h-refinements
        case(33)
          call global_href
          call update_gdof
          call update_ddof
          call verify_orient
          call verify_neig
!
!       global p-enrichments
        case(34) ; call global_pref

!       clean up h-refinements
        case(35) ; call cleanup_href
!
!       set mesh order
        case(36)
!                
          write(*,*) 'Set mesh order'
          read( *,*) ip
          call set_mesh_order(ip)
!                
!       random h-refinements
        case(37) ; call test_4
!
!       projection test
        case(50) ; call test_1
!          
!       orientation test
        case(51) ; call test_2
!                   
!       constrained approximation test
        case(52) ; call test_3
!
!       set dof        
        case(300) ; call set_nodal_dof      
!          
!       dump shape functions to Paraview
!!!        case(301) ; call vis_shape_funct
!
        endselect
!
!     end infinite loop
      enddo
!
!     finalize library
      call finalize
!     
!
endprogram main
