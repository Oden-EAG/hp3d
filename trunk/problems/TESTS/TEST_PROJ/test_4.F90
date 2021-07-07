subroutine test_4
!       
        use environment      , only : FILE_PHYS
        use data_structure3D , only : deallocds
        use PROJ             , only : ITEST,PERC,NLEVEL
        use paraview         , only : VLEVEL,PARAVIEW_DUMP_GEOM, &
                                             PARAVIEW_DUMP_ATTR
!
        implicit none
        real*8  :: perc_loc
        integer :: nlevel_loc,idec
        character(len=2) :: VLEVEL_save
        logical :: PARAVIEW_DUMP_GEOM_save,PARAVIEW_DUMP_ATTR_save
!        
!------------------------------------------------------------------------------        
!
!       save Paraview parameters
        VLEVEL_save             = VLEVEL
        PARAVIEW_DUMP_GEOM_save = PARAVIEW_DUMP_GEOM
        PARAVIEW_DUMP_ATTR_save = PARAVIEW_DUMP_ATTR
!
!       modify Paraview parameters 
        VLEVEL             = '0'
        PARAVIEW_DUMP_GEOM = .TRUE.
        PARAVIEW_DUMP_ATTR = .FALSE. 
!
!       dump initial mesh to Paraview
        call paraview_driver
!
!       -- INTERACTIVE MODE --        
        if (ITEST == 0) then
!
          idec=1
!
!         display in infinite loop
          do while(idec /= 0)
!          
            write(*,*)'Set: percentage (0,1) ='
            read( *,*) perc_loc
            write(*,*)'Set: refinement levels ='
            read( *,*) nlevel_loc
            call random_refine(perc_loc,nlevel_loc)
!            
            write(*,*)'Clean up data structure? 0 - No ; 1 - Yes'
            read( *,*) idec
!
!           clean up
            if (idec /= 0) then
              call deallocds
              call hp3gen(trim(FILE_PHYS))
            endif
!            
            write(*,*)'0 - Return ; 1 - Continue testing'
            read( *,*) idec
!
!         end infinite loop
          enddo
! 
!       -- TESTING MODE --                
        else ; call random_refine(PERC,NLEVEL)
        endif
!        
!       reset Paraview parameters
        VLEVEL             = VLEVEL_save
        PARAVIEW_DUMP_GEOM = PARAVIEW_DUMP_GEOM_save
        PARAVIEW_DUMP_ATTR = PARAVIEW_DUMP_ATTR_save
!
!
endsubroutine test_4
