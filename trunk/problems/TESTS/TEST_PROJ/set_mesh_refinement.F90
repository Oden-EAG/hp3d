subroutine set_mesh_refinement
!
      use data_structure3D , only : NRELIS,NODES
!
      implicit none
      integer :: nel,ndom
!-----------------------------------------------------------------      
!
!  ...loop over initial mesh elements
      do nel=1,NRELIS
        select case(NODES(nel)%type)
!       PYRAMID          
!!        case('mdld') ; NODES(nel)%ref_filter=10
        endselect        
      enddo        
!
!      
endsubroutine set_mesh_refinement
