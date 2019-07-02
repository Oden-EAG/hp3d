!========================================================================================
!  subroutine - upgrade2double                                                          |
!========================================================================================
!  LATEST REVISION: Jul 09                                                              |
!
!  PURPOSE: routine upgrades to double precision geometry using Newton
!           method iterations
!
!  REMARKS: new general routine!
!========================================================================================
subroutine upgrade2double_NEW
!----------------------------------------------------------------------------------------
! MODULES
  use U2D
!----------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 2  
!----------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------
! VARIABLES
  integer  :: np
#if I_PRINT >= 2  
  integer  :: n_0surf,n_1surf,n_2surf,n_3surf
  real(DP), dimension(3) :: aux_v
  real(DP)               :: aux_old,aux_new
#endif
!----------------------------------------------------------------------------------------
! printing flag (0,1,2)
!
#if I_PRINT >= 1
    write(*,*)'upgrade2double_NEW: upgrading to double precision geometry...'
#endif
    call generate_POINT_TYPE
#if I_PRINT >= 2
    write(*,*)'upgrade2double_NEW: POINT_TYPE generated.'
#endif    
    call allocate_S_PROD
#if I_PRINT >= 2
    write(*,*)'upgrade2double_NEW: S_PROD allocated.'
#endif    
    call allocate_M_PROD
#if I_PRINT >= 2
    write(*,*)'upgrade2double_NEW: M_PROD allocated.'
#endif    
    call allocate_EXTRA_PLANES
#if I_PRINT >= 2
    write(*,*)'upgrade2double_NEW: EXTRA_PLANES allocated.'
#endif    
!
#if I_PRINT >= 2
    n_0surf = 0;  n_1surf = 0;  n_2surf = 0;  n_3surf = 0
    aux_old = 0.d0
#endif    
! ..loop through points
    do np = 1, NRPOINT
#if I_PRINT >= 2
      write(*,*)'****************************************************************'
      write(*,2) np, POINT_TYPE(1,np)
2     format(' upgrade2double_NEW: np  = ',I6,' ; number of surfaces = ',I2)
      aux_v(1:3) = POINTS(np)%Rdata(1:3)
#endif  
      select case (POINT_TYPE(1,np))
! ......nothing to do!      
        case (0)
#if I_PRINT >= 2
          n_0surf = n_0surf + 1       
          write(*,*)'upgrade2double_NEW: 0_surf point, n_0surf = ',n_0surf
#endif
! ......point lies on ONE surface        
        case (1)
#if I_PRINT >= 2
          n_1surf = n_1surf + 1       
          write(*,*)'upgrade2double_NEW: 1_surf point, n_1surf = ',n_1surf
#endif
          call point_1surf(np)
! ......point lies on TWO surfaces            
        case (2)      
#if I_PRINT >= 2
          n_2surf = n_2surf + 1
          write(*,*)'upgrade2double_NEW: 2_surf point, n_2surf = ',n_2surf
#endif
          call point_2surf(np)
! ......point lies on THREE OR MORE surfaces            
        case default
#if I_PRINT >= 2
          n_3surf = n_3surf + 1
          write(*,*)'upgrade2double_NEW: point lies on 3 surfaces or more, n_3surf = ',n_3surf
#endif
          call point_3surf(np)
      end select
#if I_PRINT >= 2
      aux_v(1:3) = aux_v(1:3) - POINTS(np)%Rdata(1:3)
      call norm(aux_v, aux_new)
      if (aux_new .gt. aux_old) aux_old = aux_new
#endif
! ..end of loop through points
    enddo
!
! ..deallocate
    call deallocate_EXTRA_PLANES
    call deallocate_ARRAYS
#if I_PRINT >= 2
    write(*,*)'****************************************************************'
    write(*,1) n_0surf,n_1surf,n_2surf,n_3surf
1   format(' upgrade2double_NEW: n_0surf, n_1surf, n_2surf, n_3surf = '4(6I,1X))
    write(*,*)'upgrade2double_NEW: max shift = ',aux_old
#endif  
#if I_PRINT >= 1
write(*,*)'upgrade2double_NEW: done!'
#endif  
!
end subroutine upgrade2double_NEW
