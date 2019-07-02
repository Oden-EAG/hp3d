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
  use kinds
  use control
  use GMP
  use U2D
!----------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 2  
!----------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------
! VARIABLES
  real(DP)                 :: geom_tol_orig
  real(DP), dimension(3)   :: x_0, void_1
  real(DP), dimension(4)   :: void_2
  real(DP), dimension(3,4) :: GRAD
  integer                  :: status
  integer                  :: nt, iv, np, nc, np_s, i, j, k, l
  integer, dimension(3)    :: surfs
#if I_PRINT >= 2  
  integer                  :: n_1surf, n_2surf, n_3surf
  real(DP) :: norm_old,norm_new
  real(DP), dimension(3) :: aux
#endif
! FUNCTIONS
  integer                  :: mod3,norm
!----------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 2
!
! ..set GEOM_TOL to a fairly big value in order not to run into warning messages
    geom_tol_orig = GEOM_TOL;  GEOM_TOL = 0.1D-3
#if I_PRINT >= 1
    write(*,*)'upgrade2double_NEW: upgrading to double precision geometry...'
#endif
#if I_PRINT >= 2
    write(*,*)'upgrade2double_NEW: setting GEOM_TOL = ',GEOM_TOL
#endif
! ..update points type
    call update_point_type
#if I_PRINT >= 2
    write(*,*)'upgrade2double: points type updated.'
#endif
    call determine_N_SURFS
    call allocate_S_PROD
    call allocate_M_PROD
!
#if I_PRINT >= 2
    n_1surf = 0;  n_2surf = 0;  n_3surf = 0; norm_old = 0.d0
#endif    
! ..loop through triangles
    do nt = 1, NRTRIAN
! ....cycle if not on a surface
      if (TRIANGLES(nt)%Type .ne. 'PTITri') cycle
#if I_PRINT >= 2
      write(*,*)'****************************************************************'
      write(*,*)'upgrade2double: nt, type =',nt,TRIANGLES(nt)%Type
#endif  
! ....loop through vertices
      do iv = 1, 3
! ......get point number
        np = TRIANGLES(nt)%VertNo(iv)
! ......cycle if point has already been visited
        if (POINTS(np)%Type .eq. 'Regular') cycle
#if I_PRINT >= 2
        aux = POINTS(Np)%Rdata(1:3)
#endif        
! ......select number of surfaces
#if I_PRINT >= 2
        write(*,*)'upgrade2double: first visit to point np =',np
#endif  
        select case (POINTS(np)%Idata(1))
! ........point lies on ONE surface        
          case (1)
#if I_PRINT >= 2
            n_1surf = n_1surf + 1       
            write(*,*)'upgrade2double: 1_surf point, n_1surf = ',n_1surf
#endif
            call point_1surf(np,iv,nt)
! ........point lies on TWO surfaces            
          case (2)      
#if I_PRINT >= 2
            n_2surf = n_2surf + 1
            write(*,*)'upgrade2double: 2_surf point, n_2surf = ',n_2surf
#endif
            surfs(1:2) = POINTS(np)%Idata(2:3)
            call point_2surf(np,iv,nt,surfs(1:2))
! ........point lies on THREE OR MORE surfaces            
          case default
#if I_PRINT >= 2
            n_3surf = n_3surf + 1
            write(*,*)'upgrade2double: point lies on 3 surfaces or more, n_3surf = ',n_3surf
#endif
            call point_3surf(np,iv,nt)
#if I_PRINT >= 2
        aux = aux -  POINTS(Np)%Rdata(1:3)
        norm_new = norm(aux)
        if (norm_new .gt. norm_old) then
          write(*,*)'largest shift so far = ',norm_new
          norm_old = norm_new
        endif        
#endif        
! 
        end select
! ....end of loop through vertices
      enddo
! ..end of loop through triangles 
    enddo
!
! ..reset GEOM_TOL to original value
    GEOM_TOL = geom_tol_orig
#if I_PRINT >= 2
    write(*,*)'****************************************************************'
    write(*,1) n_1surf,n_2surf,n_3surf
1   format(' upgrade2double_NEW: n_1surf, n_2surf, n_3surf = '3(6I,1X))
    write(*,*)'upgrade2double_NEW: largest shift = ',norm_old
#endif  
#if I_PRINT >= 1
write(*,*)'upgrade2double_NEW: done!'
#endif  
!
end subroutine upgrade2double_NEW
