!-------------------------------------------------------------------------
subroutine point_2surf(Np)
!-------------------------------------------------------------------------
! MODULES
  use control
  use U2D
!-------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------
! DUMMY ARGUMENTS
  integer, intent(in)      :: Np
!-------------------------------------------------------------------------
! VARIABLES
  integer, dimension(3)  :: surfs
  real*8, dimension(3,3) :: GRAD
  real*8                 :: prod
!-------------------------------------------------------------------------
! printing flag (0,1,2)  
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'point_2surf: reprojecting point Np = ',Np
#endif    
! ..store surface numbers
    surfs(1:2) = POINT_TYPE(2:3,Np)
    surfs(3)   = NRSURFS
! ..compute normals    
    call surf(surfs(1),POINTS(Np)%Rdata(1:3), VOID_1,GRAD(1:3,1))
    call surf(surfs(2),POINTS(Np)%Rdata(1:3), VOID_1,GRAD(1:3,2))
    call normalize(GRAD(1:3,1))
    call normalize(GRAD(1:3,2))
    call cross_product(GRAD(1:3,1),GRAD(1:3,2), GRAD(1:3,3))
    call normalize(GRAD(1:3,3))
    call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,3), prod)
    prod = abs(prod)
! ..if prod is close to 0, redefine point as 1_surf   
    if (prod .lt. GEOM_TOL) then
      POINT_TYPE(1,Np) = 1
      call point_1surf(Np)
! ..treat point as 2_surf      
    else  
      SURFACES(surfs(3))%Rdata(1:3) = POINTS(Np)%Rdata(1:3)
      SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,3)
      VOID_1 = 0.d0;  VOID_2 = 0.d0
      call mnewt(1,surfs,VOID_1,VOID_2,POINTS(Np)%Rdata(1:3),VOID_2, POINTS(Np)%Rdata(1:3))
#if I_PRINT >= 2
      write(*,*)'point_2surf: Newton method applied.'
#endif
    endif
#if I_PRINT >= 1
    write(*,*)'point_2surf: done!'
#endif
!
end subroutine point_2surf
