!-------------------------------------------------------------------------
subroutine point_1surf(Np)
!-------------------------------------------------------------------------
! MODULES
  use GMP
  use U2D
!-------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------
! DUMMY ARGUMENTS
  integer, intent(in) :: Np
!-------------------------------------------------------------------------
! VARIABLES
  integer, dimension(3)  :: surfs
  real*8, dimension(3,3) :: GRAD
  integer                :: i,i1,i2
!-------------------------------------------------------------------------
! FUNCTIONS
  integer :: my_mod
!-------------------------------------------------------------------------
! printing flag (0,1,2)  
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'point_1surf: projecting on 1 surface Np = ',Np
#endif    
! ..store point surface number and extra planes surface number
    surfs(1) = POINT_TYPE(2,Np)
    surfs(2) = NRSURFS - 1
    surfs(3) = NRSURFS    
#if I_PRINT >= 2
    write(*,*)'point_1surf: surface numbers = ',surfs(1:3)
#endif    
! ..compute normal
    call surf(surfs(1),POINTS(Np)%Rdata(1:3), VOID_1,GRAD(1:3,1))
    call normalize(GRAD(1:3,1))
! ..rotate normal by 90 deg    
    do i = 1, 3
      i1 = my_mod(i + 1,3)
      i2 = my_mod(i + 2,3)
! ....choose direction of rotation based on NORMAL(i)                    
      if (GRAD(i,1) .ne. 0d0) then 
        GRAD(i1,2) =  GRAD(i,1)
        GRAD(i,2)  = -GRAD(i1,1)
        GRAD(i2,2) = 0.d0
        exit
      endif
    enddo
    call normalize(GRAD(1:3,2))
! ..determine 3rd vector
    call cross_product(GRAD(1:3,1),GRAD(1:3,2), GRAD(1:3,3))
! ..create 2 extra planes
    SURFACES(surfs(2))%Rdata(1:3) = POINTS(Np)%Rdata(1:3)
    SURFACES(surfs(2))%Rdata(4:6) = GRAD(1:3,2)
!
    SURFACES(surfs(3))%Rdata(1:3) = POINTS(Np)%Rdata(1:3)
    SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,3)
! ..apply Newton-Rapson method; use point coordinats as initial guess
    VOID_1 = 0.d0;  VOID_2 = 0.d0
    call mnewt(1,surfs,VOID_1,VOID_2,POINTS(Np)%Rdata(1:3),VOID_2, POINTS(Np)%Rdata(1:3))
#if I_PRINT >= 2
    write(*,*)'point_1surf: Newton method applied.'
#endif
#if I_PRINT >= 1
    write(*,*)'point_1surf: done!'
#endif    
!
end subroutine point_1surf
