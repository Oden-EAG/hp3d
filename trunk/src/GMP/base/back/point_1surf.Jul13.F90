!-------------------------------------------------------------------------
subroutine point_1surf(Np,Iv,Nt)
!-------------------------------------------------------------------------
! MODULES
  use kinds
  use control
  use GMP
!-------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------
! DUMMY ARGUMENTS
  integer, intent(in) :: Np
  integer, intent(in) :: Nt
  integer, intent(in) :: Iv
!-------------------------------------------------------------------------
! VARIABLES
  integer, dimension(3)    :: surfs
  integer                  :: nc,np_s
  real(DP), dimension(3,3) :: GRAD
  real(DP)                 :: m_prod
  integer                  :: status
  real(DP), dimension(3)   :: void_1
  real(DP), dimension(4)   :: void_2
!-------------------------------------------------------------------------
! FUNCTIONS
  integer :: mod3
!-------------------------------------------------------------------------
! PARAMETERS
  real(DP), parameter :: eps = 0.3d0
!-------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'point_1surf: projecting on 1 surface Np = ',Np
#endif
! ..store point surface number
    surfs(1) = POINTS(Np)%Idata(2)
#if I_PRINT >= 2
    write(*,*)'point_1surf: surface number = ',surfs(1)
#endif
! ..create 2 extra planes and store their surface number
    NRSURFS = NRSURFS + 1;  surfs(2) = NRSURFS
    NRSURFS = NRSURFS + 1;  surfs(3) = NRSURFS
! ..get 1st curve number
    nc = abs(TRIANGLES(Nt)%EdgeNo(Iv))
#if I_PRINT >= 2
    write(*,*)'point_1surf: nc_1 =',nc
#endif
! ..get curve start point
    np_s = CURVES(nc)%EndPoNo(1)
! ..account for orientation and compute gradient
    if (Np .eq. np_s) then
      call curve(nc,0.d0, void_1,GRAD(1:3,1))
    else
      call curve(nc,1.d0, void_1,GRAD(1:3,1))
    endif
    call normalize(GRAD(1:3,1))
#if I_PRINT >= 2
    write(*,1) GRAD(1:3,1)
1   format(' point_1surf: grad_1 =',3(E12.5,2X))
#endif
! ..get 2nd curve number
    nc = abs(TRIANGLES(Nt)%EdgeNo(mod3(Iv + 2)))
#if I_PRINT >= 2
    write(*,*)'point_1surf: nc_2 =',nc
#endif
! ..get curve start point
    np_s = CURVES(nc)%EndPoNo(1)
! ..compute gradient
    if (Np .eq. np_s) then
      call curve(nc,0.d0, void_1,GRAD(1:3,2))
    else
      call curve(nc,1.d0, void_1,GRAD(1:3,2))
    endif
    call normalize(GRAD(1:3,2))
#if I_PRINT >= 2
    write(*,2) GRAD(1:3,2)
2   format(' point_1surf: grad_2 =',3(E12.5,2X))
#endif
! ..check orthogonality
    call surf(surfs(1),POINTS(Np)%Rdata(1:3), void_1,GRAD(1:3,3))
    call normalize(GRAD(1:3,3))
    call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,3), m_prod)
    m_prod = abs(m_prod)
    if (m_prod .lt. eps) then
      write(*,*)'point_1surf: warning, m_prod = ',m_prod
    endif
! ..create 2 extra planes
    SURFACES(surfs(2))%Type = 'VecPt'
    allocate (SURFACES(surfs(2))%Rdata(6), STAT = status)
    if (status .ne. 0) then
      write(*,*)'point_1surf: 1st plane not allocated.'
      stop
    endif
    SURFACES(surfs(2))%Rdata(1:3) = POINTS(Np)%Rdata(1:3)
    SURFACES(surfs(2))%Rdata(4:6) = GRAD(1:3,1)
#if I_PRINT >= 2
    write(*,*)'point_1surf: 1st extra plane allocated.'
#endif
    SURFACES(surfs(3))%Type = 'VecPt'
    allocate (SURFACES(surfs(3))%Rdata(6), STAT = status)
    if (status .ne. 0) then
      write(*,*)'point_1surf: 2nd plane not allocated.'
    endif

!    if (Np .eq. 207)  call print_GMP
!    write(*,*)POINTS(Np)%Rdata(1:3)

    SURFACES(surfs(3))%Rdata(1:3) = POINTS(Np)%Rdata(1:3)
    SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,2)
#if I_PRINT >= 2
    write(*,*)'point_1surf: 2nd extra plane allocated.'
#endif
! ..apply Newton-Rapson method; use point coordinats as initial guess
    void_1 = 0.d0; void_2 = 0.d0
    call mnewt(1,surfs,void_1,void_2,POINTS(Np)%Rdata(1:3),void_2, POINTS(Np)%Rdata(1:3))
#if I_PRINT >= 2
    write(*,*)'point_1surf: Newton method applied.'
#endif
! ..delete extra surfaces
    SURFACES(surfs(2))%Type = 'Void'
    SURFACES(surfs(3))%Type = 'Void'
    deallocate(SURFACES(surfs(2))%Rdata, STAT = status)
    if (status .ne. 0) then
      write(*,*)'point_1surf: 1st plane not deallocated.'
      stop
    endif
#if I_PRINT >= 2
    write(*,*)'point_1surf: 1st extra plane deallocated.'
#endif
    deallocate(SURFACES(surfs(3))%Rdata, STAT = status)
    if (status .ne. 0) then
      write(*,*)'point_1surf: 2nd plane not deallocated.'
      stop
    endif
#if I_PRINT >= 2
    write(*,*)'point_1surf: 2nd extra plane deallocated.'
#endif
    NRSURFS = NRSURFS - 2
#if I_PRINT >= 1
    write(*,*)'point_1surf: done!'
#endif
!
end subroutine point_1surf
