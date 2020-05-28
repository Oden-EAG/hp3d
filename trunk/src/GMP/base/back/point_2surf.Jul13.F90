subroutine point_2surf(Np,Iv,Nt,NSURFS)
!
! MODULES
  use kinds
  use control
  use GMP
!-------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------
! DUMMY ARGUMENTS
  integer, intent(in)               :: Np
  integer, intent(in)               :: Nt
  integer, intent(in)               :: Iv
  integer, dimension(2), intent(in) :: NSURFS
!-------------------------------------------------------------------------
! VARIABLES
  integer, dimension(3) :: surfs
  integer :: nc,np_s
  real(DP), dimension(3,3) :: GRAD
  real(DP), dimension(2) :: M_PROD
  integer :: status
  real(DP), dimension(3) :: void_1
  real(DP), dimension(4) :: void_2
  real(DP), dimension(1) :: temp
!-------------------------------------------------------------------------
! FUNCTIONS
  integer :: mod3
!-------------------------------------------------------------------------
!
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'point_2surf: reprojecting point np = ',Np
#endif
! ..store first 2 surface numbers
    surfs(1:2) = NSURFS(1:2)
! ..create 1 extra plane and store its surface number
    NRSURFS = NRSURFS + 1;  surfs(3) = NRSURFS
! ..get 1st curve number
    nc = abs(TRIANGLES(Nt)%EdgeNo(Iv))
#if I_PRINT >= 2
    write(*,*)'point_2surf: nc_1 = ',nc
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
1   format(' point_2surf: grad_1 =',3(E12.5,2X))
#endif
! ..get 2nd curve number
    nc = abs(TRIANGLES(Nt)%EdgeNo(mod3(Iv + 2)))
! ..get curve start point
    np_s = CURVES(nc)%EndPoNo(1)
! ..account for orientation, compute gradient and normalize it
    if (np .eq. np_s) then
      call curve(nc,0.d0, void_1,GRAD(1:3,2))
    else
      call curve(nc,1.d0, void_1,GRAD(1:3,2))
    endif
    call normalize(GRAD(1:3,2))
#if I_PRINT >= 2
    write(*,2) GRAD(1:3,2)
2   format(' point_2surf: grad_2 =',3(E12.5,2X))
#endif
! ..compute normals to 2 sufaces
    call surf(surfs(1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,3))
    call surf(surfs(2),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,4))
! ..compute absolute value of 2  mixed products
    call mixed_product(GRAD(1:3,1),GRAD(1:3,3),GRAD(1:3,4), M_PROD(1))
    M_PROD(1) = abs(M_PROD(1))
    call mixed_product(GRAD(1:3,2),GRAD(1:3,3),GRAD(1:3,4), M_PROD(2))
    M_PROD(2) = abs(M_PROD(2))
! ..compute maximum and check if it is an appropriate value
    if (maxval(M_PROD(1:2)) .lt. 0.3d0) then
      write(*,*)'point_2surf: warning, mixed product less than 0.03.'
    endif
! ..create plane
    SURFACES(surfs(3))%Type = 'VecPt'
    allocate (SURFACES(surfs(3))%Rdata(6), STAT = status)
    if (status .ne. 0 ) then
      write(*,*)'point_2surf: extra plane not allocated.'
      stop
    endif
    SURFACES(surfs(3))%Rdata(1:3) = POINTS(np)%Rdata(1:3)
! ..choose appropriate normal
    temp = maxloc(M_PROD(1:2))
    if (temp(1) .eq. 1) then
     SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,1)
    elseif (temp(1) .eq. 2) then
      SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,2)
    else
      write(*,*)'point_2surf: something is wrong...'
      stop
    endif
#if I_PRINT >= 2
    write(*,*)'point_2surf: extra plane defined.'
#endif
! ..apply Newton-Rapson method; use points coordinates as initial guess
    void_1 = 0.d0;  void_2 = 0.d0
    call mnewt(1,surfs,void_1,void_2,POINTS(np)%Rdata(1:3),void_2, POINTS(np)%Rdata(1:3))
#if I_PRINT >= 2
    write(*,*)'point_2surf: Newton method applied.'
#endif
! ..delete extra surface
    SURFACES(surfs(3))%Type = 'Void'
    deallocate(SURFACES(surfs(3))%Rdata, STAT = status)
    if (status .ne. 0 ) then
      write(*,*)'point_2surf: plane not deallocated.'
      stop
    endif
    NRSURFS = NRSURFS - 1
#if I_PRINT >= 2
    write(*,*)'point_2surf: extra surface deallocated.'
#endif
! ..delete point type
    POINTS(Np)%Type = 'Regular'
    deallocate(POINTS(Np)%Idata, STAT = status)
    if (status .ne. 0) then
      write(*,*)'point_2surf: point Idata not deallocated.'
      stop
    endif
#if I_PRINT >= 2
    write(*,*)'point_2surf: point type redefined.'
#endif
#if I_PRINT >= 1
    write(*,*)'point_2surf: done!'
#endif
!
end subroutine point_2surf
