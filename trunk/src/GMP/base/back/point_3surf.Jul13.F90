subroutine point_3surf(Np,Iv,Nt)
!
! MODULES
  use kinds
  use control
  use GMP
  use U2D
!-------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!-------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------
! DUMMY ARGUMENTS
  integer, intent(in)               :: Np
  integer, intent(in)               :: Nt
  integer, intent(in)               :: Iv
!-------------------------------------------------------------------------
! VARIABLES
  integer, dimension(3) :: surfs
  integer :: nc,np_s,i,j,k,l
  real(DP), dimension(3,4) :: GRAD
  integer :: status
  real(DP), dimension(3) :: void_1
  real(DP), dimension(4) :: void_2
  integer, dimension(1) :: temp
#if I_PRINT >= 2
  integer :: combination
#endif
!-------------------------------------------------------------------------
!
#if I_PRINT >= 1
  write(*,*)'point_3surf: reprojecting point on 3 surfaces.'
#endif
#if I_PRINT >= 2
  write(*,1) POINTS(Np)%Idata(1)
1 format(' point_3surf: point lies on ',I2,' surfaces.')
  write(*,*)'point_3surf: surface numbers = ',POINTS(Np)%Idata(2:POINTS(Np)%Idata(1) + 1)
#endif
! loop through all possible choices of 3 surfaces
  l = 0
  do i = 1, (POINTS(Np)%Idata(1) - 2)
    do j = (i + 1), (POINTS(Np)%Idata(1) - 1)
      do k = (j + 1), POINTS(Np)%Idata(1)
! ......increment counter
        l = l + 1
! ......compute gradients of 3 selected surfaces and normalize them
        call surf(POINTS(Np)%Idata(i + 1),POINTS(Np)%Rdata(1:3), void_1,GRAD(1:3,1))
        call surf(POINTS(Np)%Idata(j + 1),POINTS(Np)%Rdata(1:3), void_1,GRAD(1:3,2))
        call surf(POINTS(Np)%Idata(k + 1),POINTS(Np)%Rdata(1:3), void_1,GRAD(1:3,3))
        call normalize(GRAD(1:3,1))
        call normalize(GRAD(1:3,2))
        call normalize(GRAD(1:3,3))
! ......compute absolute value of mixed product and store pertaining surfaces
        call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,3), M_PROD(1,l))
        M_PROD(1,l) = abs(M_PROD(1,l))
        M_PROD(2,l) = POINTS(Np)%Idata(i + 1)
        M_PROD(3,l) = POINTS(Np)%Idata(j + 1)
        M_PROD(4,l) = POINTS(Np)%Idata(k + 1)
      enddo
    enddo
  enddo
#if I_PRINT >= 2
! run additional check
  if (combination(POINTS(Np)%Idata(1),3) .ne. l) then
    write(*,*)'point_3surf: something is wrong...'
    write(*,*) 'l = ',l
    write(*,*) 'combinations = ',combination(POINTS(Np)%Idata(1),3)
    stop
  endif
#endif
! select biggest prod, which should be the closest to 1 (orthonormal tern)
  if (maxval(M_PROD(1,1:l)) .gt. 0.1) then
! ..get position of greatest mixed product
    temp = maxloc(M_PROD(1,1:l))
#if I_PRINT >= 2
    write(*,*)'point_3surf: found 3 appropriate surfaces = ',M_PROD(2:4,temp(1))
#endif
! ..apply Newton method
    void_1 = 0.d0;  void_2 = 0.d0
    call mnewt(1,int(M_PROD(2:4,temp(1))),void_1,void_2,POINTS(Np)%Rdata(1:3),void_2, POINTS(Np)%Rdata(1:3))
#if I_PRINT >= 2
    write(*,*)'point_3surf: Newton method applied.'
#endif
! ..reset point type
    POINTS(Np)%Type = 'Regular'
    deallocate(POINTS(Np)%Idata, STAT = status)
    if (status .ne. 0) then
      write(*,*)'point_3surf: point Idata not deallocated.'
      stop
    endif
#if I_PRINT >= 2
    write(*,*)'point_3surf: restored point type.'
#endif
! else treat point as a 2 surf point
  else
#if I_PRINT >= 2
    write(*,*)'point_3surf: redefining point as 2_surf point.'
#endif
! ..loop through all possible choices of 2 surfaces
    l = 0
    do i = 1, (POINTS(Np)%Idata(1) - 1)
      do j = 1, POINTS(Np)%Idata(1)
! ......increment counter
        l = l + 1
! ......compute gradients of 2 selected surfaces and normalize them
        call surf(POINTS(Np)%Idata(i + 1),POINTS(Np)%Rdata(1:3), void_1,GRAD(1:3,1))
        call surf(POINTS(Np)%Idata(j + 1),POINTS(Np)%Rdata(1:3), void_1,GRAD(1:3,2))
        call normalize(GRAD(1:3,1))
        call normalize(GRAD(1:3,2))
! ......compute absolute value of scalar product and store pertainig surfaces
        call scalar_product(GRAD(1:3,1),GRAD(1:3,2), S_PROD(1,l))
        S_PROD(1,l) = abs(S_PROD(1,l))
        S_PROD(2,l) = POINTS(Np)%Idata(i + 1)
        S_PROD(3,l) = POINTS(Np)%Idata(j + 1)
      enddo
    enddo
! ..choose smallest scalar product in order to retrive surfaces number
    temp = minloc(S_PROD(1,1:l))
    surfs(1:2) = S_PROD(2:3,temp(1))
#if I_PRINT >= 2
    write(*,*)'point_3surf: appropriate surfaces = ',surfs(1:2)
#endif
    call point_2surf(Np,Iv,Nt,surfs(1:2))
! end treat as a 2 surf point
  endif
#if I_PRINT >= 1
  write(*,*)'point_3surf: done!'
#endif
!
end subroutine point_3surf
