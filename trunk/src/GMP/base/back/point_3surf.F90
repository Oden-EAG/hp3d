!-----------------------------------------------------------------------------------------------
subroutine point_3surf(Np)
!-----------------------------------------------------------------------------------------------
! MODULES
  use control
  use U2D
!-----------------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!-----------------------------------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------
! DUMMY ARGUMENTS
  integer, intent(in)    :: Np
!-----------------------------------------------------------------------------------------------
! VARIABLES
  integer, dimension(3)  :: surfs
  integer                :: i,j,k,l
  real*8, dimension(3,3) :: GRAD
  integer, dimension(1)  :: temp
#if I_PRINT >= 2
  integer :: combination
#endif
!-----------------------------------------------------------------------------------------------
!
#if I_PRINT >= 1
    write(*,*)'point_3surf: reprojecting point on 3 surfaces.'
#endif
#if I_PRINT >= 2
    write(*,1) POINT_TYPE(1,Np)
1   format(' point_3surf: point lies on ',I2,' surfaces.')
    write(*,*)'point_3surf: surface numbers = ',POINT_TYPE(2:(POINT_TYPE(1,Np) + 1,Np)
#endif
! ..loop through all possible choices of 3 surfaces
    l = 0
    do i = 1, (POINT_TYPE(1,Np) - 2)
      do j = (i + 1), (POINT_TYPE(1,Np) - 1)
        do k = (j + 1), POINT_TYPE(1,Np)
! ........increment counter
          l = l + 1
! ........compute gradients of 3 selected surfaces and normalize them
          call surf(POINT_TYPE(i + 1,Np),POINTS(Np)%Rdata(1:3), VOID_1,GRAD(1:3,1))
          call surf(POINT_TYPE(j + 1,Np),POINTS(Np)%Rdata(1:3), VOID_1,GRAD(1:3,2))
          call surf(POINT_TYPE(k + 1,Np),POINTS(Np)%Rdata(1:3), VOID_1,GRAD(1:3,3))
          call normalize(GRAD(1:3,1))
          call normalize(GRAD(1:3,2))
          call normalize(GRAD(1:3,3))
! ........compute absolute value of mixed product and store pertaining surfaces
          call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,3), M_PROD(1,l))
          M_PROD(1,l) = abs(M_PROD(1,l))
          M_PROD(2,l) = POINT_TYPE(i + 1,Np)
          M_PROD(3,l) = POINT_TYPE(j + 1,Np)
          M_PROD(4,l) = POINT_TYPE(k + 1,Np)
        enddo
      enddo
    enddo
#if I_PRINT >= 2
! ..run additional check
    if (combination(POINT_TYPE(1,Np),3) .ne. l) then
      write(*,*)'point_3surf: something is wrong...'
      write(*,*) 'l = ',l
      write(*,*) 'combinations = ',combination(POINT_TYPE(1,Np),3)
      stop
    endif
#endif
! ..select biggest prod, which should be the closest to 1 (orthonormal tern)
    if (maxval(M_PROD(1,1:l)) .gt. 0.1) then
! ....get position of greatest mixed product
      temp = maxloc(M_PROD(1,1:l))
#if I_PRINT >= 2
      write(*,*)'point_3surf: found 3 appropriate surfaces = ',M_PROD(2:4,temp(1))
#endif
! ....apply Newton method
      VOID_1 = 0.d0;  VOID_2 = 0.d0
      call mnewt(1,int(M_PROD(2:4,temp(1))),VOID_1,VOID_2,POINTS(Np)%Rdata(1:3),VOID_2, &
                                                                      POINTS(Np)%Rdata(1:3))
#if I_PRINT >= 2
      write(*,*)'point_3surf: Newton method applied.'
#endif
! ..else treat point as a 2 surf point
    else
#if I_PRINT >= 2
      write(*,*)'point_3surf: redefining point as 2_surf point.'
#endif
      POINT_TYPE(1,Np) = 2
! ....loop through all possible choices of 2 surfaces
      l = 0
      do i = 1, (POINT_TYPE(1,Np) - 1)
        do j = 1, POINT_TYPE(1,Np)
! ........increment counter
          l = l + 1
! ........compute gradients of 2 selected surfaces and normalize them
          call surf(POINT_TYPE(i + 1,Np),POINTS(Np)%Rdata(1:3), VOID_1,GRAD(1:3,1))
          call surf(POINT_TYPE(j + 1,Np),POINTS(Np)%Rdata(1:3), VOID_1,GRAD(1:3,2))
          call normalize(GRAD(1:3,1))
          call normalize(GRAD(1:3,2))
! ........compute absolute value of scalar product and store pertainig surfaces
          call scalar_product(GRAD(1:3,1),GRAD(1:3,2), S_PROD(1,l))
          S_PROD(1,l) = abs(S_PROD(1,l))
          S_PROD(2,l) = POINT_TYPE(i + 1,Np)
          S_PROD(3,l) = POINT_TYPE(j + 1,Np)
        enddo
      enddo
! ....choose smallest scalar product in order to retrive surfaces number
      temp = minloc(S_PROD(1,1:l))
! ....store surface numbers
      POINT_TYPE(2:3,Np) = S_PROD(2:3,temp(1))
#if I_PRINT >= 2
      write(*,*)'point_3surf: appropriate surfaces = ',S_PROD(2:3,temp(1))
#endif
      call point_2surf(Np)
! ..end treat as a 2 surf point
    endif
#if I_PRINT >= 1
    write(*,*)'point_3surf: done!'
#endif
!
end subroutine point_3surf
!-----------------------------------------------------------------------------------------------
