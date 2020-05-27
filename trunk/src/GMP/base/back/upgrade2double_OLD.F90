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
!
subroutine upgrade2double_NEW
!
!----------------------------------------------------------------------------------------
! MODULES
  use kinds
  use control
  use GMP
  use U2D
!----------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------
! VARIABLES
  real(DP)                 :: geom_tol_orig
  real(DP), dimension(3)   :: x_0, void_1
  real(DP), dimension(4)   :: void_2
  real(DP), dimension(3,4) :: GRAD
  integer, dimension(1) :: temp
  integer                               :: status
  integer                  :: max_surfs, nt, iv, np, nc, np_s, i, j, k, l
  integer, dimension(3)    :: surfs
! FUNCTIONS
  integer                  :: mod3, factorial, combinations
!----------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 2
!
! ..set GEOM_TOL to a fairly big value in order not to run into warning messages
    geom_tol_orig = GEOM_TOL;  GEOM_TOL = 0.1D-3
#if I_PRINT >= 1
    write(*,*)'upgrade2double: upgrading to double precision geometry...'
    write(*,*)'upgrade2double: setting GEOM_TOL = ',GEOM_TOL
#endif
! ..initialize void_1, void_2
    void_1 = 0.d0;  void_2 = 0.d0
! ..update points type
    call update_point_type
#if I_PRINT >= 2
    write(*,*)'upgrade2double: points type updated.'
#endif
! ..determine max number of surfaces a point in the mesh lies on
    max_surfs = 0
    do np = 1, NRPOINT
      if (POINTS(np)%Type .eq. 'Regular') cycle
      if (POINTS(np)%Idata(1) .gt. max_surfs) max_surfs = POINTS(np)%Idata(1)
    enddo
! ..allocate S_PROD and M_PROD accordingly
    if (max_surfs .lt. 3) then
      allocate(S_PROD(3,2), STAT = status)
! ....if allocation not successful
      if (status .ne. 0) then
        write(*,*)'upgrade2double: S_PROD not allocated.'
        stop
      endif
      allocate(M_PROD(4,2), STAT = status)
! ....if allocation not successful
      if (status .ne. 0) then
        write(*,*)'upgrade2double: M_PROD not allocated.'
        stop
      endif
! ..at least 3 surfaces
    else
      allocate(S_PROD(3,combinations(max_surfs,2)), STAT = status)
! ....if allocation not successful
      if (status .ne. 0) then
        write(*,*)'upgrade2double: S_PROD not allocated.'
        stop
      endif
      allocate(M_PROD(4,combinations(max_surfs,3)), STAT = status)
! ....if allocation not successful
      if (status .ne. 0) then
        write(*,*)'upgrade2double: M_PROD not allocated.'
        stop
      endif
    endif
#if I_PRINT >= 2
    write(*,*)'upgrade2double: S_PROD, M_PROD allocated.'
#endif
!
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
! ......select number of surfaces
#if I_PRINT >= 2
        write(*,*)'upgrade2double: first visit to point np =',np
#endif
        select case (POINTS(np)%Idata(1))
! ........point lies on ONE surface
          case (1)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: 1_surf point.'
#endif
! ..........store surface number
            surfs(1) = POINTS(np)%Idata(2)
! ..........create 2 extra planes and store their surface number
            NRSURFS = NRSURFS + 1;  surfs(2) = NRSURFS
            NRSURFS = NRSURFS + 1;  surfs(3) = NRSURFS
! ..........get 1st curve number
            nc = abs(TRIANGLES(nt)%EdgeNo(iv))
#if I_PRINT >= 2
  write(*,*)'upgrade2double: nc_1 =',nc
#endif
! ..........get 1st curve start point
            np_s = CURVES(nc)%EndPoNo(1)
! ..........account for orientation and compute gradient
            if (np .eq. np_s) then
              call curve(nc,0.d0, void_1,GRAD(1:3,1))
            else
              call curve(nc,1.d0, void_1,GRAD(1:3,1))
            endif
            call normalize(GRAD(1:3,1))
#if I_PRINT >= 2
            write(*,1) GRAD(1:3,1)
1           format(' upgrade2double: grad_1 =',3(E12.5,2X))
#endif
! ..........get 2nd curve number
            nc = abs(TRIANGLES(nt)%EdgeNo(mod3(iv + 2)))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: nc_2 =',nc
#endif
! ..........get 2nd curve start point
            np_s = CURVES(nc)%EndPoNo(1)
! ..........compute gradient
            if (np .eq. np_s) then
              call curve(nc,0.d0, void_1,GRAD(1:3,2))
            else
              call curve(nc,1.d0, void_1,GRAD(1:3,2))
            endif
            call normalize(GRAD(1:3,2))
#if I_PRINT >= 2
            write(*,1) GRAD(1:3,2)
2           format(' upgrade2double: grad_2 =',3(E12.5,2X))
#endif
! ..........check orthogonality
            call surf(surfs(1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,3))
            call normalize(GRAD(1:3,3))
            call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,3), M_PROD(1,1))
            M_PROD(1,1) = abs(M_PROD(1,1))
            if (M_PROD(1,1) .lt. 0.3d0) then
              write(*,*)'upgrade2double: warning, M_PROD = ',M_PROD(1,1)
            endif
! ..........create 2 extra planes
            SURFACES(surfs(2))%Type = 'VecPt';                      SURFACES(surfs(3))%Type = 'VecPt'
            allocate (SURFACES(surfs(2))%Rdata(6));                 allocate (SURFACES(surfs(3))%Rdata(6))
            SURFACES(surfs(2))%Rdata(1:3) = POINTS(np)%Rdata(1:3);  SURFACES(surfs(3))%Rdata(1:3) = POINTS(np)%Rdata(1:3)
            SURFACES(surfs(2))%Rdata(4:6) = GRAD(1:3,1);            SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,2)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: 2 extra planes defined.'
#endif
! ..........apply Newton-Rapson method; use point coordinats as initial guess
            void_1 = 0.d0
            call mnewt(1,surfs,void_1,void_2,POINTS(np)%Rdata(1:3),void_2, POINTS(np)%Rdata(1:3))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: Newton-Rapson applied.'
#endif
! ..........delete extra surfaces
            SURFACES(surfs(2))%Type = 'Void';     SURFACES(surfs(3))%Type = 'Void'
            deallocate(SURFACES(surfs(2))%Rdata); deallocate(SURFACES(surfs(3))%Rdata)
            NRSURFS = NRSURFS - 2
! ........point lies on TWO surfaces
          case (2)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: 2_surf point.'
#endif
! ..........store surfaces number
            surfs(1:2) = POINTS(np)%Idata(2:3)
! ..........create 1 extra plane and store its surface number
            NRSURFS = NRSURFS + 1;  surfs(3) = NRSURFS
! ..........get 1st curve number
            nc = abs(TRIANGLES(nt)%EdgeNo(iv))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: nc_1 = ',nc
#endif
! ..........get curve start point
            np_s = CURVES(nc)%EndPoNo(1)
! ..........account for orientation and compute gradient
            if (np .eq. np_s) then
              call curve(nc,0.d0, void_1,GRAD(1:3,1))
            else
              call curve(nc,1.d0, void_1,GRAD(1:3,1))
            endif
            call normalize(GRAD(1:3,1))
#if I_PRINT >= 2
            write(*,1) GRAD(1:3,1)
#endif
! ..........get 2nd curve number
            nc = abs(TRIANGLES(nt)%EdgeNo(mod3(iv + 2)))
! ..........get curve start point
            np_s = CURVES(nc)%EndPoNo(1)
! ..........account for orientation, compute gradient and normalize it
            if (np .eq. np_s) then
              call curve(nc,0.d0, void_1,GRAD(1:3,2))
            else
              call curve(nc,1.d0, void_1,GRAD(1:3,2))
            endif
            call normalize(GRAD(1:3,2))
#if I_PRINT >= 2
            write(*,2) GRAD(1:3,2)
#endif
! ..........compute normals to 2 sufaces
            call surf(surfs(1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,3))
            call surf(surfs(2),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,4))
! ..........compute absolute value of 2  mixed products
            call mixed_product(GRAD(1:3,1),GRAD(1:3,3),GRAD(1:3,4), M_PROD(1,1))
            M_PROD(1,1) = abs(M_PROD(1,1))
            call mixed_product(GRAD(1:3,2),GRAD(1:3,3),GRAD(1:3,4), M_PROD(1,2))
            M_PROD(1,2) = abs(M_PROD(1,2))
! ..........compute maximum and check if it is an appropriate value
            if (maxval(M_PROD(1,1:2)) .lt. 0.3d0) then
              write(*,*)'upgrade2double: warning, mixed product less than 0.03.'
            endif
! ..........create plane
            SURFACES(surfs(3))%Type = 'VecPt'
            allocate (SURFACES(surfs(3))%Rdata(6))
            SURFACES(surfs(3))%Rdata(1:3) = POINTS(np)%Rdata(1:3)
! ..........choose appropriate normal
            temp = maxloc(M_PROD(1,1:2))
            if (temp(1) .eq. 1) then
              SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,1)
            elseif (temp(1) .eq. 2) then
              SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,2)
            else
              write(*,*)'upgrade2double: something is wrong...'
              stop
            endif

#if I_PRINT >= 2
            write(*,*)'upgrade2double: extra plane defined.'
#endif
! ..........apply Newton-Rapson method; use points coordinates as initial guess
            void_1 = 0.d0
            call mnewt(1,surfs,void_1,void_2,POINTS(np)%Rdata(1:3),void_2, POINTS(np)%Rdata(1:3))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: Newton method applied.'
#endif
! ..........delete extra surface
            SURFACES(surfs(3))%Type = 'Void'
            deallocate(SURFACES(surfs(3))%Rdata)
            NRSURFS = NRSURFS - 1
! ........point lies on THREE OR MORE surfaces
          case default
#if I_PRINT >= 2
            write(*,*)'upgrade2double: point lies on 3 surfaces or more.'
#endif
! ..........loop through all possible choices of 3 surfaces
            l = 0
            do i = 1, (POINTS(np)%Idata(1) - 2)
              do j = 1, (POINTS(np)%Idata(1) - 1)
                do k = 1, POINTS(np)%Idata(1)
! ................increment counter
                  l = l + 1
! ................compute gradients of 3 selected surfaces and normalize them
                  call surf(POINTS(np)%Idata(i + 1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,1))
                  call surf(POINTS(np)%Idata(j + 1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,2))
                  call surf(POINTS(np)%Idata(k + 1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,3))
                  call normalize(GRAD(1:3,1))
                  call normalize(GRAD(1:3,2))
                  call normalize(GRAD(1:3,3))
! ................compute absolute value of mixed product and store pertaining surfaces
                  call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,3), M_PROD(1,l))
                  M_PROD(1,l) = abs(M_PROD(1,l))
                  M_PROD(2,l) = POINTS(np)%Idata(i + 1)
                  M_PROD(3,l) = POINTS(np)%Idata(j + 1)
                  M_PROD(4,l) = POINTS(np)%Idata(k + 1)
                enddo
              enddo
            enddo
! ..........select biggest prod, which should be the closest to 1 (orthonormal tern)
            if (maxval(M_PROD(1,1:l)) .gt. 0.1) then
#if I_PRINT >= 2
              write(*,*)'upgrade2double: found 3 appropriate surfaces.'
#endif
! ............get position of greates mixed product
              temp = maxloc(M_PROD(1,1:l))
! ............apply Newton method
              void_1 = 0.d0;  void_2 = 0.d0
              call mnewt(1,M_PROD(2:4,temp(1)),void_1,void_2,POINTS(np)%Rdata(1:3),void_2, POINTS(np)%Rdata(1:3))
#if I_PRINT >= 2
              write(*,*)'upgrade2double: Newton method applied.'
#endif
! ..........else treat point as a 2 surf point
            else
#if I_PRINT >= 2
              write(*,*)'upgrade2double: redefining point as 2_surf point.'
#endif
! ............loop through all possible choices of 2 surfaces
              l = 0
              do i = 1, (POINTS(np)%Idata(1) - 1)
                do j = 1, POINTS(np)%Idata(1)
! ................increment counter
                  l = l + 1
! ................compute gradients of 2 selected surfaces and normalize them
                  call surf(POINTS(np)%Idata(i + 1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,1))
                  call surf(POINTS(np)%Idata(j + 1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,2))
                  call normalize(GRAD(1:3,1))
                  call normalize(GRAD(1:3,2))
! ................compute absolute value of scalar product and store pertainig surfaces
                  call scalar_product(GRAD(1:3,1),GRAD(1:3,2), S_PROD(1,l))
                  S_PROD(1,l) = abs(S_PROD(1,l))
                  S_PROD(2,l) = POINTS(np)%Idata(i + 1)
                  S_PROD(3,l) = POINTS(np)%Idata(j + 1)
                enddo
              enddo
! ............choose smallest scalar product in order to retrive surfaces number
              temp = minloc(S_PROD(1,1:l))
              surfs(1:2) = S_PROD(2:3,temp(1))
#if I_PRINT >= 2
              write(*,*)'upgrade2double: appropriate surfaces = ',surfs(1:2)
#endif
! ............create 1 extra plane and store its surface number
              NRSURFS = NRSURFS + 1
              if (NRSURFS .gt. MAXSU) then
                write(*,*)'upgrade2double: increase MAXSU = ',MAXSU
                stop
              endif
              surfs(3) = NRSURFS
! ............get 1st curve number
              nc = abs(TRIANGLES(nt)%EdgeNo(iv))
! ............get curve start point
              np_s = CURVES(nc)%EndPoNo(1)
! ............account for orientation, compute gradient and normalize it
              if (np .eq. np_s) then
                call curve(nc,0.d0, void_1,GRAD(1:3,3))
              else
                call curve(nc,1.d0, void_1,GRAD(1:3,3))
              endif
              call normalize(GRAD(1:3,3))
! ............get 2nd curve number
              nc = abs(TRIANGLES(nt)%EdgeNo(mod3(iv + 2)))
! ............get curve start point
              np_s = CURVES(nc)%EndPoNo(1)
! ............account for orientation, compute gradient and normalize it
              if (np .eq. np_s) then
                call curve(nc,0.d0, void_1,GRAD(1:3,4))
              else
                call curve(nc,1.d0, void_1,GRAD(1:3,4))
              endif
              call normalize(GRAD(1:3,4))
! ............compute absolute value of 2  mixed products
              call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,3), M_PROD(1,1))
              M_PROD(1,1) = abs(M_PROD(1,1))
              call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,4), M_PROD(1,2))
              M_PROD(1,2) = abs(M_PROD(1,2))

              write(*,*)'here1'

! ............compute maximum and check if it is an appropriate value
              if (maxval(M_PROD(1,1:2)) .lt. 0.2d0) then
                write(*,*)'upgrade2double: mixed product is small'
              endif

              write(*,*)'here2'
! ............create plane
              SURFACES(surfs(3))%Type = 'VecPt'
              write(*,*)'ns = ',surfs(3)

              deallocate (SURFACES(surfs(3))%Rdata, STAT = status)
              write(*,*)'here4'
              allocate (SURFACES(surfs(3))%Rdata(6), STAT = status)
              write(*,*)'here5'
              if (status .ne. 0) then
                write(*,*)'upgrade2double: allocation failed.'
                stop
              endif
              write(*,*)'here3'
              SURFACES(surfs(3))%Rdata(1:3) = POINTS(np)%Rdata(1:3)

              write(*,*)'here'

! ............choose most appropriate normal vector
              temp = maxloc(M_PROD(1,1:2))
              if (temp(1) .eq. 1) then
                SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,3)
              elseif (temp(1) .eq. 2) then
                SURFACES(surfs(3))%Rdata(4:6) = GRAD(1:3,4)
              else
                write(*,*)'upgrade2double: something is wrong...'
                stop
              endif
#if I_PRINT >= 2
              write(*,*)'upgrade2double: appropriate plane defined.'
#endif
! ............apply Newton-Rapson method; use points coordinates as initial guess
              void_1 = 0.d0;  void_2 = 0.d0
              call mnewt(1,surfs,void_1,void_2,POINTS(np)%Rdata(1:3),void_2, POINTS(np)%Rdata(1:3))
#if I_PRINT >= 2
              write(*,*)'upgrade2double: Newton method applied.'
#endif
! ............delete extra surface
              SURFACES(surfs(3))%Type = 'Void'
              deallocate(SURFACES(surfs(3))%Rdata)
              NRSURFS = NRSURFS - 1
! ..........end treat as a 2 surf point
            endif
! ......end select number of surfaces
        end select
! ......reset point type
        POINTS(np)%Type = 'Regular'
        deallocate(POINTS(np)%Idata)
#if I_PRINT >= 2
        write(*,*)'upgrade2double: resetting point type.'
#endif
! ....end of loop through vertices
      enddo
! ..end of loop through triangles
    enddo
!
! ..reset GEOM_TOL to original value
    GEOM_TOL = geom_tol_orig
#if I_PRINT >= 1
    write(*,*)'upgrade2double: done!'
#endif
!
end subroutine upgrade2double_NEW


!-----------------------------------------------------------------
integer function factorial(n)
! factorial of n
  integer :: i,n
!
  if (n .eq. 0) then
    factorial = 1
    return
  endif
  factorial = 1
  do i = 2, n
    factorial = factorial*i
  enddo
end function factorial


integer function combinations(n, m)
! how many ways can we choose m objects from n objects?
  integer :: factorial
!
  if (n .lt. m ) then
    write(*,*)'combinations: n should be greater than m; n, m = ',n,m
    stop
  endif
  combinations = factorial(n)/(factorial(m)*factorial(n - m))
end function combinations
