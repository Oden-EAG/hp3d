!--------------------------------------------------------------------------------------
subroutine allocate_NEW_POINT
!--------------------------------------------------------------------------------------
! LATEST REVISION: Jul 09                                                             |
!                                                                                     |
! PURPOSE: routine classifies points into the following groups:                       |
!                                                                                     |
!    -3 = INT       points in the interior of the split surface                       |
!    -2 = POS_SIDE  points on the positive side of the split surface                  |
!    -1 = NEG_SIDE  points on the negative side of the split surface                  |
!     0 = INT^c     points not in the interior of split surface, i.e. complement      |
!                   of INT                                                            |
!--------------------------------------------------------------------------------------
  use SPLIT_SURF
  use GMP
  use control
!--------------------------------------------------------------------------------------
  IMPLICIT NONE
!--------------------------------------------------------------------------------------
! LOCAL VARIABLES
  real*8               :: fval
  real*8, dimension(3) :: xp,dfdx
  real*8               :: pos_shift, neg_shift
  integer                :: np,np_pos,np_neg
  integer                :: status
!--------------------------------------------------------------------------------------
! FUNCTIONS
  integer                :: if_bound
!--------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'allocate_NEW_POINT: allocating NEW_POINT.'
#endif
! ..allocate new_point and initialize it to 0
    allocate(new_point(MAXNP), STAT = status)
    if (status .ne. 0) then
      write(*,*)'allocate_NEW_POINT: new_point not allocated.'
      stop
    endif
    new_point(1:MAXNP) = 0
! ..set maximum allowed positive and negative shifts to big values
    pos_shift = 1.d8;  neg_shift = 1.d8
    np_pos = 0;  np_neg = 0
! ..loop through points
    do np = 1, NRPOINT
      xp(1:3) = POINTS(np)%Rdata(1:3)
      call surf(Nsplit_MOD,xp, fval,dfdx)
!
! ******************************  CLASSIFY POINT  *********************************** |
! ....point is on the split surface
      if (abs(fval) .le. GEOM_TOL) then
        new_point(np) = -3
! ......point is on the split surface, but not inside INT, i.e. OUTSIDE or BOUDARY
        if (if_bound(xp,Nr_bound_MOD,Ns_bound_MOD) .le. 0)  new_point(np) = 0
! ....point is on NEG_SIDE
      elseif (fval .lt. -GEOM_TOL) then
! ......point is inside the bounding surfaces
        if (if_bound(xp,Nr_bound_MOD,Ns_bound_MOD) .eq. 1) then
          if (abs(fval) .lt. neg_shift) then
            np_neg = np
            neg_shift = abs(fval)
#if I_PRINT >= 2
            write(*,*)'allocate_NEW_POINT: neg_shift = ',neg_shift
#endif
          endif
        endif
        new_point(np) = -1
! ....point is on POS_SIDE
      elseif (fval .gt. GEOM_TOL) then
! ......point is inside the bounding surfaces
        if (if_bound(xp,Nr_bound_MOD,Ns_bound_MOD) .eq. 1) then
          if (abs(fval) .lt. pos_shift) then
            np_pos = np
            pos_shift = abs(fval)
#if I_PRINT >= 2
            write(*,*)'allocate_NEW_POINT: pos_shift = ',pos_shift
#endif
          endif
        endif
        new_point(np) = -2
      endif
! ..end of loop through points
    enddo
! ..check positive and negative shifts
    status = 0
    if (pos_shift .lt. Dh2_MOD) then
      status = 1
      write(*,*)'allocate_NEW_POINT: max allowed value Dh2 = ',pos_shift
      write(*,*)'                                      np  = ',np_pos
      call print_GMP
    endif
    if (neg_shift .lt. Dh1_MOD) then
      status = 1
      write(*,*)'allocate_NEW_POINT: max allowed value Dh1 = ',neg_shift
      write(*,*)'                                      np  = ',np_neg
      call print_GMP
    endif
    if (status .eq. 1)  stop
!
#if I_PRINT >= 2
    write(*,*)'allocate_NEW_POINT: max allowed value Dh1 = ',neg_shift
    write(*,*)'                                      np  = ',np_neg
    call print_GMP
    write(*,*)'allocate_NEW_POINT: max allowed value Dh2 = ',pos_shift
    write(*,*)'                                      np  = ',np_pos
    call print_GMP
#endif
!
#if I_PRINT >= 1
    write(*,*)'allocate_NEW_POINT: done!'
#endif
!
end subroutine allocate_NEW_POINT
!--------------------------------------------------------------------------------------
