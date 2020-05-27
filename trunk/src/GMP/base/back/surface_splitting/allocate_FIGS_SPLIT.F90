!----------------------------------------------------------------------------------
subroutine allocate_FIGS_SPLIT
!----------------------------------------------------------------------------------
! LATEST REVISION: Jul 09
!
! PURPOSE: routine identifies triangles and rectangles to be split.
!
! REMARKS: first iteration is needed to determine number of figures to split,
!          second interation fills up figs_split
!----------------------------------------------------------------------------------
  use GMP
  use SPLIT_SURF
  use control
!----------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------
! VARIABLES
  integer :: i,iter,ifig,nt,ion,iout,iv,np,nr
  integer :: status
!----------------------------------------------------------------------------------
! ..printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'allocate_FIGS_SPLIT: allocating figures to be splitted...'
#endif
!
! ..loop through iterations
    do iter = 1, 2
! ....set figure counter to 0
      ifig = 0
!
! **************************  T R I A N G L E S  ******************************** |
! ....loop through triangles
      do nt = 1, NRTRIAN
!==================================================================================
!  ion:  number of vertices on INT                                                |
!  iout: number of vertices on POS/NEG_SIDE                                       |
!==================================================================================
        ion = 0;  iout = 0
! ......loop through vertices
        do iv = 1, 3
          np = TRIANGLES(nt)%VertNo(iv)
! ........vertex on INT
          if (new_point(np) .eq. -3) then
          ion = ion + 1
! ........vertex on POS/NEG_SIDE
          elseif (new_point(np) .lt. 0) then
            iout = iout + 1
          endif
        enddo
!=================================================================================
!  IF (at least one vertex is in INT) .AND. (no vertex is on POS/NEG_SIDE) THEN  |
!    figure needs to be split                                                    |
!=================================================================================
        if ((ion .ge. 1) .and. (iout .eq. 0)) then
! ........increment figure counter
          ifig = ifig + 1
! ........if 2nd iteration
          if (iter .eq. 2) then
! ..........store figure type (1 = trian,  2 = recta)
            figs_split(1,ifig) = 1
! ..........store triangle number
            figs_split(2,ifig) = nt
! ..........check if triangle is on the split surface
            if (.not. associated(TRIANGLES(nt)%Idata)) then
              write(*,3) nt
              stop
            elseif (TRIANGLES(nt)%Idata(1) .ne. Nsplit_MOD) then
              write(*,3) nt
3             format(' allocate_FIGS_SPLIT: nt = ',i5,' is not on the split surface.')
              stop
            endif
! ........end if 2nd iteration
          endif
        endif
! ....end of loop through triangles
      enddo
!
! ***************************  R E C T A N G L E S  **************************** |
! ....loop through rectangles
      do nr = 1, NRRECTA
        ion = 0;  iout = 0
! ......loop through vertices
        do iv = 1, 4
          np = RECTANGLES(nr)%VertNo(iv)
          if (new_point(np) .eq. -3) then
            ion = ion + 1
          elseif (new_point(np) .lt. 0) then
            iout = iout + 1
          endif
! ......end of loop through vertices
        enddo
! ......if rectangle needs to be split
        if ((ion .eq. 2) .and. (iout .eq. 0)) then
          ifig = ifig + 1
! ........if 2nd interation store figure type and number
          if (iter .eq. 2) then
            figs_split(1,ifig) = 2
            figs_split(2,ifig) = nr
! ..........double check that rectangle is on the split surface
            if (.not. associated(RECTANGLES(nr)%Idata)) then
              write(*,4) nr
              stop
            elseif (RECTANGLES(nr)%Idata(1) .ne. Nsplit_MOD) then
              write(*,4) nr
4             format(' allocate_FIGS_SPLIT: nr = ',i5,' is not on the split surface.')
              stop
            endif
! ........end if 2nd iteration
          endif
! ......end rectangle needs to be split
        endif
! ....end of loop through rectangles
      enddo
! ....save number of figures to split
      nr_figs_to_split = ifig
      if (nr_figs_to_split .eq. 0) then
        write(*,*) 'allocate_FIGS_SPLIT: number of figures to split  = ',nr_figs_to_split
        stop
      endif
! ....allocate figs_to_split if 1st interations
      if (iter .eq. 1)  allocate(figs_split(2,nr_figs_to_split), STAT = status)
      if (status .ne. 0) then
        write(*,*)'allocate_FIGS_SPLIT: figs_split not allocated.'
        stop
#if I_PRINT >= 2
        write(*,*)'allocate_FIGS_SPLIT: figs_split allocated successfully.'
#endif
      endif
! ..end of loop through iterations
    enddo
#if I_PRINT >= 2
    write(*,1) nr_figs_to_split
1   format(' allocate_FIGS_SPLIT: nr_figs_to_split = ',i4)
    do i = 1, nr_figs_to_split
      if (figs_split(1,i) .eq. 1) then
        write(*,5) figs_split(2,i),TRIANGLES(figs_split(2,i))%VertNo(1), &
                                   TRIANGLES(figs_split(2,i))%VertNo(2), &
                                   TRIANGLES(figs_split(2,i))%VertNo(3)
5       format('  triangle: ',I8,';  vertices = ',3(I8))
      else
        write(*,6) figs_split(2,i),RECTANGLES(figs_split(2,i))%VertNo(1), &
                                   RECTANGLES(figs_split(2,i))%VertNo(2), &
                                   RECTANGLES(figs_split(2,i))%VertNo(3), &
                                   RECTANGLES(figs_split(2,i))%VertNo(4)
6       format(' rectangle: ',I8,';  vertices = ',4(I8))
      endif
    enddo
    call pause
#endif
#if I_PRINT >= 1
    write(*,*)'allocate_FIGS_TO_SPLIT: done!'
#endif
!
end subroutine allocate_FIGS_SPLIT
