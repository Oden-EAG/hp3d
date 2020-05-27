!---------------------------------------------------------------------------
subroutine duplicate_POINTS_AND_CURVES
!---------------------------------------------------------------------------
! LATEST REVISION: Jul 09
!
! PURPOSE: routine generates twin points and curves connecting twin points.
!
! REMARK: classifications of points is changed!
!---------------------------------------------------------------------------
  use SPLIT_SURF
  use U2D
!---------------------------------------------------------------------------
! VARIABLES
  integer :: i,ifig,nt,nr,nvrt,iv,np
!---------------------------------------------------------------------------
! printing flag (0,1)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'duplicate_POINTS_AND_CURVES: generating twin points and connecting curves.'
#endif
! ..allocate extra 2 planes for reprojecting points
    call allocate_EXTRA_PLANES
! ..loop through figures to split
    do ifig = 1, nr_figs_to_split
      select case(figs_split(1,ifig))
        case(1)
          nt = figs_split(2,ifig);  nvrt = 3
        case(2)
          nr = figs_split(2,ifig);  nvrt = 4
      end select
! ....loop through figure vertex points
      do iv = 1, nvrt
        select case(figs_split(1,ifig))
          case(1);  np = TRIANGLES(nt)%VertNo(iv)
          case(2);  np = RECTANGLES(nr)%VertNo(iv)
        end select
!
! *****************  CREATE TWIN POINT & CONN CURVE  ********************* |
!===========================================================================
!  REMARK: new_point(np) = -3  INT                                         |
!                          -2  POS_SIDE                                    |
!                          -1  NEG_SIDE                                    |
!                           0  INT^c                                       |
!===========================================================================
! ......if 1st visit to point
        if (new_point(np) .le. 0) then
          call create_TW_POINT_AND_CONN_CURVE(np)
        endif
!===========================================================================
!  REMARK: new_point(np) = twin_point  INT or INT^c = SPLIT SURFACE        |
!                          -2          POS_SIDE                            |
!                          -1          NEG_SIDE                            |
!===========================================================================
! .....end loop through figure vertices
      enddo
! ..end of loop through figures to split
    enddo
! ..deallocate extra 2 planes
    call deallocate_EXTRA_PLANES
#if I_PRINT >= 1
    write(*,*)'duplicate_POINTS_AND_CURVES: done!'
#endif
!
end subroutine duplicate_POINTS_AND_CURVES
