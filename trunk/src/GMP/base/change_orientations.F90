!----------------------------------------------------------------------
!
!     routine name      - change_orientations
!
!----------------------------------------------------------------------
!
!     latest revision:  - Jul 21
!
!     purpose:          - test implementation of orientations in shape
!                         functions and constrained approximation
!                         the routine should be called from input_DEFAULT
!                         just BERFORE the call to 'connect'
!
!     arguments (through module test_param):
!        in:
!          ORIENT_DEC   = 0 RETURN
!                       = 1 change CURVE to POINT connectivity
!                       = 2 change TRIANGLE to POINT connectivity
!                       = 3 change RECTANGLE to POINT connectivity
!          ORIENT_NO    - number of the object (CURVE,TRIANGLE or RECTANGLE)
!          ORIENT_OR    - orientation to be modified
!                       = 0,1               for ORIENT_DEC=1
!                         0,1,2,3,4,5       for ORIENT_DEC=2
!                         0,1,2,3,4,5,6,7   for ORIENT_DEC=2
!                         
!        out:             changes in GMP module

!----------------------------------------------------------------------
!
      subroutine change_orientations
!
      use GMP
      use test_param
!
      implicit none
!
      integer :: itemp(8)
!
!----------------------------------------------------------------------
!
      select case(ORIENT_DEC)
!
!  ...changing CURVE to POINT connectivity
      case(0)
        return
      case(1)
        if ((ORIENT_NO.lt.1).or.(ORIENT_NO.gt.NRCURVE)) then
          write(*,7010) ORIENT_NO,NRCURVE
 7010     format('change_orientations: ORIENT_NO,NRCURVE = ',2i10)
          stop 1
        endif
        itemp(1:2) = CURVES(ORIENT_NO)%EndPoNo(1:2)
        select case(ORIENT_OR)
        case(0)
        case(1)
          CURVES(ORIENT_NO)%EndPoNo(1) = itemp(2)
          CURVES(ORIENT_NO)%EndPoNo(2) = itemp(1)
        case default
          write(*,7020) ORIENT_OR
 7020     format('change_orientations: CURVE to POINT orientation = ',i10)
          stop 1
        end select
!
!  ...changing TRIANGLE to POINT connectivity
      case(2)
        if ((ORIENT_NO.lt.1).or.(ORIENT_NO.gt.NRTRIAN)) then
          write(*,7030) ORIENT_NO,NRTRIAN
 7030     format('change_orientations: ORIENT_NO,NRTRIAN = ',2i10)
          stop 1
        endif
        itemp(1:3) = TRIANGLES(ORIENT_NO)%VertNo(1:3)
        select case(ORIENT_OR)
        case(0)
        case(1)
          TRIANGLES(ORIENT_NO)%VertNo(1) = itemp(2)
          TRIANGLES(ORIENT_NO)%VertNo(2) = itemp(3)
          TRIANGLES(ORIENT_NO)%VertNo(3) = itemp(1)
        case(2)
          TRIANGLES(ORIENT_NO)%VertNo(1) = itemp(3)
          TRIANGLES(ORIENT_NO)%VertNo(2) = itemp(1)
          TRIANGLES(ORIENT_NO)%VertNo(3) = itemp(2)
        case(3)
          TRIANGLES(ORIENT_NO)%VertNo(1) = itemp(1)
          TRIANGLES(ORIENT_NO)%VertNo(2) = itemp(3)
          TRIANGLES(ORIENT_NO)%VertNo(3) = itemp(2)
        case(4)
          TRIANGLES(ORIENT_NO)%VertNo(1) = itemp(3)
          TRIANGLES(ORIENT_NO)%VertNo(2) = itemp(2)
          TRIANGLES(ORIENT_NO)%VertNo(3) = itemp(1)
        case(5)
          TRIANGLES(ORIENT_NO)%VertNo(1) = itemp(2)
          TRIANGLES(ORIENT_NO)%VertNo(2) = itemp(1)
          TRIANGLES(ORIENT_NO)%VertNo(3) = itemp(3)
        case default
          write(*,7040) ORIENT_OR
 7040     format('change_orientations: TRIANGLE to POINT orientation = ',i10)
          stop 1
        end select
!
!  ...changing RECTANGLE to POINT connectivity
      case(3)
        if ((ORIENT_NO.lt.1).or.(ORIENT_NO.gt.NRRECTA)) then
          write(*,7050) ORIENT_NO,NRRECTA
 7050     format('change_orientations: ORIENT_NO,NRRECTA = ',2i10)
          stop 1
        endif
        itemp(1:4) = RECTANGLES(ORIENT_NO)%VertNo(1:4)
        select case(ORIENT_OR)
        case(0)
        case(1)
          RECTANGLES(ORIENT_NO)%VertNo(1) = itemp(2)
          RECTANGLES(ORIENT_NO)%VertNo(2) = itemp(3)
          RECTANGLES(ORIENT_NO)%VertNo(3) = itemp(4)
          RECTANGLES(ORIENT_NO)%VertNo(4) = itemp(1)
        case(2)
          RECTANGLES(ORIENT_NO)%VertNo(1) = itemp(3)
          RECTANGLES(ORIENT_NO)%VertNo(2) = itemp(4)
          RECTANGLES(ORIENT_NO)%VertNo(3) = itemp(1)
          RECTANGLES(ORIENT_NO)%VertNo(4) = itemp(2)
        case(3)
          RECTANGLES(ORIENT_NO)%VertNo(1) = itemp(4)
          RECTANGLES(ORIENT_NO)%VertNo(2) = itemp(1)
          RECTANGLES(ORIENT_NO)%VertNo(3) = itemp(2)
          RECTANGLES(ORIENT_NO)%VertNo(4) = itemp(3)
        case(4)
          RECTANGLES(ORIENT_NO)%VertNo(1) = itemp(1)
          RECTANGLES(ORIENT_NO)%VertNo(2) = itemp(4)
          RECTANGLES(ORIENT_NO)%VertNo(3) = itemp(3)
          RECTANGLES(ORIENT_NO)%VertNo(4) = itemp(2)
        case(5)
          RECTANGLES(ORIENT_NO)%VertNo(1) = itemp(4)
          RECTANGLES(ORIENT_NO)%VertNo(2) = itemp(3)
          RECTANGLES(ORIENT_NO)%VertNo(3) = itemp(2)
          RECTANGLES(ORIENT_NO)%VertNo(4) = itemp(1)
        case(6)
          RECTANGLES(ORIENT_NO)%VertNo(1) = itemp(3)
          RECTANGLES(ORIENT_NO)%VertNo(2) = itemp(2)
          RECTANGLES(ORIENT_NO)%VertNo(3) = itemp(1)
          RECTANGLES(ORIENT_NO)%VertNo(4) = itemp(4)
        case(7)
          RECTANGLES(ORIENT_NO)%VertNo(1) = itemp(2)
          RECTANGLES(ORIENT_NO)%VertNo(2) = itemp(1)
          RECTANGLES(ORIENT_NO)%VertNo(3) = itemp(4)
          RECTANGLES(ORIENT_NO)%VertNo(4) = itemp(3)
        case default
          write(*,7060) ORIENT_OR
 7060     format('change_orientations: RECTANGLE to POINT orientation = ',i10)
          stop 1
        end select
      case default
        write(*,7070) ORIENT_DEC
 7070   format('change_orientations: ORIENT_DEC = ',i10)
        stop 1
      end select
!
!
      end subroutine change_orientations
