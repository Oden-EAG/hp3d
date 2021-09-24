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
!     arguments:
!        in:
!             Idec      = 1 change CURVE to POINT connectivity
!                       = 2 change TRIANGLE to POINT connectivity
!                       = 3 change RECTANGLE to POINT connectivity
!             No        - number of the object (CURVE,TRIANGLE or RECTANGLE)
!             Nor       - orientation to be modified
!                       = 0,1               for Idec=1
!                         0,1,2,3,4,5       for Idec=2
!                         0,1,2,3,4,5,6,7   for Idec=2
!                         
!        out:             changes in GMP module

!----------------------------------------------------------------------
!
      subroutine change_orientations(Idec,No,Nor)
!
      use GMP
!
      implicit none
!
      integer, intent(in)  :: Idec,No,Nor
!
      integer :: itemp(8)
!
!----------------------------------------------------------------------
!
      select case(Idec)
!
!  ...changing CURVE to POINT connectivity
      case(1)
        if ((No.lt.1).or.(No.gt.NRCURVE)) then
          write(*,7010) No,NRCURVE
 7010     format('change_orientations: No,NRCURVE = ',2i10)
          stop 1
        endif
        itemp(1:2) = CURVES(No)%EndPoNo(1:2)
        select case(Nor)
        case(0)
        case(1)
          CURVES(No)%EndPoNo(1) = itemp(2)
          CURVES(No)%EndPoNo(2) = itemp(1)
        case default
          write(*,7020) Nor
 7020     format('change_orientations: CURVE to POINT orientation = ',i10)
          stop 1
        end select
!
!  ...changing TRIANGLE to POINT connectivity
      case(2)
        if ((No.lt.1).or.(No.gt.NRTRIAN)) then
          write(*,7030) No,NRTRIAN
 7030     format('change_orientations: No,NRTRIAN = ',2i10)
          stop 1
        endif
        itemp(1:3) = TRIANGLES(No)%VertNo(1:3)
        select case(Nor)
        case(0)
        case(1)
          TRIANGLES(No)%VertNo(1) = itemp(2)
          TRIANGLES(No)%VertNo(2) = itemp(3)
          TRIANGLES(No)%VertNo(3) = itemp(1)
        case(2)
          TRIANGLES(No)%VertNo(1) = itemp(3)
          TRIANGLES(No)%VertNo(2) = itemp(1)
          TRIANGLES(No)%VertNo(3) = itemp(2)
        case(3)
          TRIANGLES(No)%VertNo(1) = itemp(1)
          TRIANGLES(No)%VertNo(2) = itemp(3)
          TRIANGLES(No)%VertNo(3) = itemp(2)
        case(4)
          TRIANGLES(No)%VertNo(1) = itemp(3)
          TRIANGLES(No)%VertNo(2) = itemp(2)
          TRIANGLES(No)%VertNo(3) = itemp(1)
        case(5)
          TRIANGLES(No)%VertNo(1) = itemp(2)
          TRIANGLES(No)%VertNo(2) = itemp(1)
          TRIANGLES(No)%VertNo(3) = itemp(3)
        case default
          write(*,7040) Nor
 7040     format('change_orientations: TRIANGLE to POINT orientation = ',i10)
          stop 1
        end select
!
!  ...changing RECTANGLE to POINT connectivity
      case(3)
        if ((No.lt.1).or.(No.gt.NRRECTA)) then
          write(*,7050) No,NRRECTA
 7050     format('change_orientations: No,NRRECTA = ',2i10)
          stop 1
        endif
        itemp(1:4) = RECTANGLES(No)%VertNo(1:4)
        select case(Nor)
        case(0)
        case(1)
          RECTANGLES(No)%VertNo(1) = itemp(2)
          RECTANGLES(No)%VertNo(2) = itemp(3)
          RECTANGLES(No)%VertNo(3) = itemp(4)
          RECTANGLES(No)%VertNo(4) = itemp(1)
        case(2)
          RECTANGLES(No)%VertNo(1) = itemp(3)
          RECTANGLES(No)%VertNo(2) = itemp(4)
          RECTANGLES(No)%VertNo(3) = itemp(1)
          RECTANGLES(No)%VertNo(4) = itemp(2)
        case(3)
          RECTANGLES(No)%VertNo(1) = itemp(4)
          RECTANGLES(No)%VertNo(2) = itemp(1)
          RECTANGLES(No)%VertNo(3) = itemp(2)
          RECTANGLES(No)%VertNo(4) = itemp(3)
        case(4)
          RECTANGLES(No)%VertNo(1) = itemp(1)
          RECTANGLES(No)%VertNo(2) = itemp(4)
          RECTANGLES(No)%VertNo(3) = itemp(3)
          RECTANGLES(No)%VertNo(4) = itemp(2)
        case(5)
          RECTANGLES(No)%VertNo(1) = itemp(4)
          RECTANGLES(No)%VertNo(2) = itemp(3)
          RECTANGLES(No)%VertNo(3) = itemp(2)
          RECTANGLES(No)%VertNo(4) = itemp(1)
        case(6)
          RECTANGLES(No)%VertNo(1) = itemp(3)
          RECTANGLES(No)%VertNo(2) = itemp(2)
          RECTANGLES(No)%VertNo(3) = itemp(1)
          RECTANGLES(No)%VertNo(4) = itemp(4)
        case(7)
          RECTANGLES(No)%VertNo(1) = itemp(2)
          RECTANGLES(No)%VertNo(2) = itemp(1)
          RECTANGLES(No)%VertNo(3) = itemp(4)
          RECTANGLES(No)%VertNo(4) = itemp(3)
        case default
          write(*,7060) Nor
 7060     format('change_orientations: RECTANGLE to POINT orientation = ',i10)
          stop 1
        end select
      end select
!
!
      end subroutine change_orientations
