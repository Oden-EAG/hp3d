c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine returns coordinates of a point
c
c   arguments :
c     in:
c               No     - number of the point
c     out:
c               X      - coordinates of the point
c
c----------------------------------------------------------------------
c
      subroutine pointr(No, X)
c
      use GMP
      implicit none
c
      integer :: No
      real(8) :: X(NDIM)
c
c  ...surface numbers (only first three entries are meaningful),
c     starting point for NR iterations
      integer :: nsurf(6)
      real(8) :: xs(3)
c
c  ...work space
      real(8) :: void1,void2(4),void3(4)
c
c  ...point type
      character(len=10) :: type     
c
      integer :: iprint=0
      if (iprint.eq.1) then
        write(*,7001) No
 7001   format('pointr: No = ',i4)
      endif
c
c  ...get the point's type
      type = POINTS(No)%Type
c
      select case(type)
c
c  ...regular point, or point with a normal, or point on an algebraic
c     surface
      case('Regular','CoorNrm','SharpPt')
        X(1:NDIM) = POINTS(No)%Rdata(1:NDIM)
c
c  ...implicit  point
      case('Implicit')
c
c  .....get surface numbers
        nsurf(1:3) = POINTS(No)%Idata(1:3)
c
c  .....get starting point
        xs(1:3) = POINTS(No)%Rdata(1:3)
c
c  .....call Newton solver to determine the point coordinates
        call mnewt(1,nsurf,void1,void2,xs,void3, X)
c
      case default
        write(*,*) 'pointr: type = ',type
        stop
      end select
      if (iprint.eq.1) then
        write(*,7002) X
 7002   format('pointr: X = ',3f8.3)
      endif
c
      end subroutine pointr
