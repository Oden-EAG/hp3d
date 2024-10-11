!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine returns coordinates of a point
!
!   arguments :
!     in:
!               No     - number of the point
!     out:
!               X      - coordinates of the point
!
!----------------------------------------------------------------------
!
   subroutine pointr(No, X)
!
      use GMP
      implicit none
!
      integer :: No
      real(8) :: X(NDIM)
!
!  ...surface numbers (only first three entries are meaningful),
!     starting point for NR iterations
      integer :: nsurf(6)
      real(8) :: xs(3)
!
!  ...work space
      real(8) :: void1(1),void2(4),void3(4)
!
!  ...point type
      character(len=10) :: type     
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!----------------------------------------------------------------------
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) No
 7001   format('pointr: No = ',i4)
      endif
#endif
!
!  ...get the point's type
      type = POINTS(No)%Type
!
      select case(type)
!
!  ...regular point, or point with a normal, or point on an algebraic
!     surface
      case('Regular','CoorNrm','SharpPt')
        X(1:NDIM) = POINTS(No)%Rdata(1:NDIM)
!
!  ...implicit  point
      case('Implicit')
!
!  .....get surface numbers
        nsurf(1:3) = POINTS(No)%Idata(1:3)
!
!  .....get starting point
        xs(1:3) = POINTS(No)%Rdata(1:3)
!
!  .....call Newton solver to determine the point coordinates
        call mnewt(1,nsurf,void1,void2,xs,void3, X)
!
      case default
        write(*,*) 'pointr: type = ',type
        stop
      end select
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) X
 7002   format('pointr: X = ',3f8.3)
      endif
#endif
!
   end subroutine pointr
