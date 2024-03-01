!-----------------------------------------------------------------------
!> Purpose : routine evaluates physical coordinates and derivatives of
!!           a transfinite interpolation rectangle
!!
!! @param[in]  No     - rectangle number
!! @param[in]  Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives of the physical coordinates
!!
!! @revision Nov 12
!-----------------------------------------------------------------------
!
subroutine recta_TraQua(No,Eta, X,Dxdeta)
!
      use GMP , only : RECTANGLES , SURFACES
!
      implicit none
      integer               ,intent(in ) :: No
      real(8),dimension(2  ),intent(in ) :: Eta
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,2),intent(out) :: Dxdeta
!-----------------------------------------------------------------------
!  ...edges curves numbers and orientations
      integer,dimension(4)   :: noc,norientc
!  ...point on a curve and its derivative
      real(8),dimension(3)   :: xp,dxp
!  ...vertex and edge blending functions
      real(8),dimension(8)   :: val
!  ...derivatives of blending functions
      real(8),dimension(2,8) :: dval
!  ...edge coordinate
      real(8)                :: etac
!  ...derivative of edge coordinate
      real(8),dimension(2)   :: detac
!
      integer                :: k,i,j,np,ns
      integer                :: iprint
!-----------------------------------------------------------------------
!
      iprint=0
!
!  ...check consistency
      select case(RECTANGLES(No)%Type)
      case('TraQua')
      case('PTIRec')
        ns=RECTANGLES(No)%Idata(1)
        if (SURFACES(ns)%Type.ne.'VecPt') then
          write(*,7000) No,RECTANGLES(No)%Type
 7000     format(' recta_TraQua: wrong type! No,Type = ',i7,2x,a12)
          stop
        endif
      case default
        write(*,7000) No,RECTANGLES(No)%Type
        stop
      endselect
!
!  ...get the edge curves and orientations
      noc(1:4)=abs(RECTANGLES(No)%EdgeNo(1:4))
      norientc(1:4)=0
      if (RECTANGLES(No)%EdgeNo(1).lt.0) norientc(1)=1
      if (RECTANGLES(No)%EdgeNo(2).lt.0) norientc(2)=1
      if (RECTANGLES(No)%EdgeNo(3).lt.0) norientc(3)=1
      if (RECTANGLES(No)%EdgeNo(4).lt.0) norientc(4)=1
!
!  ...printing
      if (iprint.eq.1) then
        write(*,7001) No,Eta(1:2),noc(1:4),norientc(1:4)
 7001   format(' recta_TraQua: No,Eta,noc,norientc = ',i7,2(e12.5,2x),2x,4i4,2x,4i2)
      endif
!
!  ...calculate the blending functions
      call recta_blend(Eta, val,dval)
!
!  ...initialize
      X(1:3)=0.d0 ; Dxdeta(1:3,1:2)=0.d0 ; k=0
!
!-----------------------------------------------------------------------
!  STEP 1 : subtract vertex interpolant
!-----------------------------------------------------------------------
!  ...loop over vertices
      do i=1,4
        k=k+1
        np=RECTANGLES(No)%VertNo(i)
        call pointr(np, xp)
        X(1:3) = X(1:3) - xp(1:3)*val(k)
        do j=1,2
          Dxdeta(1:3,j) = Dxdeta(1:3,j) - xp(1:3)*dval(j,k)
        enddo
      enddo
!
!  ...printing
      if (iprint.eq.1) then
        write(*,7002)
 7002   format(' recta_TraQua: AFTER VERTEX CONTRIBUTIONS')
        write(*,7003) X(1:3)
 7003   format(' X = ',3(e12.5,2x))
        do i=1,2
          write(*,7004) i,Dxdeta(1:3,i)
 7004   format(' i,Dxdeta(:,i) = ',i1,2x,3(e12.5,2x))
        enddo
      endif
!
!-----------------------------------------------------------------------
!  STEP 2 : add edge bubbles
!-----------------------------------------------------------------------
!  ...loop over edges
      do i=1,4
        k=k+1
!
!  .....local curve coordinate
        detac(1:2)=0.d0
        select case(i)
        case(1,3) ; etac=Eta(1) ; detac(1)=1.d0
        case(2,4) ; etac=Eta(2) ; detac(2)=1.d0
        endselect
!
!  .....accumulate
        call curve_local(noc(i),norientc(i),etac, xp,dxp)
        X(1:3) = X(1:3) + xp(1:3)*val(k)
        do j=1,2
          Dxdeta(1:3,j) = Dxdeta(1:3,j)                    &
                           + dxp(1:3)*detac(j)*val(k)      &
                           + xp(1:3)*dval(j,k)
        enddo
!
!  .....printing
        if (iprint.eq.1) then
          write(*,6999) i
 6999     format(' recta_TraQua: AFTER EDGE ',i1)
          write(*,7003) X(1:3)
          do j=1,2
            write(*,7004) j,Dxdeta(1:3,j)
          enddo
          call pause
        endif
      enddo
!
!
end subroutine recta_TraQua
!
!
!
!-----------------------------------------------------------------------
!
!   routine name       - recta_blend
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine evaluates bilinear and linear
!                        blending functions for the transfinite
!                        interpolation on master rectangle and their
!                        derivatives
!
!   arguments :
!     in:
!               Eta    - reference coordinates of a point in the
!                        reference rectangle
!     out:
!               Val    - value of the blending functions
!               Dval   - derivatives
!
!-----------------------------------------------------------------------
!
   subroutine recta_blend(Eta, Val,Dval)
!
      implicit none
!
      real(8) :: Eta(2),Val(8),Dval(2,8)
!
      integer :: k
!
!-----------------------------------------------------------------------
!
!  ...vertex bilinear blending functions...
!
!  ...first vertex bilinear function
      k=1
      Val(k) = (1.d0-Eta(1))*(1.d0-Eta(2))
      Dval(1,k) = - (1.d0-Eta(2))
      Dval(2,k) = - (1.d0-Eta(1))
!
!  ...second vertex bilinear function
      k=2
      Val(k) = Eta(1)*(1.d0-Eta(2))
      Dval(1,k) =   (1.d0-Eta(2))
      Dval(2,k) = - Eta(1)
!
!  ...third vertex bilinear function
      k=3
      Val(k) = Eta(1)*Eta(2)
      Dval(1,k) =   Eta(2)
      Dval(2,k) =   Eta(1)
!
!  ...fourth vertex bilinear function
      k=4
      Val(k) = (1.d0-Eta(1))* Eta(2)
      Dval(1,k) = - Eta(2)
      Dval(2,k) =   1.d0-Eta(1)
!
!-----------------------------------------------------------------------
!
!  ...edge linear blending functions
!
!  ...first edge linear function
      k=5
      Val(k) = 1.d0-Eta(2)
      Dval(1,k) =  0.d0
      Dval(2,k) = -1.d0
!
!  ...second  edge linear function
      k=6
      Val(k) = Eta(1)
      Dval(1,k) = 1.d0
      Dval(2,k) = 0.d0
!
!  ...third edge linear function
      k=7
      Val(k) = Eta(2)
      Dval(1,k) = 0.d0
      Dval(2,k) = 1.d0
!
!  ...fourth edge linear function
      k=8
      Val(k) = 1.d0-Eta(1)
      Dval(1,k) = -1.d0
      Dval(2,k) =  0.d0
!
!
   end subroutine recta_blend
