!-----------------------------------------------------------------------
!> @brief physical coordinates for a rectangle parametrization, and
!!           their derivatives wrt to a GIVEN system of coordinates
!!
!> @param[in ] No      - rectangle number
!> @param[in ] Norient - orientation of GLOBAL REFERENCE coordinates
!!                       wrt a GIVEN system of coordinates
!> @param[in ] T       - coordinates of a point in the rectangle
!!                       according to the GIVEN system
!> @param[out] X       - physical coordinates
!> @param[out] Dxdt    - derivatives of physical coordinates
!!
!> @date Nov 12
!-----------------------------------------------------------------------
subroutine recta_local(No,Norient,T, X,Dxdt)
!
      use node_types, only: QUAD
      implicit none
      integer,               intent(in ) :: No,Norient
      real(8),dimension(2  ),intent(in ) :: T
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,2),intent(out) :: Dxdt
!
      real(8),dimension(2  ) :: eta
      real(8),dimension(2,2) :: detadt
      real(8),dimension(3,2) :: dxdeta
      integer :: i
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!-----------------------------------------------------------------------
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7000) No,T,Norient
 7000   format(' recta_local: No,T,Norient = ',i6,2x,2(e12.5,2x,i1))
      endif
#endif
!
!  ...GIVEN -> GLOBAL REFERENCE
      call local2global(QUAD,T,Norient, eta,detadt)
!
!  ...GLOBAL REFERENCE -> PHYSICAL
      call recta(No,eta(1:2), X(1:3),dxdeta(1:3,1:2))
!
!  ...derivatives of PHYSICAL wrt GIVEN
      do i=1,2
        Dxdt(1:3,i) = dxdeta(1:3,1)*detadt(1,i) +  &
                      dxdeta(1:3,2)*detadt(2,i)
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        do i=1,3
          write(*,7001) i,X(i),Dxdt(i,1:2)
 7001     format(' i = ',i1,'; X(i),Dxdt(i,1:2) = ',3(e12.5,2x))
        enddo
      endif
#endif
!
end subroutine recta_local
!
!
!--------------------------------------------------------------------------
!> @brief physical coordinates and derivatives wrt to reference
!!           coordinates for a point in a GMP rectangle
!!
!> @param[in ]  No     - the rectangle number
!> @param[in ]  Eta    - reference coordinates of a point
!> @param[out]  X      - physical coordinates of the point
!> @param[out]  Dxdeta - derivatives of the physical coordinates
!!
!> @date Nov 12
!--------------------------------------------------------------------------
!
subroutine recta(No,Eta, X,Dxdeta)
!
      use GMP , only : RECTANGLES , CURVES
      integer               ,intent(in ) :: No
      real(8),dimension(2  ),intent(in ) :: Eta
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,2),intent(out) :: Dxdeta
      real(8),dimension(3,3)             :: dxdeta_aux
!--------------------------------------------------------------------------
!  ...vertex coordinates
      real(8),dimension(3,4) :: xv
!  ...workspace
      real(8),dimension(3,3) :: aux
!
      integer :: i,np,nc,j
!
!--------------------------------------------------------------------------
!
      select case(RECTANGLES(No)%Type)
!
!  ...Bilinear Rectangle...................................................
      case('BilQua')
!
!  .....check curve compatibility
        do i=1,4
          nc=abs(RECTANGLES(No)%EdgeNo(i))
          if (CURVES(nc)%Type.ne.'Seglin') then
            write(*,*)'recta: INCOMPATIBLE curve definition!'
            write(*,7002) No,i,nc,CURVES(nc)%Type
 7002       format(' No,i,nc = ',i6,2x,i1,2x,i6,'; type = ',a10)
            stop
          endif
        enddo
!
!  .....evaluate vertex coordinates
        do i=1,4
          np=RECTANGLES(No)%VertNo(i)
          call pointr(np, xv(1:3,i))
        enddo
!
!  .....loop over components
        do j=1,3
!  .......compute coefficients defining the bilinear map
          aux(j,1) = xv(j,2)-xv(j,1)
          aux(j,2) = xv(j,4)-xv(j,1)
          aux(j,3) = xv(j,1)-xv(j,2)+xv(j,3)-xv(j,4)
        enddo
!
!  .....evaluate physical coordinates and their derivatives
        do j=1,3
          Dxdeta(j,1) = aux(j,1) + Eta(2)*aux(j,3)
          Dxdeta(j,2) = aux(j,2) + Eta(1)*aux(j,3)
          X(j) = xv(j,1) + Eta(1)*aux(j,1) + Eta(2)*aux(j,2) + Eta(1)*Eta(2)*aux(j,3)
        enddo
!
!  ...Transfinite Interpolation..........................................
      case('TraQua') ; call recta_TraQua(No,Eta, X,Dxdeta)
!
!  ...Parametric Transfinite Interpolation...............................
      case('PTIRec') ; call recta_PTIRec(No,Eta, X,Dxdeta)
!
!  ...implicit rectangle.................................................
      case('ImpRec') ; call recta_ImpRec(No,Eta, X,Dxdeta)
!
!  ...conforming rectangle (LEGACY)......................................
      !case('ConfRec') ; call recta_ConfRec(No,Eta, X,Dxdeta)
!
!  ...C^1 transfinite interpolation/geometry reconstruction rectangle (GOOD LUCK!)
      !case('HermRec') ; call recta_HermRec(No,Eta, X,Dxdeta)
!
!  ...cylindrical coordinates rectangle..................................
      case('CylRec')
        call recta_CylRec(No,Eta, X,dxdeta_aux)
        Dxdeta(1:3,1:2) = dxdeta_aux(1:3,1:2)
!
      case default
        write(*,7003) RECTANGLES(No)%Type
   7003 format(' recta: unknown Type = ',a10)
      stop
!
      end select
!
!
end subroutine recta
