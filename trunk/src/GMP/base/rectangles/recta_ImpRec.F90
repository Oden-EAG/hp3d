!-----------------------------------------------------------------------
!> Purpose : routine evaluates physical coordinates and derivatives of 
!!           of an implicit rectangle
!!
!! @param[in ] No     - rectangle number
!! @param[in ] Eta    - reference coordinates of a point
!! @param[out] X      - physical coordinates of the point
!! @param[out] Dxdeta - derivatives of the physical coordinates
!!      
!! @revision Nov 12
!-----------------------------------------------------------------------
!  System of equations:
!
!                          f_1 = 0
!  (1 - eta_1)*f_5 + eta_1*f_3 = 0
!  (1 - eta_2)*f_2 + eta_2*f_4 = 0
!
!  Stretching factors s = s_iedge(eta), such that s_iedge(0)=0
!-----------------------------------------------------------------------
subroutine recta_ImpRec(No,Eta, X,Dxdeta)
!
      use GMP          , only : RECTANGLES
      use element_data , only : QUADR_EDGE_TO_VERT
      use control      , only : GEOM_TOL
      implicit none
      integer              ,intent(in ) :: No
      real*8,dimension(2  ),intent(in ) :: Eta
      real*8,dimension(3  ),intent(out) :: X
      real*8,dimension(3,2),intent(out) :: Dxdeta
!      
      real*8, dimension(  4) :: shapH
      real*8, dimension(2,4) :: dshapH
      real*8, dimension(3,4) :: xv
      real*8, dimension(3,3) :: aux,adx1,adx2
      real*8, dimension(3  ) :: der,der1,der2,xs,xmid,aux1,aux2,xc,dxcdt,de
      real*8, dimension(4  ) :: sfact,fxi,dfxi
      real*8, dimension(2  ) :: dtdeta
      integer,dimension(4  ) :: ncurv,npoint
      integer,dimension(5  ) :: nsurf

      real*8  :: fval,fval1,fval2,fval3,fval4,fval5,t
      integer :: i,j,iv1,iv2,norient,ifl1,ifl2,ii
      integer :: iprint
!-----------------------------------------------------------------------
!
      iprint=0
      if ((RECTANGLES(No)%Type.ne.'ImpRec')) then
        write(*,7000) RECTANGLES(No)%Type
 7000   format(' recta_ImpRec: WRONG CALL! Type = ',a10)
        stop
      endif
!
 5    continue
!
      if (iprint.eq.1) then
        write(*,7001) No,Eta
 7001   format(' recta_ImpRec: No,Eta = ',i4,2x,2f8.3)
      endif
!
!  ...vertex shape functions
      shapH(1)=(1.d0-Eta(1))*(1.d0-Eta(2))
      shapH(2)=(     Eta(1))*(1.d0-Eta(2))
      shapH(3)=(     Eta(1))*(     Eta(2))
      shapH(4)=(1.d0-Eta(1))*(     Eta(2))
!
      dshapH(1,1)=-(1.d0-Eta(2))
      dshapH(2,1)=-(1.d0-Eta(1))
      dshapH(1,2)= (1.d0-Eta(2))
      dshapH(2,2)=-Eta(1)
      dshapH(1,3)= Eta(2)
      dshapH(2,3)= Eta(1)
      dshapH(1,4)=-Eta(2)
      dshapH(2,4)= (1.d0-Eta(1))
!
!  ...curves, vertex numbers and coordinates
      do i=1,4
        ncurv( i)=RECTANGLES(No)%EdgeNo(i)
        npoint(i)=RECTANGLES(No)%VertNo(i)
        call pointr(npoint(i), xv(1,i))
      enddo
!
!  ...evaluate coefficients defining the bilinear map
      do j=1,3
        aux(j,1)=xv(j,2)-xv(j,1)
        aux(j,2)=xv(j,4)-xv(j,1)
        aux(j,3)=xv(j,1)-xv(j,2)+xv(j,3)-xv(j,4)
      enddo
!
!  ...evaluate starting point for NR iterations and the midpoint
      do j=1,3
        xs(  j)=xv(j,1)+Eta(1)*aux(j,1)+Eta(2)*aux(j,2)+Eta(1)*Eta(2)*aux(j,3)
        xmid(j)=xv(j,1)+0.5d0*aux(j,1)+0.5d0*aux(j,2)+0.25d0*aux(j,3)
      enddo
!
!  ...get surface numbers
      nsurf(1:5)=RECTANGLES(No)%Idata(1:5)
      if (iprint.eq.1) then
        write(*,7002) ncurv, npoint, nsurf
 7002   format(' recta_ImpRec: CURVES   = ',4i6, &
             /,'               POINTS   = ',4i6, &
             /,'               SURFACES = ',5i6) 
      endif
!
!  ...determine the sign factors to adjust orientation of surfaces
!
!                 
!       4        |        3
!        * * * * * * * * * 
!        *       4       *
!        *               *
!        *               *
!      5 * --    1     3 * --  
!        *               *
!        *               *
!        *       |       *
!        * * * * * * * * *
!       1        2        2
!
!      
      sfact(1:4) = 1.d0
      call surf(nsurf(2),xmid, fval,de)
      if (fval.lt.0) sfact(1) = -1.d0
      call surf(nsurf(3),xmid, fval,de)
      if (fval.gt.0) sfact(2) = -1.d0
      call surf(nsurf(4),xmid, fval,de)
      if (fval.gt.0) sfact(3) = -1.d0
      call surf(nsurf(5),xmid, fval,de)
      if (fval.lt.0) sfact(4) = -1.d0
      if (iprint.eq.1) then
        write(*,8003) (sfact(ii),ii=1,4)
 8003   format(' recta_ImpRec: sfact = ',4f8.3)
      endif
!
!  ...verify the compatibility of vertex data
      do i=1,4
        call surf(nsurf(1),xv(1:3,i), fval1,der1)
        if (abs(fval1).gt.GEOM_TOL) then
          write(*,7008) i,npoint(i),1,nsurf(1),fval1,No
 7008     format(' recta_ImpRec: INCOMPATIBILITY! POINT ',i2,'(',i4,')', &
                 ' DOES NOT LIE ON SURFACE ',i2,'(',i3,')',              &
                 ' VALUE = ',e12.5,' No = ',i4)
          call pause
        endif
        call surf(nsurf(1+i),xv(1:3,i), fval1,der1)
        if (abs(fval1).gt.GEOM_TOL) then
          write(*,7008) i,npoint(i),1+i,nsurf(1+i),fval1,No
          call pause
        endif
        j=mod(i+2,4)
        call surf(nsurf(2+j),xv(1:3,i), fval1,der1)
        if (abs(fval1).gt.GEOM_TOL) then
          write(*,7008) i,npoint(i),2+j,nsurf(2+j),fval1,No
          call pause
        endif
      enddo
!
!-----------------------------------------------------------------------
!  Step 1 : evaluate stretching functions and their derivatives wrt 
!           reference coordinates
!-----------------------------------------------------------------------
!
!  ...FIRST EDGE: SURFS 5 -> 3..........................................
      iv1=QUADR_EDGE_TO_VERT(1,1) ; iv2=QUADR_EDGE_TO_VERT(2,1)
      call proj_r2e(iv1,iv2,shapH,dshapH, t,dtdeta)
      norient=0 ; if (ncurv(1).lt.0) norient=1
      call curve_local(iabs(ncurv(1)),norient,t, xc,dxcdt)
!
!  ...check
      call surf(nsurf(1),xc, fval1,der1)
      call surf(nsurf(2),xc, fval2,der2)
      if ((fval1.gt.GEOM_TOL).or.(fval2.gt.GEOM_TOL)) then
        write(*,7007) iabs(ncurv(1)),nsurf(1),nsurf(2),No,fval1,fval2
 7007   format(' recta_ImpRec: INCONSISTENCY: POINT ON CURVE ',i5, &
               ' DOES NOT LIE ON SURFACES ',2i6,                   &
               ' DEFINING IMPLICIT RECTANGLE', i5,2e12.5)
        stop
      endif
!
!  ...surfaces bounding 1st edge curve
      call surf(nsurf(5),xc, fval1,der1)
      call surf(nsurf(3),xc, fval2,der2)
!  ...adjust orientations 
      fval1=fval1*sfact(4) ; der1(1:3)=der1(1:3)*sfact(4)
      fval2=fval2*sfact(2) ; der2(1:3)=der2(1:3)*sfact(2)
      if (iprint.eq.1) then
        write(*,8888) fval1,fval2
 8888   format(' recta_ImpRec: 1st EDGE, fval1,fval2 = ',2(e12.5,2x))
      endif
!
!  ...evaluate the value of the stretching function and its derivative
      fxi(1)  = fval1/(fval1-fval2)
      dfxi(1) = ((1.d0-fxi(1))*dot_product(der1,dxcdt)                &
                     + fxi(1) *dot_product(der2,dxcdt))/(fval1-fval2)
!
!  ...SECOND EDGE : SURFS 2 -> 4........................................
      iv1=QUADR_EDGE_TO_VERT(1,2) ; iv2=QUADR_EDGE_TO_VERT(2,2)
      call proj_r2e(iv1,iv2,shapH,dshapH, t,dtdeta)
      norient=0 ; if (ncurv(2).lt.0) norient=1
      call curve_local(iabs(ncurv(2)),norient,t, xc,dxcdt)
!
!  ...check
      call surf(nsurf(1),xc, fval1,der1)
      call surf(nsurf(3),xc, fval2,der2)
      if ((fval1.gt.GEOM_TOL).or.(fval2.gt.GEOM_TOL)) then
        write(*,7007) iabs(ncurv(2)),nsurf(1),nsurf(3),No,fval1,fval2
        stop
      endif
!
!  ...bounding surfaces      
      call surf(nsurf(2),xc, fval1,der1)
      call surf(nsurf(4),xc, fval2,der2)
!  ...adjust orientations 
      fval1=fval1*sfact(1) ; der1(1:3)=der1(1:3)*sfact(1)
      fval2=fval2*sfact(3) ; der2(1:3)=der2(1:3)*sfact(3)
      if (iprint.eq.1) then
        write(*,8889)fval1,fval2
 8889   format(' recta_ImpRec: 2nd EDGE, fval1,fval2 = ',2(e12.5,2x))
      endif
!
!  ...evaluate the value of the stretching function and its derivative
      fxi(2)  = fval1/(fval1-fval2)
      dfxi(2) = ((1.d0-fxi(2))*dot_product(der1,dxcdt)                &
                     + fxi(2) *dot_product(der2,dxcdt))/(fval1-fval2)
!
!
!  ...THIRD EDGE : SURFS 5 -> 3.........................................
      iv1=QUADR_EDGE_TO_VERT(1,3) ; iv2=QUADR_EDGE_TO_VERT(2,3)
      call proj_r2e(iv1,iv2,shapH,dshapH, t,dtdeta)
      norient=0 ; if (ncurv(3).lt.0) norient=1
      call curve_local(iabs(ncurv(3)),norient,t, xc,dxcdt)
!
!  ...check
      call surf(nsurf(1),xc, fval1,der1)
      call surf(nsurf(4),xc, fval2,der2)
      if ((fval1.gt.GEOM_TOL).or.(fval2.gt.GEOM_TOL)) then
        write(*,7007) iabs(ncurv(3)),nsurf(1),nsurf(4),No,fval1,fval2
        stop
      endif
!
!  ...bounding surfaces      
      call surf(nsurf(5),xc, fval1,der1)
      call surf(nsurf(3),xc, fval2,der2)
!  ...adjust orientations 
      fval1=fval1*sfact(4) ; der1(1:3)=der1(1:3)*sfact(4)
      fval2=fval2*sfact(2) ; der2(1:3)=der2(1:3)*sfact(2)
      if (iprint.eq.1) then
        write(*,8890)fval1,fval2
 8890   format(' recta_ImpRec: 3rd EDGE, fval1,fval2 = ',2(e12.5,2x))
      endif
!
      fxi(3) = fval1/(fval1-fval2)
      dfxi(3) = ((1.d0-fxi(3))*dot_product(der1,dxcdt)                & 
                      +fxi(3) *dot_product(der2,dxcdt))/(fval1-fval2)
!
!
!  ...FOURTH EDGE : SURFS 2 -> 4........................................
      iv1=QUADR_EDGE_TO_VERT(1,4) ; iv2=QUADR_EDGE_TO_VERT(2,4)
      call proj_r2e(iv1,iv2,shapH,dshapH, t,dtdeta)
      norient=0 ; if (ncurv(4).lt.0) norient=1
      call curve_local(iabs(ncurv(4)),norient,t, xc,dxcdt)
!
!  ...check consistency
      call surf(nsurf(1),xc, fval1,der1)
      call surf(nsurf(5),xc, fval2,der2)
      if ((fval1.gt.GEOM_TOL).or.(fval2.gt.GEOM_TOL)) then
        write(*,7007)iabs(ncurv(4)),nsurf(1),nsurf(5),No,fval1,fval2
        stop
      endif
!
!  ...bounding surfaces      
      call surf(nsurf(2),xc, fval1,der1)
      call surf(nsurf(4),xc, fval2,der2)
!  ...adjust orientations 
      fval1=fval1*sfact(1) ; der1(1:3)=der1(1:3)*sfact(1)
      fval2=fval2*sfact(3) ; der2(1:3)=der2(1:3)*sfact(3)
      if (iprint.eq.1) then
        write(*,8891)fval1,fval2
 8891   format(' recta_ImpRec: 4th EDGE, fval1,fval2 = ',2(e12.5,2x))
      endif
!
      fxi(4) = fval1/(fval1-fval2)
      dfxi(4) = ((1.d0-fxi(4))*dot_product(der1,dxcdt)          &
                     + fxi(4) *dot_product(der2,dxcdt))/(fval1-fval2)
!
!  ...printing
      if (iprint.eq.1) then
         write(*,*)'NCURV =',ncurv
         write(*,7003) ((xv(i,j), i=1,3), j=1,4)
 7003    format(' recta_ImpRec: POINTS =',4(3f8.3,2x))
         write(*,7004) xs
 7004    format(' recta_ImpRec: STARTING POINT',3f8.3)
         write(*,7005) Eta
 7005    format(' recta_ImpRec: Eta =',2f8.3)
         write(*,7006) (fxi(i),dfxi(i),i=1,4)
 7006    format(' recta_ImpRec: fxi,dfxi = ',4(2f8.3,3x))
      endif
!
!-----------------------------------------------------------------------
! STEP 2 : solve the linear system defining the point
!-----------------------------------------------------------------------
      call mnewt(4,nsurf,Eta,fxi,xs,sfact, X)
!
!-----------------------------------------------------------------------
! STEP 3 : compute derivatives
!-----------------------------------------------------------------------
!
!  ...FIRST EQUATION....................................................
      call surf(nsurf(1),X, fval1,der)
!
!  ...check consistency
      if (fval1.gt.GEOM_TOL) then
        write(*,8001) No, nsurf(1), fval1
 8001   format(' recta_ImpRec: RECT,SURFACE,fval1 = ',2i5,e12.5)
        stop
      endif
      aux1(1) = 0.d0
      aux2(1) = 0.d0
      do i=1,3
        adx1(1,i) = der(i)
      enddo
!
!  ...SECOND EQUATION...................................................
      call surf(nsurf(5),X, fval5,der1)
      call surf(nsurf(3),X, fval3,der2)
      fval5 = fval5*sfact(4); der1(1:3) = der1(1:3)*sfact(4)
      fval3 = fval3*sfact(2); der2(1:3) = der2(1:3)*sfact(2)
      do i=1,3
        adx1(2,i)=(1.d0-Eta(2))*((1.d0-fxi(1))*der1(i)+fxi(1)*der2(i)) &
                       +Eta(2) *((1.d0-fxi(3))*der1(i)+fxi(3)*der2(i))
      enddo
      aux1(2) = ((1.d0-Eta(2))*dfxi(1) + Eta(2)*dfxi(3))*(fval5-fval3)
      aux2(2) = (fval3-fval5)*(fxi(1)-fxi(3))
!
!  ...THIRD EQUATION....................................................
      call surf(nsurf(2),X, fval2,der1)
      call surf(nsurf(4),X, fval4,der2)
      fval2 = fval2*sfact(1); der1(1:3) = der1(1:3)*sfact(1)
      fval4 = fval4*sfact(3); der2(1:3) = der2(1:3)*sfact(3)
      do i=1,3
        adx1(3,i)=(1.d0-Eta(1))*((1.d0-fxi(4))*der1(i)+fxi(4)*der2(i)) &
                       +Eta(1) *((1.d0-fxi(2))*der1(i)+fxi(2)*der2(i))
      enddo
      aux2(3) = ((1.d0-Eta(1))*dfxi(4) + Eta(1)*dfxi(2))*(fval2-fval4)
      aux1(3) = (fval2-fval4)*(fxi(2)-fxi(4))
!
!  ...stiffness matrices are identical for both derivatives
      adx2(1:3,1:3)=adx1(1:3,1:3)
!
      if (iprint.eq.1) then
        write(*,*) 'recta_ImpRec: MATRICES FOR dx/dxi1 = '
        do i=1,3
          write(*,8002) aux1(i), (adx1(i,j),j=1,3)
 8002     format(e12.5,3x,3e12.5)
        enddo
        write(*,*) 'recta_ImpRec: MATRICES FOR dx/dxi2 = '
        do i=1,3
          write(*,8002) aux2(i), (adx2(i,j),j=1,3)
        enddo
        call pause
      endif
!
!  ...solve for the derivatives
      call saruss(adx1,aux1, Dxdeta(1,1),ifl1)
      call saruss(adx2,aux2, Dxdeta(1,2),ifl2)
      if (ifl1.ne.0.or.ifl2.ne.0) then
        iprint=1
        goto 5
      endif
!
!
endsubroutine recta_ImpRec
