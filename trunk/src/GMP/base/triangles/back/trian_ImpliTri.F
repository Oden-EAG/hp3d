c----------------------------------------------------------------------
c
c   routine name       - trian_ImpliTri
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for an implicit triangle
c
c   arguments :
c     in:
c               No     - the curve number
c               Xi     - a value of the curve parameter (between 0
c                        and 1)
c     out:
c               X      - physical coordinates of the point
c               Dxdxi  - derivatives of the physical coordinates wrt
c                        to the parameters
c
c----------------------------------------------------------------------
c
      subroutine trian_ImpliTri(No,Xi, X,Dxdxi)
c
      use element_data , only : TRIAN_EDGE_TO_VERT
      use GMP
      use control
c
      implicit none
c
      integer :: No
      real(8) :: Xi(2),X(3),Dxdxi(3,2)
c
      real(8) :: ctr50(3)
      common /sphepoint/ ctr50
c
      real(8) :: xs(3),xif(4),xifdxi(2,2),der1(3),der2(3),
     .           dxcdt(3),aux1(3),aux2(3),xc(3),aux(3,3),
     .           sfact(4),de(3),xmid(3),dtdeta(2),xv(3,3)
      integer :: nsutr(6),ncurv(3)
c
      real(8) :: shapH(3)
      real(8) :: dshapH(2,3)
c
      integer :: i,ifl1,ifl2,iv1,iv2,nc,norient,np
      real(8) :: fval,fval1,fval2,t
c
c     function dip(A,B)
      real(8) :: dip
c
      integer :: iprint=0
c
 5    continue
c
      shapH(1)=1.d0-Xi(1)-Xi(2) ; dshapH(1:2,1)=(/-1.d0,-1.d0/)
      shapH(2)=     Xi(1)       ; dshapH(1:2,2)=(/ 1.d0, 0.d0/)
      shapH(3)=           Xi(2) ; dshapH(1:2,3)=(/ 0.d0, 1.d0/)
c
      xif(1:4)=0.d0 ; nsutr(1:6)=0 ; sfact(1:4)=0.d0
      if(TRIANGLES(No)%Type.ne.'ImpliTri')then
         write(*,*)'IMPTRI: WRONG TYPE FOR IMPLICIT TRIANGLE'
      endif
c
c  ...get edge curves and vertices
      do i=1,3
        nc=TRIANGLES(No)%EdgeNo(i) ; ncurv(i)=nc
        np=TRIANGLES(No)%VertNo(i)
        call pointr(np, xv(1,i))
      enddo
c
c  ...store first vertex, compute centroid
      ctr50(1:3)=xv(1:3,1)
      xmid(1:3)=xv(1:3,1)+xv(1:3,2)+xv(1:3,3) ; xmid(1:3)=xmid(1:3)/3.d0
c
c  ...linear vertex interpolant (start point for NR iterations)
      Dxdxi(1:3,1) = xv(1:3,2)-xv(1:3,1)
      Dxdxi(1:3,2) = xv(1:3,3)-xv(1:3,1)
      xs(1:3) = xv(1:3,1) + Xi(1)*Dxdxi(1:3,1)+Xi(2)*Dxdxi(1:3,2)
c
c  ...if at first vertex, done
      if ((Xi(1)+Xi(2)).eq.0) then
        X(1:3)=xs(1:3) ; return
      endif
c
c  ...get number of surfaces
      nsutr(1:4)=TRIANGLES(No)%Idata(1:4)
c
c  ...sign factors to adjust orientation of surfaces
      sfact(1:4)=1.d0
      call surf(nsutr(2),xmid, fval,de) ; if (fval.lt.0) sfact(1)=-1.d0
      call surf(nsutr(3),xmid, fval,de) ; if (fval.lt.0) sfact(2)= 1.d0
      call surf(nsutr(4),xmid, fval,de) ; if (fval.lt.0) sfact(3)= 1.d0
c
c  ...if close to first vertex, do not be too close!
      if ((Xi(1)+Xi(2)).lt.1.d-6) then
        Xi(1)=1.d-6 ; Xi(2)=1.d-6
      endif
c
c  ...1st EDGE : SURFS sph -> 3.........................................
      iv1=TRIAN_EDGE_TO_VERT(1,1) ; iv2=TRIAN_EDGE_TO_VERT(2,1)
      call proj_t2e(iv1,iv2,shapH,dshapH, t,dtdeta)
      norient=0 ; if (ncurv(1).lt.0) norient=1
      call curve_local(abs(ncurv(1)),norient,t, xc,dxcdt)
c
c
ccc      if (ncurv(1).gt.0) then
ccc        xic = Xi(1)+Xi(2) ; sign= 1.d0
ccc      else
ccc        xic = 1.d0 - Xi(1)-Xi(2)
ccc        sign= -1.d0
ccc      endif
ccc      call curve(abs(ncurv(1)),xic, xc,dxcdt)
ccc      dxcdt(1:3) = sign*dxcdt(1:3)
c
c
c  ...surfaces bounding 1st edge
      call sphere1(xc,ctr50,0.d0, fval2,der2)
      call surf(Nsutr(3),xc,      fval1,der1)
c  ...adjust orientation
      fval1=fval1*sfact(2) ; der1(1:3)=der1(1:3)*sfact(2)
c
c  ...value of stretching function and its derivative
      xif(1) = fval2/(fval2-fval1)
      xifdxi(1,1) = ((1.d0-xif(1))*dip(der2,dxcdt)
     .                   + xif(1) *dip(der1,dxcdt))/(fval2-fval1)
c
c  ...2nd EDGE : SURFS 2 -> 4...........................................
      iv1=TRIAN_EDGE_TO_VERT(1,2) ; iv2=TRIAN_EDGE_TO_VERT(2,2)
      call proj_t2e(iv1,iv2,shapH,dshapH, t,dtdeta)
      norient=0 ; if (ncurv(2).lt.0) norient=1
      call curve_local(abs(ncurv(2)),norient,t, xc,dxcdt)
c
c
ccc      if ((Xi(1)+Xi(2)).ne.0) then
ccc        xih = ((Xi(2)-Xi(1))/(Xi(2)+Xi(1))+1)/2
ccc      else
ccc        xih = 0
ccc      endif
ccc      if (ncurv(2).gt.0) then
ccc        xic = xih
ccc        sign=1.d0
ccc      else
ccc        xic = 1.d0 - xih
ccc        sign=-1.d0
ccc      endif
ccc      call curve(abs(ncurv(2)),xic, xc,dxcdt)
ccc      dxcdt(1:3) = sign*dxcdt(1:3)
c
c
c
c  ...surfaces bounding 2nd edge
      call surf(Nsutr(2),xc, fval1,der1)
      call surf(Nsutr(4),xc, fval2,der2)
c  ...adjust orientations
      fval1=fval1*sfact(1) ; der1(1:3)=der1(1:3)*sfact(1)
      fval2=fval2*sfact(3) ; der2(1:3)=der2(1:3)*sfact(3)
c
c  ...stretching function and its derivative
      xif(3) = fval1/(fval1-fval2)
      xifdxi(1,2) = ((1.d0-xif(3))*dip(der1,dxcdt)
     .                    + xif(3)*dip(der2,dxcdt))/(fval1-fval2)
c
c  ...3rd EDGE : SURFS sph -> 3.........................................
      iv1=TRIAN_EDGE_TO_VERT(1,3) ; iv2=TRIAN_EDGE_TO_VERT(2,3)
      call proj_t2e(iv1,iv2,shapH,dshapH, t,dtdeta)
      norient=0 ; if (ncurv(3).lt.0) norient=1
      call curve_local(abs(ncurv(3)),norient,t, xc,dxcdt)
c
c
ccc      if (ncurv(3).lt.0) then
ccc        xic = Xi(2)+Xi(1)
ccc        sign=1.d0
ccc      else
ccc        xic = 1.d0 - Xi(2)-Xi(1)
ccc        sign=-1.d0
ccc      endif
ccc      call curve(abs(ncurv(3)),xic, xc,dxcdt)
ccc      dxcdt(1:3) = sign*dxcdt(1:3)
c
c
c  ...surfaces bounding 3rd edge
      call sphere1(xc,ctr50,0.d0, fval2,der2)
      call surf(nsutr(3),xc, fval1,der1)
c  ...adjust orientation
      fval1=fval1*sfact(3) ; der1(1:3)=der1(1:3)*sfact(3)
c
c  ...stretching function and its derivative
      xif(4) = fval2/(fval2-fval1)
      xifdxi(2,1) = ((1.d0-xif(4))*dip(der2,dxcdt)
     .                    + xif(4)*dip(der1,dxcdt))/(fval2-fval1)

ccc      write(*,1111)nsutr(1:4),Xi,xif,xs
ccc 1111 format(' nsutr,Xi,xif,xs = ',4(i1,2x),2x,9(e12.5,2x))
      call mnewt(3,nsutr,Xi,xif,xs,sfact, X)
c
c
      if (iprint.eq.1) then
         write(*,*)'X =',X
         write(*,*)' xif', xif
         write(*,*)' xifdxi', xifdxi

c         pause
      endif
c
      if ((Xi(1)+Xi(2)).eq.0) return
c
c    ...compute derivatives

      call surf(nsutr(1),X, fval1,der1)
      aux1(1) = 0.d0
      aux2(1) = 0.d0
      aux(1,1) = der1(1)
      aux(1,2) = der1(2)
      aux(1,3) = der1(3)
c
      call surf(nsutr(2),X, fval1,der1)
      call surf(nsutr(4),X, fval2,der2)
      fval1=fval1*sfact(1) ; der1(1:3)=der1(1:3)*sfact(1)
      fval2=fval2*sfact(3) ; der2(1:3)=der2(1:3)*sfact(3)
      do  i=1,3
        aux(2,i) = (1.d0-xif(3))*der1(i) + xif(3)*der2(i)
      enddo
      aux1(2) = (fval2-fval1)*xifdxi(1,2)*Xi(2)/(Xi(2)+Xi(1))**2
      aux2(2) = (fval1-fval2)*xifdxi(1,2)*Xi(1)/(Xi(2)+Xi(1))**2
c
      call surf(nsutr(3),X, fval1,der1)
      fval1=fval1*sfact(3) ; der1(1:3)=der1(1:3)*sfact(3)
      call sphere1(X,ctr50,0.d0, fval2,der2)
      do i=1,3
       aux(3,i) = (Xi(1)*(1-xif(1))+Xi(2)*(1-xif(4)))*der2(i)
     .          + (Xi(1)*xif(1) + Xi(2)*xif(4))*der1(i)
      enddo
      aux2(3) = (fval2-fval1)*(Xi(2)*xifdxi(2,1)+Xi(1)*xifdxi(1,1)
     .        + xif(4)) - fval2
      aux1(3) = (fval2-fval1)*(Xi(2)*xifdxi(2,1)+Xi(1)*xifdxi(1,1)
     .        + xif(1)) - fval2

      call saruss(aux,aux1, Dxdxi(1,1),ifl1)
      call saruss(aux,aux2, Dxdxi(1,2),ifl2)
      if (ifl1.ne.0.or.ifl2.ne.0) then
        iprint=1
        go to 5
      endif
c
      end subroutine trian_ImpliTri
c
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - compute inner product
c
c----------------------------------------------------------------------
      function dip(A,B)
      implicit none
c
      real(8) :: dip
      real(8) :: A(3),B(3)
c
      dip = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
c
      end function dip
