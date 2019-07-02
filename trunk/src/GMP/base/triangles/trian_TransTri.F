c----------------------------------------------------------------------
c> Purpose : transfinite interpolation for a triangle
c
c  @param[in ] No        - triangle number
c  @param[in ] Eta       - reference coordinates of a point
c  @param[out] X         - physical coordinates of the point
c  @param[out] Dxdeta    - derivatives wrt reference coordinates
c
c  @revision Nov 12      
c----------------------------------------------------------------------
c
      subroutine trian_TransTri(No,Eta, X,Dxdeta)
c
      use control
      use GMP
      use element_data
#include "syscom.blk"
      common /cTI_triangle/ iprint
c
      dimension Eta(2), X(3),Dxdeta(3,2)
c
c  ...linear shape functions
      dimension shapH(3),dshapH(2,3)
c
c  ...point on edge
      dimension dseds(2),xc(3),dxcdse(3)
c
c  ...blending function
      dimension dblend(2)
c
c-----------------------------------------------------------------------
c
      iprint=0
c
c  ...evaluate linear shape functions
      shapH(1) = 1.d0-Eta(1)-Eta(2); dshapH(1:2,1) = -1.d0
      shapH(2) = Eta(1); dshapH(1,2) = 1.d0; dshapH(2,2) = 0.d0   
      shapH(3) = Eta(2); dshapH(1,3) = 0.d0; dshapH(2,3) = 1.d0  
c
      X(1:3) = 0.d0; Dxdeta(1:3,1:2) = 0.d0 
c
c  ...compute vertex interpolant
      do iv=1,3
        np = TRIANGLES(No)%VertNo(iv)
        X(1:3) = X(1:3) + POINTS(np)%Rdata(1:3)*shapH(iv)
        do j=1,2
          Dxdeta(1:3,j) = Dxdeta(1:3,j) 
     .                  + POINTS(np)%Rdata(1:3)*dshapH(j,iv)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'trian_TransTri: AFTER VERTICES X,Dxdeta = '
        do ivar=1,3
          write(*,7035) X(ivar),Dxdeta(ivar,1:2)
 7035     format(e12.5,3x,3e12.5)
        enddo
      endif
c
c  ...add edge bubbles
      do ie=1,3
c
c  .....get the curve number
        nc = TRIANGLES(No)%EdgeNo(ie); norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
        if (CURVES(nc)%Type.eq.'Seglin') cycle
c
c  .....get the edge vertices specifying the local edge orientation
        iv1=TRIAN_EDGE_TO_VERT(1,ie) ; iv2=TRIAN_EDGE_TO_VERT(2,ie)
c
c  .....project s onto the edge
        call proj_t2e(iv1,iv2,shapH,dshapH, se,dseds)
        if ((se.lt.GEOM_TOL).or.(se.gt.1.d0-GEOM_TOL)) cycle
        if (iprint.eq.1) then
          write(*,7003) ie,nc,CURVES(nc)%Type
 7003     format(' trian_TransTri: ie,nc,Type = ',i2,i5,2x,a5)
        endif
c
c  .....compute the kernel function
        call curveK(nc,se,norient, xc,dxcdse)
c
c  .....blending function
        blend = shapH(iv1)*shapH(iv2)
        dblend(1:2) = dshapH(1:2,iv1)*shapH(iv2)
     .              + shapH(iv1)*dshapH(1:2,iv2)
c
c  .....add adge contribution
        X(1:3) = X(1:3) + xc(1:3)*blend
        do j=1,2
          Dxdeta(1:3,j) = Dxdeta(1:3,j) 
     .                  + dxcdse(1:3)*dseds(j)*blend
     .                  + xc(1:3)*dblend(j)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'trian_TransTri: AFTER EDGES X,Dxdeta = '
        do ivar=1,3
          write(*,7035) X(ivar),Dxdeta(ivar,1:2)
        enddo
        call pause
      endif
c
c
      endsubroutine trian_TransTri
