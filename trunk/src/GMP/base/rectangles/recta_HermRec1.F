c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : this is a LEGACY routine. Do not use it, unless
c                    you revise it!
c//////////////////////////////////////////////////////////////////////
c
c
c------------------------------------------------------------------------
c
c   routine name       - recta_HermRec
c
c------------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates using a biquartic parametrization
c
c   arguments :
c     in:
c               No     - the rectangle number
c               Eta    - reference coordinates of a point
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the reference coordinates
c
c------------------------------------------------------------------------
c
      subroutine recta_HermRec1(No,Eta, X,Dxdeta)
c
      use GMP
      use control
      use parameters , only : MAXP
c
#include "syscom.blk"
      common /crecta_HermRec1/ iprint
c
      dimension Eta(2),X(3),Dxdeta(3,2)
c
c  ...vertex coordinates, first and second mixed derivatives at vertices
      dimension xv(3,2,2),av(3,2,2),bv(3,2,2),dv(3,2,2)
c
c  ...dof for curves in local coordinates
      dimension xedof(1:3,6,4)
c
c  ...dof for curve normals in local coordinates
      dimension rndof(1:3,4,4)
c
c  ...dof for the curve normal tangents in local coordinates
      dimension gdof(1:3,4)
c
c  ...normal vector and its derivative, same for the unit normal vector
      dimension rn(1:3),drn(1:3),rnmod(1:3),drnmod(1:3),en(1:3),den(1:3)
c
c  ...dof for vector g along an edge
      dimension gendf_dof(4),gdf_dof(4)
c
c  ...values of edge parametrizations and their derivatives
      dimension f(1:3),df(1:3),ddf(1:3),ef(1:3),def(1:3),
     .          g(1:3),dg(1:3),gmod(1:3),dgmod(1:3)
c
c  ...Hermite shape functions and their derivatives at a point
      dimension vHshap1(MAXP+1), dvHshap1(MAXP+1), void(MAXP+1)
      dimension vHshap2(MAXP+1), dvHshap2(MAXP+1)
c
c  ...fifth order C^2-conforming shape functions
      dimension vshap(MAXP+1),dvshap(MAXP+1),ddvshap(MAXP+1)
c
c  ...work space
      dimension temp(1:3)
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
      data iprint/0/
c
c------------------------------------------------------------------------
c
      if ((RECTANGLES(No)%Type.ne.'HermRec').or.(NDIM.ne.3)) then
        write(*,7000)
 7000   format('recta_HermRec: WRONG CALL')
        return
      endif
      if (iprint.eq.1) then
        write(*,7002) No, Eta
 7002   format('recta_HermRec: No, Eta = ',i5,2x,2f8.3)
      endif
c
c  ...evaluate vertex coordinates
      call pointr(RECTANGLES(No)%VertNo(1), xv(1:3,1,1))
      call pointr(RECTANGLES(No)%VertNo(2), xv(1:3,2,1))
      call pointr(RECTANGLES(No)%VertNo(3), xv(1:3,2,2))
      call pointr(RECTANGLES(No)%VertNo(4), xv(1:3,1,2))
c
c  ...collect the data for the curves
      do i=1,4
        select case(i)
        case(1,2)
          norient = 1
        case(3,4)
          norient =-1
        end select
        nc = iabs(RECTANGLES(No)%EdgeNo(i))
c
c  .....consistent orientations
        if (norient*RECTANGLES(No)%EdgeNo(i).gt.0) then
          do j=1,2
            np = CURVES(nc)%EndPoNo(j)
            xedof(1:3,j,i) = POINTS(np)%Rdata(1:3)
            rndof(1:3,j,i) = POINTS(np)%Rdata(4:6)
          enddo
          xedof(1:3,3,i) = CURVES(nc)%Rdata(1:3)
          xedof(1:3,4,i) = CURVES(nc)%Rdata(4:6)
          xedof(1:3,5,i) = CURVES(nc)%Rdata(7:9)
          xedof(1:3,6,i) = CURVES(nc)%Rdata(10:12)
          rndof(1:3,3,i) = CURVES(nc)%Rdata(13:15)
          rndof(1:3,4,i) = CURVES(nc)%Rdata(16:18)
c
c  .....opposite orientations
        else
          do j=1,2
            j1=imod(j+1,2)
            np = CURVES(nc)%EndPoNo(j1)
            xedof(1:3,j,i) = POINTS(np)%Rdata(1:3)
            rndof(1:3,j,i) = POINTS(np)%Rdata(4:6)
          enddo
          xedof(1:3,3,i) = -CURVES(nc)%Rdata(4:6)
          xedof(1:3,4,i) = -CURVES(nc)%Rdata(1:3)
          xedof(1:3,5,i) =  CURVES(nc)%Rdata(10:12)
          xedof(1:3,6,i) =  CURVES(nc)%Rdata(7:9)
          rndof(1:3,3,i) = -CURVES(nc)%Rdata(16:18)
          rndof(1:3,4,i) = -CURVES(nc)%Rdata(13:15)
        endif
c
c  ...end of loop through edges
      enddo
c
c  ...copy the data for the derivatives at the vertices
      av(1:3,1,1) = xedof(1:3,3,1)
      bv(1:3,1,1) = xedof(1:3,3,4)
      dv(1:3,1,1) = RECTANGLES(No)%Rdata(1:3)
c
      av(1:3,2,1) = xedof(1:3,4,1)
      bv(1:3,2,1) = xedof(1:3,3,2)
      dv(1:3,2,1) = RECTANGLES(No)%Rdata(4:6)
c
      av(1:3,2,2) = xedof(1:3,4,3)
      bv(1:3,2,2) = xedof(1:3,4,2)
      dv(1:3,2,2) = RECTANGLES(No)%Rdata(7:9)
c
      av(1:3,1,2) = xedof(1:3,3,3)
      bv(1:3,1,2) = xedof(1:3,4,4)
      dv(1:3,1,2) = RECTANGLES(No)%Rdata(10:12)
c
c  ...CHECK ORIENTATIONS
      call mixed_product(av(1:3,1,1),bv(1:3,1,1),rndof(1:3,1,1), s)
      if (s.lt.0.d0) then
c
c  .....modify the dof for the normals
        do ie=1,4
          do j=1,4
            rndof(1:3,j,ie) = -rndof(1:3,j,ie)
          enddo
        enddo
      endif
      if (iprint.eq.1) then
        do i=1,4
          write(*,8018) i,rndof(1:3,1:4,i)
 8018     format('i, rndof = ',i2,4(3f8.3,2x))
        enddo
        call pause
      endif
      call mixed_product(av(1:3,2,1),bv(1:3,2,1),rndof(1:3,1,2), s)
      if (s.lt.0.d0) then
        write(*,8038) No,2,RECTANGLES(No)%VertNo(1),
     .                RECTANGLES(No)%VertNo(2)
 8038   format('recta_HermRec: BAD ORIENTATIONS No,vert,np = ',
     .          i5,i2,2i5)
        call pause
      endif
      call mixed_product(av(1:3,2,2),bv(1:3,2,2),rndof(1:3,2,3), s)
      if (s.lt.0.d0) then
        write(*,8038) No,3,RECTANGLES(No)%VertNo(1),
     .                RECTANGLES(No)%VertNo(3)
        call pause
      endif
      call mixed_product(av(1:3,1,2),bv(1:3,1,2),rndof(1:3,2,4), s)
      if (s.lt.0.d0) then
        write(*,8038) No,4,RECTANGLES(No)%VertNo(1),
     .                RECTANGLES(No)%VertNo(4)
        call pause
      endif
c
c-----------------------------------------------------------------------
c
      if (iprint.eq.1) then
        write(*,7003)
 7003   format('recta_HermRec: VERTEX DATA = ')
        do j=1,2
        do i=1,2
          write(*,7004) i,j,xv(1:3,i,j)
 7004     format('i,j = ',2i2,2x,'xv = ',3f8.3)
          write(*,7005) av(1:3,i,j),bv(1:3,i,j),dv(1:3,i,j)
 7005     format('av,bv,dv = ',3(3f8.3,2x))
        enddo
        enddo
      endif
c
c  ...initiate
      X(1:3) = 0.d0; Dxdeta(1:3,1:2) = 0.d0
c
c  ...compute Hermite polynomials at the point
      call Hshape1(3,Eta(1), vHshap1,dvHshap1,void)
      call Hshape1(3,Eta(2), vHshap2,dvHshap2,void)
c
c  ...vertex contributions (notice negative signs!)
      do i=1,2
      do j=1,2
        X(1:3) = X(1:3) - xv(1:3,i,j)*vHshap1(i)  *vHshap2(j)
     .                  - av(1:3,i,j)*vHshap1(2+i)*vHshap2(j)
     .                  - bv(1:3,i,j)*vHshap1(i)  *vHshap2(2+j)
     .                  - dv(1:3,i,j)*vHshap1(2+i)*vHshap2(2+j)
        Dxdeta(1:3,1) = Dxdeta(1:3,1)
     .                  - xv(1:3,i,j)*dvHshap1(i)  *vHshap2(j)
     .                  - av(1:3,i,j)*dvHshap1(2+i)*vHshap2(j)
     .                  - bv(1:3,i,j)*dvHshap1(i)  *vHshap2(2+j)
     .                  - dv(1:3,i,j)*dvHshap1(2+i)*vHshap2(2+j)
        Dxdeta(1:3,2) = Dxdeta(1:3,2)
     .                  - xv(1:3,i,j)*vHshap1(i)  *dvHshap2(j)
     .                  - av(1:3,i,j)*vHshap1(2+i)*dvHshap2(j)
     .                  - bv(1:3,i,j)*vHshap1(i)  *dvHshap2(2+j)
     .                  - dv(1:3,i,j)*vHshap1(2+i)*dvHshap2(2+j)
        if (iprint.eq.1) then
          write(*,*) 'i,Dxdeta(1:3,1) = ',i,Dxdeta(1:3,1)
        endif
      enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,*) 'recta_HermRec: AFTER VERTEX CONTRIBUTIONS'
        write(*,7001) X(1:3), Dxdeta(1:3,1), Dxdeta(1:3,2)
 7001   format('X, Dxdeta = ',3(3f8.3,3x))
      endif
c
c  ...edge contributions
      do ie=1,4
        f(1:3) = 0.d0;  rn(1:3) = 0.d0; g(1:3) = 0.d0
        df(1:3) = 0.d0; drn(1:3) = 0.d0; dg(1:3) = 0.d0
        ddf(1:3) = 0.d0; gendf = 0.d0; gdf = 0.d0
        select case(ie)
        case(1)
          etac = Eta(1)
          gdof(1:3,1) = bv(1:3,1,1); gdof(1:3,2) = bv(1:3,2,1)
          gdof(1:3,3) = dv(1:3,1,1); gdof(1:3,4) = dv(1:3,2,1)
c
c  .......dof for component gendf...
          call get_gcomp(gdof(1:3,1),gdof(1:3,3),
     .         rndof(1:3,1,1),rndof(1:3,3,1),
     .         xedof(1:3,3,1),xedof(1:3,5,1),
     .         gendf_dof(1),gendf_dof(3),gdf_dof(1),gdf_dof(3))
          call get_gcomp(gdof(1:3,2),gdof(1:3,4),
     .         rndof(1:3,2,1),rndof(1:3,4,1),
     .         xedof(1:3,4,1),xedof(1:3,6,1),
     .         gendf_dof(2),gendf_dof(4),gdf_dof(2),gdf_dof(4))
          it=1; in=2
        case(2)
          etac = Eta(2)
          gdof(1:3,1) = av(1:3,2,1); gdof(1:3,2) = av(1:3,2,2)
          gdof(1:3,3) = dv(1:3,2,1); gdof(1:3,4) = dv(1:3,2,2)
          it=2; in=1
        case(3)
          etac = Eta(1)
          gdof(1:3,1) = bv(1:3,1,2); gdof(1:3,2) = bv(1:3,2,2)
          gdof(1:3,3) = dv(1:3,1,2); gdof(1:3,4) = dv(1:3,2,2)
          it=1; in=2
        case(4)
          etac = Eta(2)
          gdof(1:3,1) = av(1:3,1,1); gdof(1:3,2) = av(1:3,1,2)
          gdof(1:3,3) = dv(1:3,1,1); gdof(1:3,4) = dv(1:3,1,2)
          it=2; in=1
        end select
c
c  .....evaluate the fifth order polynomials
        call Gshape1(5,etac, vshap,dvshap,ddvshap)
c
c  .....compute the edge contributions
        do k=1,6
          f(1:3)   = f(1:3)   + xedof(1:3,k,ie)*vshap(k)
          df(1:3)  = df(1:3)  + xedof(1:3,k,ie)*dvshap(k)
          ddf(1:3) = ddf(1:3) + xedof(1:3,k,ie)*ddvshap(k)
        enddo
c
c  .....compute the unit tangent vector
        call scalar_product(df,df, s)
        s = sqrt(s)
        call scalar_product(df,ddf, s1)
        ef(1:3) = df(1:3)/s
        def(1:3) = ddf(1:3)/s - s1/s**3*df(1:3)
c
        select case(ie)
        case(1,3)
          do k=1,4
            rn(1:3)  = rn(1:3) + rndof(1:3,k,ie)*vHshap1(k)
            drn(1:3) = drn(1:3) + rndof(1:3,k,ie)*dvHshap1(k)
            g(1:3)  = g(1:3) + gdof(1:3,k)*vHshap1(k)
            dg(1:3) = dg(1:3) + gdof(1:3,k)*dvHshap1(k)
            gendf = gendf + gendf_dof(k)*vHshap1(k)
            gdf = gdf + gdf_dof(k)*vHshap1(k)
          enddo
        case(2,4)
          do k=1,4
            rn(1:3)  = rn(1:3) + rndof(1:3,k,ie)*vHshap2(k)
            drn(1:3) = drn(1:3) + rndof(1:3,k,ie)*dvHshap2(k)
            g(1:3)  = g(1:3) + gdof(1:3,k)*vHshap2(k)
            dg(1:3) = dg(1:3) + gdof(1:3,k)*dvHshap2(k)
          enddo
        end select
c
        if (iprint.eq.1) then
          write(*,7048) gendf,gdf
 7048     format('recta_HermRec: gendf,gdf = ',2f8.3)
        endif
c
c  .....enforce orthogonality between the normal and the tangent
c       unit vectors
        call scalar_product(rn,ef, s)
        if (iprint.eq.1) then
          write(*,7936) ie,s
 7936     format('recta_HermRec: ie, rn \dot ef= ',i2,f8.3)
        endif
        call scalar_product(drn,ef, s1)
        call scalar_product(rn,def, s2)
ccc        write(*,*) 'old rn,drn = ',rn,drn
        rn(1:3) = rn(1:3) - s*ef(1:3)
        drn(1:3) = drn(1:3) - (s1+s2)*ef(1:3) - s*def(1:3)
ccc        write(*,*) 'new rn,drn = ',rn,drn
ccc        call pause
c
c  .....verify
        call scalar_product(rn,ef, s)
        if (abs(s).gt.GEOM_TOL) then
          write(*,8936) No,Eta,ie,s
 8936     format('recta_HermRec: No,Eta,ie,rn*df = ',i4,2f8.3,i2,e12.5)
          call pause
        endif
c
c  .....compute the unit normal vector and its derivative
        call scalar_product(rn,rn, s)
        s = sqrt(s)
        call scalar_product(rn,drn, s1)
        en(1:3) = rn(1:3)/s
        den(1:3) = drn(1:3)/s - s1/s**3*rn(1:3)
c
c  .....modify the tangent vector to enforce orthogonality with the normal
        call scalar_product(g,en, s)
        call scalar_product(dg,en, s1)
        call scalar_product(g,den, s2)
        gmod(1:3) = g(1:3) - s*en(1:3)
        dgmod(1:3) = dg(1:3) - (s1+s2)*en(1:3) - s*den(1:3)
c
c  .....blending function and its derivative
        select case(ie)
        case(1)
          b1 = vHshap2(1); db1 = dvHshap2(1)
          b2 = vHshap2(3); db2 = dvHshap2(3)
        case(2)
          b1 = vHshap1(2); db1 = dvHshap1(2)
          b2 = vHshap1(4); db2 = dvHshap1(4)
        case(3)
          b1 = vHshap2(2); db1 = dvHshap2(2)
          b2 = vHshap2(4); db2 = dvHshap2(4)
        case(4)
          b1 = vHshap1(1); db1 = dvHshap1(1)
          b2 = vHshap1(3); db2 = dvHshap1(3)
        end select
        if (iprint.eq.1) then
          write(*,*) 'ie = ',ie
          write(*,7009) rn,en
 7009     format('recta_HermRec: rn,en = ',2(3f8.3,2x))
          write(*,7054) drn,den
 7054     format('recta_HermRec: drn,den = ',2(3f8.3,2x))
          write(*,7010) f, g, gmod
 7010     format('recta_HermRec: f,g,gmod = ',3(3f8.3,2x))
          write(*,7012) df, dg, dgmod
 7012     format('recta_HermRec: df,dg,dgmod = ',3(3f8.3,2x))
          write(*,7011) b1,db1,b2,db2
 7011     format('recta_HermRec: b1,db1,b2,db2 = ',2f8.3,2x,2f8.3)
          call scalar_product(dgmod,en,s1)
          call scalar_product(gmod,den,s2)
          write(*,*) 'd(g*en) = ',s1+s2
          call pause
        endif
c
c
        X(1:3) = X(1:3) + f(1:3)*b1 + gmod(1:3)*b2
        Dxdeta(1:3,it) = Dxdeta(1:3,it) + df(1:3)*b1 + dgmod(1:3)*b2
        Dxdeta(1:3,in) = Dxdeta(1:3,in) + f(1:3)*db1 + gmod(1:3)*db2
c
c  ...end of loop through edges
      enddo
      if (iprint.eq.1) then
        write(*,*) 'recta_HermRec: AFTER EDGE CONTRIBUTIONS'
        write(*,7001) X(1:3), Dxdeta(1:3,1), Dxdeta(1:3,2)
        call cross_product(Dxdeta(1:3,1), Dxdeta(1:3,2), temp)
        call scalar_product(temp,temp, s)
        s = sqrt(s)
        write(*,7006) temp(1:3)/s
 7006   format('recta_HermRec: NORMAL = ',3f8.3)
        call pause
      endif
c
c
      end
