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
c   computer           - machine independent
c
c   latest revision    - Oct 07
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
      subroutine recta_HermRec(No,Eta, X,Dxdeta)
c
      use GMP
      use control
      use parameters , only : MAXP
c
#include "syscom.blk"
      common /crecta_HermRec/ iprint
c
      dimension Eta(2),X(3),Dxdeta(3,2)
c
c  ...vertex coordinates, first and second derivatives at vertices
      dimension xv(3,2,2),
     .          av(3,2,2),aav(3,2,2),
     .          bv(3,2,2),bbv(3,2,2),
     .          dv(3,2,2)
c
c  ...dof for curves in local coordinates
      dimension xedof(1:3,6,4)
c
c  ...normals at vertices
      dimension rnvert(1:3,2,2)
c
c  ...dof for functions psi
      dimension psidof(1:3,4,4)
c
c  ...coefficients alpha,beta,gamma for edges
      dimension alpha(2,2),beta(2,2),deldof(2,4),gamdof(4,4)
c
c  ...dual basis vectors at a vertex
      dimension adual(3),bdual(3)
c
c  ...values of edge parametrizations and their derivatives
      dimension f(3),df(3),ddf(3),g(3),dg(3),psi(3),dpsi(3)
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
c  ...sharp edge flags for the rectangle edges
      dimension nsharp_edge(4)
c
c  ...integer endpoint coordinates for edges
      dimension iedge_vert(2,2,4)
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
      data iprint/0/
c
      data iedge_vert/ 1,1,2,1, 2,1,2,2, 1,2,2,2, 1,1,1,2/
c
c------------------------------------------------------------------------
c
      if ((RECTANGLES(No)%Type.ne.'HermRec').or.(NDIM.ne.3)) then
        write(*,7000)
 7000   format('recta_HermRec: WRONG CALL')
        return
      endif
 10   continue
c
      if (iprint.eq.1) then
        write(*,7002) No, Eta
 7002   format('recta_HermRec: No, Eta = ',i5,2x,2f8.3)
      endif
c
c-----------------------------------------------------------------------
c
c  Step 1: Collect the vertex and curve data
c
c  ...evaluate vertex coordinates
      call pointr(RECTANGLES(No)%VertNo(1), xv(1:3,1,1))
      call pointr(RECTANGLES(No)%VertNo(2), xv(1:3,2,1))
      call pointr(RECTANGLES(No)%VertNo(3), xv(1:3,2,2))
      call pointr(RECTANGLES(No)%VertNo(4), xv(1:3,1,2))
c
c  ...collect the data for the edge curves
      do i=1,4
        select case(i)
        case(1,2)
          norient = 1
        case(3,4)
          norient =-1
        end select
        nc = iabs(RECTANGLES(No)%EdgeNo(i))
        call check_edge(nc, nsharp_edge(i))
        if (iprint.eq.1) then
          write(*,*) 'i,nsharp_edge(i) = ',i,nsharp_edge(i)
        endif
c
c  .....consistent orientations
        if (norient*RECTANGLES(No)%EdgeNo(i).gt.0) then
          do j=1,2
            np = CURVES(nc)%EndPoNo(j)
            xedof(1:3,j,i) = POINTS(np)%Rdata(1:3)
          enddo
          xedof(1:3,3,i) = CURVES(nc)%Rdata(1:3)
          xedof(1:3,4,i) = CURVES(nc)%Rdata(4:6)
          xedof(1:3,5,i) = CURVES(nc)%Rdata(7:9)
          xedof(1:3,6,i) = CURVES(nc)%Rdata(10:12)
c
c  .....opposite orientations
        else
          do j=1,2
            j1=imod(j+1,2)
            np = CURVES(nc)%EndPoNo(j1)
            xedof(1:3,j,i) = POINTS(np)%Rdata(1:3)
          enddo
          xedof(1:3,3,i) = -CURVES(nc)%Rdata(4:6)
          xedof(1:3,4,i) = -CURVES(nc)%Rdata(1:3)
          xedof(1:3,5,i) =  CURVES(nc)%Rdata(10:12)
          xedof(1:3,6,i) =  CURVES(nc)%Rdata(7:9)
        endif
c
c  ...end of loop through edges
      enddo
c
c  ...collect the data for the element normals
      rnvert(1:3,1,1) = RECTANGLES(No)%Rdata(13:15)
      rnvert(1:3,2,1) = RECTANGLES(No)%Rdata(16:18)
      rnvert(1:3,2,2) = RECTANGLES(No)%Rdata(19:21)
      rnvert(1:3,1,2) = RECTANGLES(No)%Rdata(22:24)
      if (iprint.eq.1) then
        write(*,8038) ((rnvert(1:3,i,j),i=1,2),j=1,2)
 8038   format('rnvert = ',4(3f8.3,2x))
      endif
c
c
c  ...copy the data for the derivatives at the vertices
      av(1:3,1,1)  = xedof(1:3,3,1)
      aav(1:3,1,1) = xedof(1:3,5,1)
      bv(1:3,1,1)  = xedof(1:3,3,4)
      bbv(1:3,1,1) = xedof(1:3,5,4)
      dv(1:3,1,1)  = RECTANGLES(No)%Rdata(1:3)
c
      av(1:3,2,1)  = xedof(1:3,4,1)
      aav(1:3,2,1) = xedof(1:3,6,1)
      bv(1:3,2,1)  = xedof(1:3,3,2)
      bbv(1:3,2,1) = xedof(1:3,5,2)
      dv(1:3,2,1)  = RECTANGLES(No)%Rdata(4:6)
c
      av(1:3,2,2)  = xedof(1:3,4,3)
      aav(1:3,2,2) = xedof(1:3,6,3)
      bv(1:3,2,2)  = xedof(1:3,4,2)
      bbv(1:3,2,2) = xedof(1:3,6,2)
      dv(1:3,2,2)  = RECTANGLES(No)%Rdata(7:9)
c
      av(1:3,1,2)  = xedof(1:3,3,3)
      aav(1:3,1,2) = xedof(1:3,5,3)
      bv(1:3,1,2)  = xedof(1:3,4,4)
      bbv(1:3,1,2) = xedof(1:3,6,4)
      dv(1:3,1,2)  = RECTANGLES(No)%Rdata(10:12)
c
c  ...check orientations
      do i=1,2
      do j=1,2
        call mixed_product(av(1:3,i,j),bv(1:3,i,j),rnvert(1:3,i,j), s)
        if (s.le.0.d0) then
          rnvert(1:3,i,j) = - rnvert(1:3,i,j)
ccc          write(*,3078) i,j,s
ccc 3078     format('recta:HermRec: i,j,[av,bv,rn] = ',2i2,2x,e12.5)
ccc          call pause
        endif
      enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,7003) RECTANGLES(No)%VertNo(1:4)

 7003   format('recta_HermRec: VERTEX DATA = ',4i6)
        do j=1,2
        do i=1,2
          write(*,7004) i,j,xv(1:3,i,j),rnvert(1:3,i,j)
 7004     format('i,j = ',2i2,2x,'xv = ',3f8.3,'rnvert = ',3f8.3)
          write(*,7005) av(1:3,i,j),bv(1:3,i,j),dv(1:3,i,j)
 7005     format('av,bv,dv = ',3(3f8.3,2x))
        enddo
        enddo
      endif
c
c-----------------------------------------------------------------------
c
c  Step 2: Generate the data for edge functions psi, coefficients
c          alpha,beta and tangential components of the mixed derivatives
c          at the vertices
c
c  ...compute coefficients alpha for horizontal edges and beta for
c     vertical edges
      do j=1,2
      do i=1,2
        call scalar_product(av(1:3,i,j),av(1:3,i,j), s1)
        call scalar_product(bv(1:3,i,j),bv(1:3,i,j), s2)
        call scalar_product(av(1:3,i,j),bv(1:3,i,j), s3)
        alpha(i,j) = s3/s1; beta(i,j) = s3/s2
      enddo
      enddo
c
c  ...compute the in-plane components of the mixed derivatives at
c     the vertices
      do j=1,2
      do i=1,2
        i1 = imod(i+1,2); j1 = imod(j+1,2)
        dalpha = (alpha(i1,j) - alpha(i,j))*(-1)**i1
        dbeta  = (beta(i,j1) - beta(i,j))*(-1)**j1
        call scalar_product(aav(1:3,i,j),av(1:3,i,j), s1)
        call scalar_product(av(1:3,i,j),av(1:3,i,j), s2)
        acomp = s1*alpha(i,j)+s2*dalpha
        call scalar_product(bbv(1:3,i,j),bv(1:3,i,j), s1)
        call scalar_product(bv(1:3,i,j),bv(1:3,i,j), s2)
        bcomp = s1*beta(i,j)+s2*dbeta
c
c  .....dual basis vectors
        call cross_product(rnvert(1:3,i,j),bv(1:3,i,j), adual)
        call scalar_product(adual, av(1:3,i,j), s)
        adual(1:3) = adual(1:3)/s
        call cross_product(rnvert(1:3,i,j),av(1:3,i,j), bdual)
        call scalar_product(bdual, bv(1:3,i,j), s)
        bdual(1:3) = bdual(1:3)/s
c
c  .....check whether tangential components of the mixed derivative
c       vanish
        call scalar_product(dv(1:3,i,j),rnvert(1:3,i,j), s)
        temp(1:3) = dv(1:3,i,j) - s*rnvert(1:3,i,j)
        call norm(temp, s)
        if (s.gt.GEOM_TOL) then
          write(*,7057) i,j,s
 7057     format('recta_HermRec: WRONG dv, i,j,s = ',2i2,2x,e12.5)
          call pause
          iprint=1
          go to 10
        endif
c
c  .....adjust the mixed derivative
        dv(1:3,i,j) = dv(1:3,i,j) + acomp*adual(1:3) + bcomp*bdual(1:3)
      enddo
      enddo
c
c  ...function psi and coefficient gamma for horizontal edges
      do ie=1,3,2
        i1 = iedge_vert(1,1,ie); i2 = iedge_vert(1,2,ie)
        j1 = iedge_vert(2,1,ie); j2 = iedge_vert(2,2,ie)
        dalpha = alpha(i2,j1) - alpha(i1,j1)
        do iv=1,2
          i = iedge_vert(1,iv,ie); j = iedge_vert(2,iv,ie)
          deldof(iv,ie) = alpha(i,j)
c
c  .......compute vector psi and coefficient gamma
          temp(1:3) = av(1:3,i,j); call normalize(temp)
          call cross_product(rnvert(1:3,i,j),temp, psidof(1:3,iv,ie))
          call scalar_product(bv(1:3,i,j),psidof(1:3,iv,ie),
     .                        gamdof(iv,ie))
c
c  .......compute derivative of coefficient gamma
          call scalar_product(dv(1:3,i,j),psidof(1:3,iv,ie), s1)
          call scalar_product(aav(1:3,i,j),psidof(1:3,iv,ie), s2)
          gamdof(2+iv,ie) = s1 - s2*alpha(i,j)
c
c  .......compute derivative of vector psi
          psidof(1:3,2+iv,ie) = (dv(1:3,i,j) - aav(1:3,i,j)*alpha(i,j)
     .         - av(1:3,i,j)*dalpha - psidof(1:3,iv,ie)*gamdof(2+iv,ie))
     .           /gamdof(iv,ie)
c
c  .......check consistency, the in-plane components of the derivative
c         should be zero
          call scalar_product(psidof(1:3,2+iv,ie),rnvert(1:3,i,j), s)
          temp(1:3) = psidof(1:3,2+iv,ie) - s*rnvert(1:3,i,j)
          call norm(temp, s)
          if (s.gt.GEOM_TOL) then
            write(*,6028) i,j,s
 6028       format('recta_HermRec: WRONG HORIZONTAL dpsi, s = ',
     .              2i2,2x,e12.5)
            call pause
            iprint=1
            go to 10
          endif
        enddo
      enddo
c
c  ...function psi and coefficient gamma for vertical edges
      do ie=2,4,2
        i1 = iedge_vert(1,1,ie); i2 = iedge_vert(1,2,ie)
        j1 = iedge_vert(2,1,ie); j2 = iedge_vert(2,2,ie)
        dbeta = beta(i1,j2) - beta(i1,j1)
        do iv=1,2
          i = iedge_vert(1,iv,ie); j = iedge_vert(2,iv,ie)
          deldof(iv,ie) = beta(i,j)
c
c  .......compute vector psi and coefficient gamma
          temp(1:3) = bv(1:3,i,j); call normalize(temp)
          call cross_product(rnvert(1:3,i,j),temp, psidof(1:3,iv,ie))
          call scalar_product(av(1:3,i,j),psidof(1:3,iv,ie),
     .                        gamdof(iv,ie))
c
c  .......compute derivative of coefficient gamma
          call scalar_product(dv(1:3,i,j),psidof(1:3,iv,ie), s1)
          call scalar_product(bbv(1:3,i,j),psidof(1:3,iv,ie), s2)
          gamdof(2+iv,ie) = s1 - s2*beta(i,j)
c
c  .......compute derivative of vector psi
          psidof(1:3,2+iv,ie) = (dv(1:3,i,j) - bbv(1:3,i,j)*beta(i,j)
     .         - bv(1:3,i,j)*dbeta - psidof(1:3,iv,ie)*gamdof(2+iv,ie))
     .           /gamdof(iv,ie)
c
c  .......check consistency, the in-plane components of the derivative
c         should be zero
          call scalar_product(psidof(1:3,2+iv,ie),rnvert(1:3,i,j), s)
          temp(1:3) = psidof(1:3,2+iv,ie) - s*rnvert(1:3,i,j)
          call norm(temp, s)
          if (s.gt.GEOM_TOL) then
            write(*,6029) i,j,s
 6029       format('recta_HermRec: WRONG VERTICAL dpsi, s = ',
     .              2i2,2x,e12.5)
            call pause
            iprint=1
            go to 10
          endif
        enddo
      enddo
      if (iprint.eq.1) then
        do ie=1,4
          write(*,4001) ie
 4001     format('PSI,DELTA,GAMMA DOF FOR EDGE = ',i1)
          write(*,4002) (psidof(1:3,i,ie),i=1,4)
 4002     format('psi = ',4(3f8.3,2x))
          write(*,4003) (deldof(i,ie),i=1,2)
 4003     format('delta = ',2(3f8.3,2x))
          write(*,4004) (gamdof(i,ie),i=1,4)
 4004     format('gamma = ',4(3f8.3,2x))
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  Step 3: Compute the transformation and its derivatives
c
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
        f(1:3) = 0.d0;  g(1:3) = 0.d0
        df(1:3) = 0.d0; dg(1:3) = 0.d0
        ddf(1:3) = 0.d0
        select case(ie)
        case(1,3)
          etac = Eta(1)
          it=1; in=2
        case(2,4)
          etac = Eta(2)
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
        if (iprint.eq.1) then
          write(*,7047) df,ddf
 7047     format('df,ddf = ',3f8.3,3x,3f8.3)
        endif
c
c  .....compute the normal derivative (vector g)
        del = deldof(1,ie)*(1.d0-etac) + deldof(2,ie)*etac
        ddel = deldof(2,ie) - deldof(1,ie)
        psi(1:3) = 0.d0; dpsi(1:3) = 0.d0;
        gamma = 0.d0; dgamma = 0.d0
c
        select case(ie)
        case(1,3)
          do k=1,4
            psi(1:3)  = psi(1:3) + psidof(1:3,k,ie)*vHshap1(k)
            dpsi(1:3)  = dpsi(1:3) + psidof(1:3,k,ie)*dvHshap1(k)
            gamma = gamma + gamdof(k,ie)*vHshap1(k)
            dgamma = dgamma + gamdof(k,ie)*dvHshap1(k)
          enddo
        case(2,4)
          do k=1,4
            psi(1:3)  = psi(1:3) + psidof(1:3,k,ie)*vHshap2(k)
            dpsi(1:3)  = dpsi(1:3) + psidof(1:3,k,ie)*dvHshap2(k)
            gamma = gamma + gamdof(k,ie)*vHshap2(k)
            dgamma = dgamma + gamdof(k,ie)*dvHshap2(k)
          enddo
        end select
c
        if (iprint.eq.1) then
          write(*,5027) ie,gamma
 5027     format('ie,gamma = ',i2,2x,f8.3)
        endif
        g(1:3) = df(1:3)*del + psi(1:3)*gamma
        dg(1:3) = ddf(1:3)*del + df(1:3)*ddel
     .          + dpsi(1:3)*gamma + psi(1:3)*dgamma
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
          write(*,7010) f, g
 7010     format('recta_HermRec: f,g = ',2(3f8.3,2x))
        endif
c
c
        X(1:3) = X(1:3) + f(1:3)*b1 + g(1:3)*b2
        Dxdeta(1:3,it) = Dxdeta(1:3,it) + df(1:3)*b1 + dg(1:3)*b2
        Dxdeta(1:3,in) = Dxdeta(1:3,in) + f(1:3)*db1 + g(1:3)*db2
        if (iprint.eq.1) then
          write(*,*) 'recta_HermRec: AFTER CONTRIBUTION FROM EDGE',ie
          write(*,7001) X(1:3), Dxdeta(1:3,1), Dxdeta(1:3,2)
          call pause
        endif
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

