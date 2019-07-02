c-----------------------------------------------------------------------
c
c   routine name       - tetra
c
c-----------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for a point in reference
c                        tetrahedron
c
c   arguments :
c     in:
c               No     - a GMP tetrahedral block number
c               Eta    - reference coordinates of a point
c                        in the tetrahedron
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the parameters
c
c-----------------------------------------------------------------------
c
      subroutine tetra(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
c
      dimension Eta(3),X(3),Dxdeta(3,3)
c
c  ...vertex coordinates, edge curves
      dimension xv(3,4)
c
c  ...linear shape functions for a tetrahedron
      dimension vshap(8),dvshap(3,8)
c
      if (No.eq.7) then
        iprint=0
      else
        iprint=0
      endif
c
      if (iprint.eq.1) then
        write(*,7001) No,Eta,TETRAS(No)%Type
 7001   format('tetra: No,Eta,Type = ',i8,2x,3f8.3,2x,a10)
      endif
c
c  ...get the vertex coordinates
      do i=1,4
        np = TETRAS(No)%VertNo(i)
        call pointr(np, xv(1:3,i))
      enddo
c
c
      select case(TETRAS(No)%Type)
c
c  ...linear (affine) tetrahedron
      case('Linear')
c
c  .....check compatibility
        do ie=1,6
          nc = iabs(TETRAS(No)%EdgeNo(ie))
          if (CURVES(nc)%Type.ne.'Seglin') then
            write(*,7002) No,ie,nc,CURVES(nc)%Type
 7002       format('tetra: INCOMPATIBLE DEFINITIONS: No,ie,nc=',
     .              i4,2x,i2,i4,' type = ',a10)
            stop
          endif
        enddo
        do if=1,4
          if (TETRAS(No)%FigNo(if).eq.0) cycle
          call decode(TETRAS(No)%FigNo(if), nt,lab)
          if (TRIANGLES(nt)%Type.ne.'PlaneTri') then
            write(*,7003) No,if,nt,TRIANGLES(nt)%Type
 7003       format('tetra: INCOMPATIBLE DEFINITIONS: No,if,nt=',
     .              i4,2x,i2,i4,' type = ',a10)
            stop
          endif
        enddo
c  .....get linear shape functions
        call vshape3('tetr',Eta, vshap,dvshap)
c
        X(1:3)=0.d0; Dxdeta(1:3,1:3) = 0.d0
        do k=1,4
          X(1:3) = X(1:3) + xv(1:3,k)*vshap(k)
          do i=1,3
            Dxdeta(1:3,i) = Dxdeta(1:3,i) + xv(1:3,k)*dvshap(i,k)
          enddo
        enddo
c
c  ...transfinite interpolation tetrahedron
      case('TraTet')
        call tetra_TraTet(No,Eta, X,Dxdeta)
c
c  ...cylindrical coordinates tetrahedron
      case('CylTet')
        call tetra_CylTet(No,Eta, X,Dxdeta)
c
      case default
        write(*,7010) TETRAS(No)%Type
 7010   format('tetra: Type = ',a10)
        stop
c
      end select
      if (iprint.eq.1) then
        do ivar=1,3
          write(*,7004) ivar,X(ivar),Dxdeta(ivar,1:3)
 7004     format('tetra: ivar,X,Dxdeta = ',i2,2x,e12.5,2x,3e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine tetra
c
c
c-----------------------------------------------------------------------
c
c   routine name       - tetra_TraTet
c
c-----------------------------------------------------------------------
c
c   latest revision    - Dec 08
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for the transfinite interpolation
c                        tetrahedron
c
c   arguments :
c     in:
c               No     - a GMP tetrahedral block number
c               Eta    - reference coordinates of a point
c                        in the triangle
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the parameters
c
c-----------------------------------------------------------------------
c
      subroutine tetra_TraTet(No,Eta, X,Dxdeta)
c
      use control
      use GMP
      use element_data
#include "syscom.blk"
      common /ctrianK/ iprint_trianK
      common /ctetra_TraTet/ iprint
c
      dimension Eta(3),X(3),Dxdeta(3,3)
      dimension xold(3),dxolddeta(3,3)
c
c  ...vertex shape functions
      dimension shapH(8),dshapH(3,8)
c
c  ...derivatives of edge coordinate
      dimension dtedeta(3)
c
c  ...face coordinates
      dimension tf(2),dtfdeta(2,3)
c
c  ...edge kernels
      dimension xe(3),dxedt(3)
c
c  ...face kernels
      dimension xf(3),dxfdtf(3,2)
c
c  ...blending function
      dimension dblend(1:3)
c
c------------------------------------------------------------------------
c
ccc      iprint=0
 10   continue
      iprint_trianK=iprint
ccc      iprint_TraTet_OLD = iprint
      if (iprint.eq.1) then
        write(*,*) '************************************************'
        write(*,7001) No,Eta(1:3)
 7001   format('tetra_TraTet: No,Eta = ',i5,2x,3e12.5)
      endif
c
c  ...check Eta
      if ( (Eta(1)              .lt.-GEOM_TOL).or.
     .     (Eta(2)              .lt.-GEOM_TOL).or.
     .     (Eta(3)              .lt.-GEOM_TOL).or.
     .     (Eta(1)+Eta(2)+Eta(3).gt.1.d0+GEOM_TOL)    ) then
        write(*,7001) No,Eta(1:3)
        call pause
      endif
c
c  ...vertex shape functions (affine coordinates)
      call vshape3('tetr',Eta, shapH,dshapH)
c
ccc      k=1
ccc      shapH(k) = 1.d0 - Eta(1) - Eta(2) - Eta(3)
ccc      dshapH(1:3,k) = -1.d0
ccc      do i=1,3
ccc        k=k+1
ccc        shapH(k) = Eta(i)
ccc        dshapH(1:3,k) = 0.d0; dshapH(i,k) = 1.d0
ccc      enddo
c
c------------------------------------------------------------------------
c  ...start with linear interpolant
c------------------------------------------------------------------------
      X(1:3) = 0.d0; Dxdeta(1:3,1:3) = 0.d0
      do iv=1,4
        np = TETRAS(No)%VertNo(iv)
        X(1:3) = X(1:3) + POINTS(np)%Rdata(1:3)*shapH(iv)
        do ivar=1,3
          Dxdeta(ivar,1:3) = Dxdeta(ivar,1:3)
     .                     + POINTS(np)%Rdata(ivar)*dshapH(1:3,iv)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'tetra_TraTet: VERTEX INTERPOLANT = '
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
 7011     format(e12.5,3x,3e12.5)
        enddo
      endif
c
c------------------------------------------------------------------------
c  ...add edge bubbles
c------------------------------------------------------------------------
      do ie=1,6
c
        nc = TETRAS(No)%EdgeNo(ie); norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
        if (CURVES(nc)%Type.eq.'Seglin') cycle
c
c  .....get the edge vertices specifying the local edge orientation
        iv1=TETRA_EDGE_TO_VERT(1,ie) ; iv2=TETRA_EDGE_TO_VERT(2,ie)
c
c  .....project Eta onto the edge
        call proj_n2e(iv1,iv2,shapH,dshapH, te,dtedeta)
        if ((abs(te).lt.GEOM_TOL).or.(abs(1.d0-te).lt.GEOM_TOL)) cycle
c
        if (iprint.eq.1) then
          write(*,7012) ie,nc,CURVES(nc)%Type
 7012     format('tetra_TraTet: ie,nc,Type = ',i2,i8,2x,a5)
        endif
c
c  .....evaluate edge kernel function
        call curveK(nc,te,norient, xe,dxedt)
        if (iprint.eq.1) then
          write(*,*) 'xe = ',xe
          write(*,*) 'dxedt = ',dxedt
        endif
c
c  .....blending function
        blend = shapH(iv1)*shapH(iv2)
        dblend(1:3) = dshapH(1:3,iv1)*shapH(iv2)
     .              + shapH(iv1)*dshapH(1:3,iv2)
c
c  .....add edge contribution
        X(1:3) = X(1:3) + xe(1:3)*blend
        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar)
     .                     + dxedt(1:3)*dtedeta(ivar)*blend
     .                     + xe(1:3)*dblend(ivar)
        enddo
      if (iprint.eq.1) then
        write(*,*) 'tetra_TraTet: AFTER EDGE ie = ',ie
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
        call pause
      endif
      enddo
c
c------------------------------------------------------------------------
c  ...add face bubbles
c------------------------------------------------------------------------
      do if=1,4
        call decode(TETRAS(No)%FigNo(if), nt,norient)
        if ((TRIANGLES(nt)%Type.eq.'TransTri').or.
     .      (TRIANGLES(nt)%Type.eq.'PlaneTri')) cycle
c
c  .....get the vertices for the face defining its local orientation
        iv1=TETRA_FACE_TO_VERT(1,if)
        iv2=TETRA_FACE_TO_VERT(2,if)
        iv3=TETRA_FACE_TO_VERT(3,if)
c
c  .....project Eta onto the face
        call proj_n2f(iv1,iv2,iv3,shapH,dshapH, tf,dtfdeta)
c
c  .....if point is on the edge, then the bubble contribution is zero
        if ((abs(tf(2)).lt.GEOM_TOL).or.(abs(tf(1)).lt.GEOM_TOL).or.
     .      (abs(1.d0-tf(1)-tf(2)).lt.GEOM_TOL)) cycle
        if (iprint.eq.1) then
          write(*,7013) if,nt,TRIANGLES(nt)%Type
 7013     format('tetra_TraTet: if,nt,Type = ',i2,i8,2x,a5)
        endif
c
c  .....compute the face kernel
        call trianK(nt,tf,norient, xf,dxfdtf)
c
c  .....blending function
        blend = shapH(iv1)*shapH(iv2)*shapH(iv3)
        dblend(1:3) = dshapH(1:3,iv1)*shapH(iv2)*shapH(iv3)
     .              + shapH(iv1)*dshapH(1:3,iv2)*shapH(iv3)
     .              + shapH(iv1)*shapH(iv2)*dshapH(1:3,iv3)
c
c  .....add face contribution
        X(1:3) = X(1:3) + xf(1:3)*blend
        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar)
     .                     + (dxfdtf(1:3,1)*dtfdeta(1,ivar)
     .                       +dxfdtf(1:3,2)*dtfdeta(2,ivar))*blend
     .                     + xf(1:3)*dblend(ivar)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'tetra_TraTet: AFTER FACES = '
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
        call pause
      endif
c
ccc      call tetra_TraTet_OLD(No,Eta, xold,dxolddeta)
ccc      if (iprint.eq.1) then
ccc        write(*,*) 'tetra_TraTet: OLD RESULTS = '
ccc        do ivar=1,3
ccc          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
ccc        enddo
ccc        call pause
ccc      endif
ccc      smax=0.d0
ccc      do ivar=1,3
ccc        smax=max(smax,abs(X(ivar)-xold(ivar)))
ccc        do ieta=1,3
ccc          smax = max(smax,abs(Dxdeta(ivar,ieta)-dxolddeta(ivar,ieta)))
ccc        enddo
ccc      enddo
ccc      if (smax.gt.GEOM_TOL) then
ccc        write(*,*) 'tetra_TraTet: DEBUGGING ?(1/0)'
ccc        read(*,*) ians
ccc        if (ians.eq.1) then
ccc          iprint=1
ccc          go to 10
ccc        endif
ccc      endif
c
      end subroutine tetra_TraTet
c
c-----------------------------------------------------------------------
c
c   routine name       - tetra_CylTet
c
c-----------------------------------------------------------------------
c
c   latest revision    - Dec 15
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        coordinates for a point on the image
c                        of a linear tet through a global system
c                        of coordinates: x,y=rcos(\theta),z=rsin(\theta)
c                        and their  derivative wrt to reference
c                        coordinates

c
c   arguments :
c     in:
c               No     - a GMP tetrahedral block number
c               Eta    - reference coordinates of a point
c                        in the triangle
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the parameters
c
c-----------------------------------------------------------------------
c
      subroutine tetra_CylTet(No,Eta, X,Dxdeta)
c
      use control
      use GMP
#include "syscom.blk"
      common /ctetra_CylTet/ iprint
c
      dimension Eta(3),X(3),Dxdeta(3,3)
c
c  ...cylindrical coordinates of the endpoints of the triangle
      dimension xp(3,4)
c
      iprint=0
c
      if ((TETRAS(No)%Type.ne.'CylTet'.or.(NDIM.ne.3))) then
        write(*,7001) TETRAS(No)%Type
 7001   format('tria_CylTet: WRONG TETRA TYPE = ',a10)
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('trian_CylTet: No,Eta = ',i4,2x,3f8.3)
      endif
c
      pi = acos(-1.d0)
      twopi = pi*2.d0
c
c  ...compute the cylindrical coordinates of the endpoints
      do iv=1,4
        np = TETRAS(No)%VertNo(iv)
        xp(1,iv) = POINTS(np)%Rdata(1)
        call coord_cart2polar(POINTS(np)%Rdata(2:3), r,theta)
        xp(2,iv) = r
c
c  .....adjust the angle, if necessary
        select case(iv)
        case(2,3,4)
          if (theta-xp(3,1).gt.pi) theta = theta - twopi
          if (theta-xp(3,1).lt.-pi) theta = theta + twopi
        end select
        xp(3,iv) = theta
      enddo
c
c  ...interpolate
      X(1)  = (1.d0-Eta(1)-Eta(2)-Eta(3))*xp(1,1)
     .      + Eta(1)*xp(1,2) + Eta(2)*xp(1,3) + Eta(3)*xp(1,4)
      r     = (1.d0-Eta(1)-Eta(2)-Eta(3))*xp(2,1)
     .      + Eta(1)*xp(2,2) + Eta(2)*xp(2,3) + Eta(3)*xp(2,4)
      theta = (1.d0-Eta(1)-Eta(2)-Eta(3))*xp(3,1)
     .      + Eta(1)*xp(3,2) + Eta(2)*xp(3,3) + Eta(3)*xp(3,4)
      dr_deta1    = xp(2,2) - xp(2,1)
      dr_deta2    = xp(2,3) - xp(2,1)
      dr_deta3    = xp(2,4) - xp(2,1)
      dthet_deta1 = xp(3,2) - xp(3,1)
      dthet_deta2 = xp(3,3) - xp(3,1)
      dthet_deta3 = xp(3,4) - xp(3,1)
      costhet = cos(theta); sinthet = sin(theta)
      X(2) = r*costhet
      X(3) = r*sinthet
      Dxdeta(1,1) = xp(1,2) - xp(1,1)
      Dxdeta(1,2) = xp(1,3) - xp(1,1)
      Dxdeta(1,3) = xp(1,4) - xp(1,1)
      Dxdeta(2,1) = dr_deta1*costhet - r*sinthet*dthet_deta1
      Dxdeta(2,2) = dr_deta2*costhet - r*sinthet*dthet_deta2
      Dxdeta(2,3) = dr_deta3*costhet - r*sinthet*dthet_deta3
      Dxdeta(3,1) = dr_deta1*sinthet + r*costhet*dthet_deta1
      Dxdeta(3,2) = dr_deta2*sinthet + r*costhet*dthet_deta2
      Dxdeta(3,3) = dr_deta3*sinthet + r*costhet*dthet_deta3
c
c
      end subroutine tetra_CylTet


