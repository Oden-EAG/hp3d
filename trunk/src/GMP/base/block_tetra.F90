!-----------------------------------------------------------------------
!
!   routine name       - tetra
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine evaluates physical coordinates
!                        and its derivatives wrt to reference
!                        coordinates for a point in reference
!                        tetrahedron
!
!   arguments :
!     in:
!               No     - a GMP tetrahedral block number
!               Eta    - reference coordinates of a point
!                        in the tetrahedron
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the parameters
!
!-----------------------------------------------------------------------
!
   subroutine tetra(No,Eta, X,Dxdeta)
!
      use GMP
      use node_types, only: TETR
      implicit none
!
      integer, intent(in)  :: No
      real(8), intent(in)  :: Eta(3)
      real(8), intent(out) :: X(3),Dxdeta(3,3)
!
!  ...vertex coordinates, edge curves
      real(8) :: xv(3,4)
!
!  ...linear shape functions for a tetrahedron
      real(8) :: vshap(8),dvshap(3,8)
!
      integer :: i,ie,ifc,ivar,k,lab,np,nc,nt
!
      integer :: iprint
      iprint=0
!
!-----------------------------------------------------------------------
!  ...initialize
      X = 0; Dxdeta = 0
!
      if (No.eq.7) iprint = 0
!
      if (iprint.eq.1) then
        write(*,7001) No,Eta,TETRAS(No)%Type
 7001   format('tetra: No,Eta,Type = ',i8,2x,3f8.3,2x,a10)
      endif
!
!  ...get the vertex coordinates
      do i=1,4
        np = TETRAS(No)%VertNo(i)
        call pointr(np, xv(1:3,i))
      enddo
!
!
      select case(TETRAS(No)%Type)
!
!  ...linear (affine) tetrahedron
      case('Linear')
!
!  .....check compatibility
        do ie=1,6
          nc = abs(TETRAS(No)%EdgeNo(ie))
          if (CURVES(nc)%Type.ne.'Seglin') then
            write(*,7002) No,ie,nc,CURVES(nc)%Type
 7002       format('tetra: INCOMPATIBLE DEFINITIONS: No,ie,nc=', &
                    i4,2x,i2,i4,' type = ',a10)
            stop
          endif
        enddo
        do ifc=1,4
          if (TETRAS(No)%FigNo(ifc).eq.0) cycle
          call decode(TETRAS(No)%FigNo(ifc), nt,lab)
          if (TRIANGLES(nt)%Type.ne.'PlaneTri') then
            write(*,7003) No,ifc,nt,TRIANGLES(nt)%Type
 7003       format('tetra: INCOMPATIBLE DEFINITIONS: No,ifc,nt=', &
                    i4,2x,i2,i4,' type = ',a10)
            stop
          endif
        enddo
!  .....get linear shape functions
        call vshape3(TETR,Eta, vshap,dvshap)
!
        X(1:3)=0.d0; Dxdeta(1:3,1:3) = 0.d0
        do k=1,4
          X(1:3) = X(1:3) + xv(1:3,k)*vshap(k)
          do i=1,3
            Dxdeta(1:3,i) = Dxdeta(1:3,i) + xv(1:3,k)*dvshap(i,k)
          enddo
        enddo
!
!  ...transfinite interpolation tetrahedron
      case('TraTet')
        call tetra_TraTet(No,Eta, X,Dxdeta)
!
!  ...cylindrical coordinates tetrahedron
      case('CylTet')
        call tetra_CylTet(No,Eta, X,Dxdeta)
!
      case default
        write(*,7010) TETRAS(No)%Type
 7010   format('tetra: Type = ',a10)
        stop
!
      end select
      if (iprint.eq.1) then
        do ivar=1,3
          write(*,7004) ivar,X(ivar),Dxdeta(ivar,1:3)
 7004     format('tetra: ivar,X,Dxdeta = ',i2,2x,e12.5,2x,3e12.5)
        enddo
        call pause
      endif
!
!
   end subroutine tetra
!
!
!-----------------------------------------------------------------------
!
!   routine name       - tetra_TraTet
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine evaluates physical coordinates
!                        and its derivatives wrt to reference
!                        coordinates for the transfinite interpolation
!                        tetrahedron
!
!   arguments :
!     in:
!               No     - a GMP tetrahedral block number
!               Eta    - reference coordinates of a point
!                        in the triangle
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the parameters
!
!-----------------------------------------------------------------------
!
   subroutine tetra_TraTet(No,Eta, X,Dxdeta)
!
      use control
      use GMP
      use element_data
      use node_types, only: TETR
      implicit none
!
      integer, intent(in)  :: No
      real(8), intent(in)  :: Eta(3)
      real(8), intent(out) :: X(3),Dxdeta(3,3)
!     real(8) :: xold(3),dxolddeta(3,3)
!
!  ...vertex shape functions
      real(8) :: shapH(8),dshapH(3,8)
!
!  ...derivatives of edge coordinate
      real(8) :: dtedeta(3)
!
!  ...face coordinates
      real(8) :: tf(2),dtfdeta(2,3)
!
!  ...edge kernels
      real(8) :: xe(3),dxedt(3)
!
!  ...face kernels
      real(8) :: xf(3),dxfdtf(3,2)
!
!  ...blending function
      real(8) :: dblend(1:3)
!
      real(8) :: blend,te
      integer :: ie,ifc,iv,iv1,iv2,iv3,ivar,nc,nt,norient,np
!
      integer :: iprint
      iprint=0
!
!------------------------------------------------------------------------
!  ...initialize
      X = 0; Dxdeta = 0
!
 10   continue
      if (iprint.eq.1) then
        write(*,*) '************************************************'
        write(*,7001) No,Eta(1:3)
 7001   format('tetra_TraTet: No,Eta = ',i5,2x,3e12.5)
      endif
!
!  ...check Eta
      if ( (Eta(1)              .lt.-GEOM_TOL).or. &
           (Eta(2)              .lt.-GEOM_TOL).or. &
           (Eta(3)              .lt.-GEOM_TOL).or. &
           (Eta(1)+Eta(2)+Eta(3).gt.1.d0+GEOM_TOL)    ) then
        write(*,7001) No,Eta(1:3)
        call pause
      endif
!
!  ...vertex shape functions (affine coordinates)
      call vshape3(TETR,Eta, shapH,dshapH)
!
!!!      k=1
!!!      shapH(k) = 1.d0 - Eta(1) - Eta(2) - Eta(3)
!!!      dshapH(1:3,k) = -1.d0
!!!      do i=1,3
!!!        k=k+1
!!!        shapH(k) = Eta(i)
!!!        dshapH(1:3,k) = 0.d0; dshapH(i,k) = 1.d0
!!!      enddo
!
!------------------------------------------------------------------------
!  ...start with linear interpolant
!------------------------------------------------------------------------
      X(1:3) = 0.d0; Dxdeta(1:3,1:3) = 0.d0
      do iv=1,4
        np = TETRAS(No)%VertNo(iv)
        X(1:3) = X(1:3) + POINTS(np)%Rdata(1:3)*shapH(iv)
        do ivar=1,3
          Dxdeta(ivar,1:3) = Dxdeta(ivar,1:3) &
                           + POINTS(np)%Rdata(ivar)*dshapH(1:3,iv)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'tetra_TraTet: VERTEX INTERPOLANT = '
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
 7011     format(e12.5,3x,3e12.5)
        enddo
      endif
!
!------------------------------------------------------------------------
!  ...add edge bubbles
!------------------------------------------------------------------------
      do ie=1,6
!
        nc = TETRAS(No)%EdgeNo(ie); norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
        if (CURVES(nc)%Type.eq.'Seglin') cycle
!
!  .....get the edge vertices specifying the local edge orientation
        iv1=TETRA_EDGE_TO_VERT(1,ie) ; iv2=TETRA_EDGE_TO_VERT(2,ie)
!
!  .....project Eta onto the edge
        call proj_n2e(iv1,iv2,shapH,dshapH, te,dtedeta)
        if ((abs(te).lt.GEOM_TOL).or.(abs(1.d0-te).lt.GEOM_TOL)) cycle
!
        if (iprint.eq.1) then
          write(*,7012) ie,nc,CURVES(nc)%Type
 7012     format('tetra_TraTet: ie,nc,Type = ',i2,i8,2x,a5)
        endif
!
!  .....evaluate edge kernel function
        call curveK(nc,te,norient, xe,dxedt)
        if (iprint.eq.1) then
          write(*,*) 'xe = ',xe
          write(*,*) 'dxedt = ',dxedt
        endif
!
!  .....blending function
        blend = shapH(iv1)*shapH(iv2)
        dblend(1:3) = dshapH(1:3,iv1)*shapH(iv2) &
                    + shapH(iv1)*dshapH(1:3,iv2)
!
!  .....add edge contribution
        X(1:3) = X(1:3) + xe(1:3)*blend
        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar) &
                           + dxedt(1:3)*dtedeta(ivar)*blend &
                           + xe(1:3)*dblend(ivar)
        enddo
      if (iprint.eq.1) then
        write(*,*) 'tetra_TraTet: AFTER EDGE ie = ',ie
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
        call pause
      endif
      enddo
!
!------------------------------------------------------------------------
!  ...add face bubbles
!------------------------------------------------------------------------
      do ifc=1,4
        call decode(TETRAS(No)%FigNo(ifc), nt,norient)
        if ((TRIANGLES(nt)%Type.eq.'TransTri').or. &
            (TRIANGLES(nt)%Type.eq.'PlaneTri')) cycle
!
!  .....get the vertices for the face defining its local orientation
        iv1=TETRA_FACE_TO_VERT(1,ifc)
        iv2=TETRA_FACE_TO_VERT(2,ifc)
        iv3=TETRA_FACE_TO_VERT(3,ifc)
!
!  .....project Eta onto the face
        call proj_n2f(iv1,iv2,iv3,shapH,dshapH, tf,dtfdeta)
!
!  .....if point is on the edge, then the bubble contribution is zero
        if ((abs(tf(2)).lt.GEOM_TOL).or.(abs(tf(1)).lt.GEOM_TOL).or. &
            (abs(1.d0-tf(1)-tf(2)).lt.GEOM_TOL)) cycle
        if (iprint.eq.1) then
          write(*,7013) ifc,nt,TRIANGLES(nt)%Type
 7013     format('tetra_TraTet: ifc,nt,Type = ',i2,i8,2x,a5)
        endif
!
!  .....compute the face kernel
        call trianK(nt,tf,norient, xf,dxfdtf)
!
!  .....blending function
        blend = shapH(iv1)*shapH(iv2)*shapH(iv3)
        dblend(1:3) = dshapH(1:3,iv1)*shapH(iv2)*shapH(iv3) &
                    + shapH(iv1)*dshapH(1:3,iv2)*shapH(iv3) &
                    + shapH(iv1)*shapH(iv2)*dshapH(1:3,iv3)
!
!  .....add face contribution
        X(1:3) = X(1:3) + xf(1:3)*blend
        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar) &
                           + (dxfdtf(1:3,1)*dtfdeta(1,ivar) &
                             +dxfdtf(1:3,2)*dtfdeta(2,ivar))*blend &
                           + xf(1:3)*dblend(ivar)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'tetra_TraTet: AFTER FACES = '
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
        call pause
      endif
!
!!!      call tetra_TraTet_OLD(No,Eta, xold,dxolddeta)
!!!      if (iprint.eq.1) then
!!!        write(*,*) 'tetra_TraTet: OLD RESULTS = '
!!!        do ivar=1,3
!!!          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
!!!        enddo
!!!        call pause
!!!      endif
!!!      smax=0.d0
!!!      do ivar=1,3
!!!        smax=max(smax,abs(X(ivar)-xold(ivar)))
!!!        do ieta=1,3
!!!          smax = max(smax,abs(Dxdeta(ivar,ieta)-dxolddeta(ivar,ieta)))
!!!        enddo
!!!      enddo
!!!      if (smax.gt.GEOM_TOL) then
!!!        write(*,*) 'tetra_TraTet: DEBUGGING ?(1/0)'
!!!        read(*,*) ians
!!!        if (ians.eq.1) then
!!!          iprint=1
!!!          goto 10
!!!        endif
!!!      endif
!
   end subroutine tetra_TraTet
!
!-----------------------------------------------------------------------
!
!   routine name       - tetra_CylTet
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine evaluates physical coordinates
!                        and its derivatives wrt to reference
!                        coordinates for a point on the image
!                        of a linear tet through a global system
!                        of coordinates: x,y=rcos(\theta),z=rsin(\theta)
!                        and their  derivative wrt to reference
!                        coordinates

!
!   arguments :
!     in:
!               No     - a GMP tetrahedral block number
!               Eta    - reference coordinates of a point
!                        in the triangle
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the parameters
!
!-----------------------------------------------------------------------
!
   subroutine tetra_CylTet(No,Eta, X,Dxdeta)
!
      use control
      use GMP
!
      implicit none
!
      integer, intent(in)  :: No
      real(8), intent(in)  :: Eta(3)
      real(8), intent(out) :: X(3),Dxdeta(3,3)
!
!  ...cylindrical coordinates of the endpoints of the triangle
      real(8) :: xp(3,4)
!
!  ...miscellanea
      real(8) :: pi,twopi,r,dr_deta1,dr_deta2,dr_deta3
      real(8) :: theta,sinthet,costhet
      real(8) :: dthet_deta1,dthet_deta2,dthet_deta3
      integer :: np,iv
!
      integer :: iprint
      iprint=0
!
!-----------------------------------------------------------------------
!  ...initialize
      X = 0; Dxdeta = 0
!
      if ((TETRAS(No)%Type.ne.'CylTet'.or.(NDIM.ne.3))) then
        write(*,7001) TETRAS(No)%Type
 7001   format('tria_CylTet: WRONG TETRA TYPE = ',a10)
        stop 1
      endif
!
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format('trian_CylTet: No,Eta = ',i4,2x,3f8.3)
      endif
!
      pi = acos(-1.d0)
      twopi = pi*2.d0
!
!  ...compute the cylindrical coordinates of the endpoints
      do iv=1,4
        np = TETRAS(No)%VertNo(iv)
        xp(1,iv) = POINTS(np)%Rdata(1)
        call coord_cart2polar(POINTS(np)%Rdata(2:3), r,theta)
        xp(2,iv) = r
!
!  .....adjust the angle, if necessary
        select case(iv)
        case(2,3,4)
          if (theta-xp(3,1).gt. pi) theta = theta - twopi
          if (theta-xp(3,1).lt.-pi) theta = theta + twopi
        end select
        xp(3,iv) = theta
      enddo
!
!  ...interpolate
      X(1)  = (1.d0-Eta(1)-Eta(2)-Eta(3))*xp(1,1) &
            + Eta(1)*xp(1,2) + Eta(2)*xp(1,3) + Eta(3)*xp(1,4)
      r     = (1.d0-Eta(1)-Eta(2)-Eta(3))*xp(2,1) &
            + Eta(1)*xp(2,2) + Eta(2)*xp(2,3) + Eta(3)*xp(2,4)
      theta = (1.d0-Eta(1)-Eta(2)-Eta(3))*xp(3,1) &
            + Eta(1)*xp(3,2) + Eta(2)*xp(3,3) + Eta(3)*xp(3,4)
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
!
!
   end subroutine tetra_CylTet


