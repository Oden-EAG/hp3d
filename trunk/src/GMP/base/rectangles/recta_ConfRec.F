c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : this is a LEGACY routine. Do not use it, unless
c                    you revise it!
c//////////////////////////////////////////////////////////////////////
c
c
c----------------------------------------------------------------------
c
c   routine name       - recta_ConfRec
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 06
c
c   purpose            - routine evaluates physical coordinates
c                        and its derivatives wrt to reference
c                        rectangle parameters for a point in
c                        a rectangle conforming to a surface
c
c   arguments :
c     in:
c               No     - the rectangle number
c               Xi     - coordinates of a point in the reference
c                        rectrangle
c     out:
c               X      - physical coordinates of the point
c               Dxdxi  - derivatives of the physical coordinates wrt
c                        to the reference rectrangle coordinates
c
c----------------------------------------------------------------------
c
      subroutine recta_ConfRec(No,Xi, X,Dxdxi)
c
      use GMP
      use control
#include "syscom.blk"
      common /cCylGeod/ iprint_CylGeod
      common /ccurve/ iprint_curve
c
      dimension Xi(2), X(3),Dxdxi(3,2)
c
c  ...transformation matrix from global to local system of coordinates
      dimension transf(3,3)
c
c  ...cylinder origin and axis, auxiliary space
      dimension xcenter(3),xaxis(3),xvert(3,4),aux(3)
c
c  ...global coordinates of a point on a curve and their derivatives
      dimension xc(3),dxcdxi(3)
c
c  ...coordinates of a point and their derivatives in the local
c     Cartesian system
      dimension xp(3),dxpdxi(3,2)
c
c  ...u coordinates of the rectangle vertices
      dimension uv(2,4)
c
c  ...parametrizations of the edges of the rectangle in terms
c     of the rectangle explicit parametrization parameters
      dimension ue(2,4),duedxi(2,4)
c
c  ...coordinates of a point and its derivatives in the surface
c     curvilinear coordinates
      dimension u(2),dudxi(2,2)
c
c
      iprint=0
 10   continue
      iprint_CylGeod = iprint
      iprint_curve = iprint
      if (iprint.eq.1) then
        write(*,7000) No, Xi
 7000   format('recta_conf: DEBUGGING for No, Xi = ',i4,2f8.3)
      endif
c
      if (RECTANGLES(No)%Type.ne.'ConfRec') then
        write(*,7001) RECTANGLES(No)%Type
 7001   format('recta_conf: WRONG RECTANGLE TYPE = ',a10)
        stop 1
      endif
c
c  ...get the surface to conform to
      ns = RECTANGLES(No)%Idata(1)
c
c
      select case(SURFACES(ns)%Type)
c
c  .....a cylinder
        case('Cylinder')
c
c  .....collect cylinder data
        xcenter(1:3) = SURFACES(ns)%Rdata(1:3)
        xaxis(1:3)   = SURFACES(ns)%Rdata(4:6)
        radius   = SURFACES(ns)%Rdata(7)
c
c  .....check if the rectangle points are indeed on the cylinder
        do iv=1,4
          nv = RECTANGLES(No)%VertNo(iv)
          call pointr(nv, xvert(1:3,iv))
          call surf(ns,xvert(1:3,iv), fval,aux)
          if (abs(fval).gt.GEOM_TOL) then
            write(*,7002) nv,Ns,fval
 7002       format('recta_conf: POINT ',i4' NOT ON SURFACE ',i4,
     .             'fval = ',e12.5)
            stop 1
          endif
        enddo
c
c  .....compute the transformation matrix
c
c  .....the third local cylinder versor coincides with its axis
        s = sqrt(xaxis(1)**2+xaxis(2)**2+xaxis(3)**2)
        transf(1:3,3) = xaxis(1:3)/s
c
c  .....compute cross product of the cylinder vector with the vector
c       connecting its origin with the first curve endpoint
        aux(1:3) = xvert(1:3,1) - xcenter(1:3)
        transf(1,2) = xaxis(2)*aux(3) - xaxis(3)*aux(2)
        transf(2,2) = xaxis(3)*aux(1) - xaxis(1)*aux(3)
        transf(3,2) = xaxis(1)*aux(2) - xaxis(2)*aux(1)
c
c  .....normalize
        s = sqrt(transf(1,2)**2+transf(2,2)**2+transf(3,2)**2)
        transf(1:3,2) = transf(1:3,2)/s
c
c  .....compute the first axis as the cross product of the second and third
c       axes versors
        transf(1,1) = transf(2,2)*transf(3,3) - transf(3,2)*transf(2,3)
        transf(2,1) = transf(3,2)*transf(1,3) - transf(1,2)*transf(3,3)
        transf(3,1) = transf(1,2)*transf(2,3) - transf(2,2)*transf(1,3)
        if (iprint.eq.1) then
          do i=1,3
            write(*,7003) i,transf(1:3,i)
 7003       format('recta_conf: ',i1,'-TH VERSOR = ',3f8.3)
          enddo
        endif
c
c  Step 1: compute r,theta coordinates of the rectangle vertices
        do iv=1,4
c
c  .......transform to the local Cartesian system
          aux(1:3) = xvert(1:3,iv) - xcenter(1:3)
          call transform(12,transf,aux,xp)
c
c  .......compute r,theta coordinates
          call cart_to_polar(xp(1:2), void,uv(1,iv))
          if (iprint.eq.1) then
            write(*,7004) iv, xvert(1:3,iv), xp(1:3), void, uv(1,iv)
 7004       format('recta_conf: iv,xvert,xp,r,theta = ',
     .              i2,2x,3e12.5,2x,3e12.5,2x,2e12.5)
          endif
          if (abs(void-radius).gt.GEOM_TOL) then
            write(*,*) 'recta_conf: vertices: void - radius = ',
     .                  void-radius
            call pause
            iprint=1
            go to 10
          endif
          uv(2,iv) = xp(3)
        enddo
c
c  Step 2: compute theta,z coordinates of the edge curves and their
c          derivatives wrt reference coordinates
c
c  .....loop through edges
        do ie=1,4
c
c  .......compute edge reference coordinate and starting value for
c         NR iterations
          select case(ie)
          case(1)
            xil = Xi(1)
          case(2)
            xil = Xi(2)
          case(3)
            xil = 1.d0-Xi(1)
          case(4)
            xil = 1.d0-Xi(2)
          end select
c
c  .......find the global Cartesian coordinates of the point on the curve
          nc = RECTANGLES(No)%EdgeNo(ie)
          if (nc.gt.0) then
            call curve(nc,xil, xc,dxcdxi)
          else
            call curve(-nc,1.d0-xil, xc,dxcdxi)
            select case(ie)
            case(1,2)
              dxcdxi(1:3) = - dxcdxi(1:3)
            case(3,4)
            end select
          endif
          if (iprint.eq.1) write(*,*) 'nc,xc = ',nc,xc
c
c  .......transform to the local Cartesian system
          aux(1:3) = xc(1:3) - xcenter(1:3)
          call transform(12,transf,aux,xp)
          if (iprint.eq.1) then
            write(*,*) 'xp = ',xp
          endif
          call transform(12,transf,dxcdxi,dxpdxi(1:3,1))
c
c  .......compute theta,z coordinates
          call cart_to_polar(xp(1:2), void,ue(1,ie))
          if (abs(void-radius).gt.GEOM_TOL) then
            write(*,*) 'recta_conf: edges: void-radius = ',void-radius
            call pause
            iprint=1
            go to 10
          endif
          ue(2,ie) = xp(3)
c
c  .......determine derivatives
          sintheta = sin(ue(1,ie)); costheta = cos(ue(1,ie))
          if (abs(sintheta).gt.0.1d0) then
            duedxi(1,ie) = dxpdxi(1,1)/(-radius*sintheta)
          else
            duedxi(1,ie) = dxpdxi(2,1)/(radius*costheta)
          endif
          duedxi(2,ie) = dxpdxi(3,1)
c
          if (iprint.eq.1) then
            write(*,7005) ie, xc(1:3), ue(1:2,ie), duedxi(1:2,ie)
 7005       format('recta_conf: ie, xc, theta,z, duedxi = ',
     .              i2,2x,3e12.5,2x,2e12.5,3x,2e12.5)
          endif

c
c  .....end of loop through edges
        enddo
c
c  Step 3: use transfinite interpolation with linear blending functions
c       to compute the parametrization in surfaces (theta,z) coordinates
c
        u(1:2) = -( uv(1:2,1)*(1.d0-Xi(1))*(1.d0-Xi(2))
     .           +  uv(1:2,2)*      Xi(1) *(1.d0-Xi(2))
     .           +  uv(1:2,3)*      Xi(1) *      Xi(2)
     .           +  uv(1:2,4)*(1.d0-Xi(1))*      Xi(2))
     .           +  ue(1:2,1)*(1.d0-Xi(2))
     .           +  ue(1:2,2)*      Xi(1)
     .           +  ue(1:2,3)*      Xi(2)
     .           +  ue(1:2,4)*(1.d0-Xi(1))
c
c  .....formulas for the derivatives below are obtained
c       by differentiating the formula for the transfinite
c       interpolation above
        dudxi(1:2,1) = -( - uv(1:2,1)*(1.d0-Xi(2))
     .                 +    uv(1:2,2)*(1.d0-Xi(2))
     .                 +    uv(1:2,3)*      Xi(2)
     .                 -    uv(1:2,4)*      Xi(2))
     .                 +duedxi(1:2,1)*(1.d0-Xi(2))
     .                 +    ue(1:2,2)
     .                 +duedxi(1:2,3)*      Xi(2)
     .                 -    ue(1:2,4)
        dudxi(1:2,2) = -( - uv(1:2,1)*(1.d0-Xi(1))
     .                 -   uv(1:2,2)*      Xi(1)
     .                 +    uv(1:2,3)*      Xi(1)
     .                 +    uv(1:2,4)*(1.d0-Xi(1)))
     .                 -    ue(1:2,1)
     .                 +duedxi(1:2,2)*      Xi(1)
     .                 +    ue(1:2,3)
     .                 +duedxi(1:2,4)*(1.d0-Xi(1))
c
      if (iprint.eq.1) then
        write(*,7006) u,dudxi(1:2,1),dudxi(1:2,2)
 7006   format('recta_conf: u, dudxi = ',3(2f8.3,2x))
      endif
c
c  .....compute the local Cartesian coordinates and their derivatives
        sintheta = sin(u(1)); costheta = cos(u(1))
        xp(1) = radius*costheta
        xp(2) = radius*sintheta
        xp(3) = u(2)
        dxpdxi(1,1:2) = -xp(2)*dudxi(1,1:2)
        dxpdxi(2,1:2) =  xp(1)*dudxi(1,1:2)
        dxpdxi(3,1:2) =  dudxi(2,1:2)
c
c  .....switch to the global system of coordinates
        call transform(21,transf,X,xp)
        X = X + xcenter
        call transform(21,transf,Dxdxi(1:3,1),dxpdxi(1:3,1))
        call transform(21,transf,Dxdxi(1:3,2),dxpdxi(1:3,2))
c
      case default
        write(*,*) 'recta_conf: UNFINISHED !!'
        stop 1
      end select
c
      if (iprint.eq.1) then
        write(*,7010) xp, X
 7010   format('recta_conf: xp,X = ',2(3f8.3,2x))
        call pause
      endif
c
c
      end
