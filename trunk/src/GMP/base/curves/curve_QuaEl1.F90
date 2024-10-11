!----------------------------------------------------------------------
!
!   routine name       - curve_QuaEl1
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 18
!
!   purpose            - routine defines a parametrization for
!                        a quarter of ellipse between 0 and pi/2
!
!   arguments :
!     in:
!               No     - the curve number
!               Eta    - reference coordinate  (between 0 and 1)
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        reference coordinate
!
!----------------------------------------------------------------------
   subroutine curve_QuaEl1(No,Eta, X,Dxdeta)
!
      use GMP
!
      implicit none
!
      integer,                 intent(in)  :: No
      real(8),                 intent(in)  :: Eta
      real(8), dimension(3),   intent(out) :: X,Dxdeta
!
      real(8), dimension(3)   :: aG,bG,center,aC,bC,xL,dxLdeta
      real(8), dimension(3,3) :: rotM
      real(8) :: rx,ry,tt,pihalf
!
#if HP3D_DEBUG
      integer :: i,iprint
      iprint=0
#endif
!----------------------------------------------------------------------
!
      if (CURVES(No)%Type.ne.'QuaEl1') then
        write(*,7001)
 7001   format(' curve_QuaEl1: WRONG CALL')
        stop
      endif
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) No,Eta
 7002   format(' curve_QuaEl1: No,Eta = ',i4,2x,f8.3)
      endif
#endif
!
!  ...get the endpoints physical (global) coordinates
      call pointr(CURVES(No)%EndPoNo(1), aG)
      call pointr(CURVES(No)%EndPoNo(2), bG)
!
!  ...get coordinates of the center and relevant powers
      center(1:3)=CURVES(No)%Rdata(1:3)
!
!  ...get relative endpoints
      aC(1:3)=aG(1:3)-center(1:3)
      bC(1:3)=bG(1:3)-center(1:3)
!
!  ...calculate the semiaxes radii using the fact that the points
!     are the semiaxes of the superellipse
      call norm(aC, rx)
      call norm(bC, ry)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) ' '
        write(*,7011) aG(1:3)
 7011   format(' aG  = ',3(e12.5,2x))
        write(*,7012) bG(1:3)
 7012   format(' bG  = ',3(e12.5,2x))
        write(*,7008) center(1:3)
 7008   format(' center  = ',3(e12.5,2x))
        write(*,7009) rx,ry
 7009   format(' rx,ry = ',2(e12.5,2x))
      endif
#endif
!
!  ...now calculate the coordinates in 2D x-y plane
!  ...first rescale Eta to between 0 and pi/2
      pihalf = acos(0.d0)
      tt = Eta*pihalf
!  ...next evaluate the superellipse (and eta derivative)
      xL(1) = rx*cos(tt)
      xL(2) = ry*sin(tt)
      xL(3) = 0.d0
      dxLdeta(1) = pihalf*rx*(-sin(tt))
      dxLdeta(2) = pihalf*ry*(+cos(tt))
      dxLdeta(3) = 0.d0
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) ' '
        write(*,7013) tt
 7013   format(' tt  = ',e12.5)
        write(*,7006) xL(1:3)
 7006   format(' xL  = ',3(e12.5,2x))
        write(*,7007) dxLdeta(1:3)
 7007   format(' dxLdeta = ',3(e12.5,2x))
      endif
#endif
!
!  ...the idea is to map these coordinates to the global physical
!     coordinates
!  ...this is a rotation (about the origin) plus a trivial translation
!                       X=rotM*xL+center
!  ...first calculate the rotation (the hard part) by using that:
!                  aL=[rx;0;0],    bL=[0;ry;0],
!           rotM*aL(1:3)      = aC(1:3)
!           rotM*bL(1:3)      = bC(1:3)
!           rotM*cross(aL,bL) = cross(aC,bC)
!     after some algebra, it follows
!     rotM=[(1/rx)*aC,(1/ry)*bC,(1/(rx*ry))*cross(aC,bC)]
      rotM(1:3,1) = (1.d0/rx)*aC(1:3)
      rotM(1:3,2) = (1.d0/ry)*bC(1:3)
      call cross_product(aC(1:3),bC(1:3), rotM(1:3,3))
      rotM(1:3,3) = (1.d0/(rx*ry))*rotM(1:3,3)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        do i=1,3
          write(*,7005) i,rotM(i,1:3)
 7005     format(' i,rotM(i,:) = ',i1,2x,3(e12.5,2x))
        enddo
      endif
#endif
!
!  ...finally compute X and Dxdeta
      X(1:3) = rotM(1:3,1)*xL(1) + rotM(1:3,2)*xL(2) + &
               rotM(1:3,3)*xL(3) + center(1:3)
      Dxdeta(1:3) = rotM(1:3,1)*dxLdeta(1) + rotM(1:3,2)*dxLdeta(2) + &
                    rotM(1:3,3)*dxLdeta(3)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) ' '
        write(*,7003) X(1:3)
        write(*,7004) Dxdeta(1:3)
 7003   format(' X      = ',3(e12.5,2x))
 7004   format(' Dxdeta = ',3(e12.5,2x))
      endif
#endif
!
   end subroutine curve_QuaEl1
