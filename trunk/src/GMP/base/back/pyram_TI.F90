!-----------------------------------------------------------------------
!
!   routine name       - pyram_TI
!
!-----------------------------------------------------------------------
!
!   latest revision    - Dec 08
!
!   purpose            - transfinite interpolation pyramid
!
!   arguments :
!     in:
!               No     - a GMP pyramid number
!               Eta    - reference coordinates of a point
!                        in the triangle
!     out:
!               Xp      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the parameters
!
!-----------------------------------------------------------------------
!
      subroutine pyram_TI(No,Eta, Xp,Dxdeta)
!
      use control
      use GMP
      use element_data
#include "syscom.blk"
!
      dimension Eta(3),Xp(3),Dxdeta(3,3)
!
!  ...pyramid element order and shape functions
      dimension norder(14),shapH(5),dshapH(3,5)
!
!  ...blending functions
      dimension dblend(3),blend_edge(4),dblend_edge(3,4)
!
!  ...projections
      dimension dtedeta(3),tf(2),dtfdeta(2,3)
!
!  ...edge kernels
      dimension xe(3),dxedt(3)
!
!  ...face kernels
      dimension xf(3),dxfdtf(3,2)
!
!  ...work space
      dimension nvoid(12)
!
      double precision x,y,z,xz1,yz1,z1
!
!----------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) No,Eta(1:3)
 7001   format('pyram_TI: No,Eta = ',i5,2x,3e12.5)
      endif
!
      x = Eta(1); y = Eta(2); z = Eta(3)
      if (abs(z-1.d0).lt.1.d-12) z = 1.d0-1.d-12
      xz1 = 1.d0-x-z
      yz1 = 1.d0-y-z
      z1 = 1.d0-z
!
!  ...compute vertex shape functions
      norder(1:14)=1; norder(9)=11; nvoid(1:12)=0
      call shapeHd(Eta,norder,nvoid,nvoid, nrdofH,shapH,dshapH)
!
!  ...start with linear interpolant
      Xp(1:3) = 0.d0; Dxdeta(1:3,1:3) = 0.d0
      do iv=1,5
        np = PYRAMIDS(No)%VertNo(iv)
        Xp(1:3) = Xp(1:3) + POINTS(np)%Rdata(1:3)*shapH(iv)
        do ivar=1,3
          Dxdeta(ivar,1:3) = Dxdeta(ivar,1:3)                    &
                        + POINTS(np)%Rdata(ivar)*dshapH(1:3,iv)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'pyram_TI: VERTEX INTERPOLANT = '
        do ivar=1,3
          write(*,7011) Xp(ivar),Dxdeta(ivar,1:3)
 7011     format(e12.5,3x,3e12.5)
        enddo
      endif
!
!-----------------------------------------------------------------------
!
!  ...edge contributions
      do ie=1,8
!
        nc = PYRAMIDS(No)%EdgeNo(ie); norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
        if (CURVES(nc)%Type.eq.'Seglin') cycle
!
!  .....project x,y,z onto the edge
        call proj_d2e(Eta,ie,ShapH,DshapH, te,dtedeta)
        if ((abs(te).lt.GEOM_TOL).or.(abs(1.d0-te).lt.GEOM_TOL)) cycle
        if (iprint.eq.1) then
          write(*,7012) ie,nc,CURVES(nc)%Type
 7012     format('pyram_TI: ie,nc,Type = ',i2,i5,2x,a5)
        endif
!
!  .....blending function
        select case(ie)
        case(1)
          blend = x*xz1*yz1/z1
          dblend(1) = (xz1-x)*yz1/z1
          dblend(2) = -x*xz1/z1
          dblend(3) = x*(-yz1-xz1)/z1 + x*xz1*yz1/z1**2
        case(2)
          blend = x*y*yz1/z1
          dblend(1) = y*yz1/z1
          dblend(2) = x*(yz1-y)/z1
          dblend(3) = -x*y/z1 + x*y*yz1/z1**2
        case(3)
          blend = x*xz1*y/z1
          dblend(1) = (xz1-x)*y/z1
          dblend(2) = x*xz1/z1
          dblend(3) = -x*y/z1 + x*xz1*y/z1**2
        case(4)
          blend = xz1*y*yz1/z1
          dblend(1) = -y*yz1/z1
          dblend(2) = xz1*(yz1-y)/z1
          dblend(3) = y*(-yz1-xz1)/z1 + xz1*y*yz1/z1**2
        case(5,6,7,8)
          iv = ie-4
          blend = ShapH(iv)*z
          dblend(1:2) = DshapH(1:2,iv)*z
          dblend(3) = DshapH(3,iv)*z + ShapH(iv)
        end select
        if (ie.le.4) then
          blend_edge(ie) = blend
          dblend_edge(1:3,ie) = dblend(1:3)
        endif
!
!  .....add edge contribution
        Xp(1:3) = Xp(1:3) + xe(1:3)*blend
        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar)                  &
                           + dxedt(1:3)*dtedeta(ivar)*blend    &
                           + xe(1:3)*dblend(ivar)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'pyram_TI: AFTER EDGES = '
        do ivar=1,3
          write(*,7011) Xp(ivar),Dxdeta(ivar,1:3)
        enddo
      endif
!
!-----------------------------------------------------------------------
!
!  ...bottom face
      if=1
      call decode(PYRAMIDS(No)%FigNo(if), nr,norient)
      if ((RECTANGLES(nr)%Type .eq. 'BilQua') .or. &
          (RECTANGLES(nr)%Type .eq. 'TraQua') .or. &
          (RECTANGLES(nr)%Type .eq. 'PTIRec')) go to 20
      if (iprint.eq.1) then
        write(*,7014) if,nr,RECTANGLES(nr)%Type
 7014   format('pyram_TI: if,nr,Type = ',i2,i5,2x,a5)
      endif
!
!  ...blending function
      blend = ShapH(1)*ShapH(3)
      dblend(1:3) = DshapH(1:3,1)*ShapH(3) + ShapH(1)*DshapH(1:3,3)
!
!  ...project onto the bottom face
      xz = x/z1; yz = y/z1
!
!  ...finish
      write(*,*) 'pyram_TI: unfinished'; stop
 20   continue
!
!-----------------------------------------------------------------------
!
!  ...lateral mid-face nodes shape functions
      do if=2,5
        call decode(PYRAMIDS(No)%FigNo(if), nt,norient)
        if ((TRIANGLES(nt)%Type.eq.'TransTri') .or. &
            (TRIANGLES(nt)%Type.eq.'PlaneTri')) cycle
!
!  .....project xyz onto the face
        call proj_d2f(Eta,if, tf,dtfdeta)
!
!  .....if point is on the edge, then the bubble contribution is zero
        if ((abs(tf(2)).lt.GEOM_TOL).or.(abs(tf(1)).lt.GEOM_TOL) .or. &
            (abs(1.d0-tf(1)-tf(2)).lt.GEOM_TOL)) cycle
        if (iprint.eq.1) then
          write(*,7013) if,nt,TRIANGLES(nt)%Type
 7013     format('pyram_TI: if,nt,Type = ',i2,i5,2x,a5)
        endif
!
!  .....compute the face kernel
        call trianK(nt,tf,norient, xf,dxfdtf)
!
!  .....blending function
        blend = blend_edge(if-1)*z
        dblend(1:3) = dblend_edge(1:3,if-1)*z
        dblend(3) = dblend(3) + blend_edge(if-1)
!
!  .....add face contribution
        Xp(1:3) = Xp(1:3) + xf(1:3)*blend
        do ivar=1,3
          Dxdeta(1:3,ivar) = Dxdeta(1:3,ivar)                       &
                           + (dxfdtf(1:3,1)*dtfdeta(1,ivar)         &
                           + dxfdtf(1:3,2)*dtfdeta(2,ivar))*blend   &
                           + xf(1:3)*dblend(ivar)
        enddo
      enddo
!
!
      end subroutine pyram_TI
