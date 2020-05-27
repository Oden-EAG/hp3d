c----------------------------------------------------------------------
c
c   routine name       - invmap
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 09
c
c   purpose            - given physical coordinates Xp of a point in 
c                        space and geometry dof for an element,
c                        routine determines whether the point lies within
c                        the element and returns its master element
c                        coordinates
c
c   arguments :
c     in:
c             Type     - element type
c             Xp       - physical coordinates of a point
c             Norder   - element order
c             Nedge_orient, Nface_orient - edge and face orientations
c             Xnod     - element geometry dof
c     out:
c             Idec     = 1 if the point is within the element
c                        0 otherwise or NR iterations have not converged
c             Xi       - master element coordinates of Xp
c
c----------------------------------------------------------------------
c
      subroutine invmap(Type,Xp,Norder,Nedge_orient,Nface_orient,
     .                  Xnod, Idec,Xi)
c  
      use data_structure3D
      use element_data
      use control
#include "syscom.blk"
      common /cinvmap/ iprint_invmap
c
      character(len=4) :: Type
      dimension Norder(19),Nedge_orient(12),Nface_orient(6)
      dimension Xp(3),Xnod(3,MAXbrickH)
c
c  ...shape functions and their derivatives wrt master coordinates
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
c
c  ...geometry
      dimension xi(3),dxdxi(3,3),dxidx(3,3)
c
c  ...increment in xi
      dimension dxi(3)
c
      dimension x(3)
      parameter (eps = 1.d-8)
c
c
c-----------------------------------------------------------------------
c
      iprint=iprint_invmap
ccc      iprint=0
c
c  ...initiate the NR iterations
      if (iprint.eq.1) then
        write(*,*)'***************************'
        write(*,*)'invmap: NEWTON-RAPHSON LOOP'
      endif
c
c  ...use element centroid as initial guess
      xi(1:3)=0.d0
      do i=1,nvert(Type)
        select case(Type) 
        case('mdln','tetr') ; xi(1:3) = xi(1:3) + TETRA_COORD(1:3,i)
        case('mdlb','bric') ; xi(1:3) = xi(1:3) + BRICK_COORD(1:3,i)
        case('mdlp','pris') ; xi(1:3) = xi(1:3) + PRISM_COORD(1:3,i)
        case('mdld','pyra') ; xi(1:3) = xi(1:3) + PYRAM_COORD(1:3,i)
        endselect
      enddo
      xi(1:3)=xi(1:3)/nvert(Type)

      do iter=1,15
        call shape3DH(Type,xi,Norder, Nedge_orient,Nface_orient,
     .                nrdofH,shapH,gradH)
c
c  .....determine physical coordinates and the derivatives of 
c       the physical coordinates wrt master element coordinates
        x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
        do k=1,nrdofH
          x(1:3) = x(1:3) + Xnod(1:3,k)*shapH(k)
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + Xnod(1:3,k)*gradH(i,k)
          enddo
        enddo
        if (iprint.eq.1) then
          write(*,*)'----------------------------'
          write(*,7001) iter,xi(1:3),x(1:3)
 7001     format('invmap: iter,xi,x = ',i4,3(e12.5,2x),2x,3(e12.5,2x))
        endif
c
c  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag)
c
c  .....negative jacobian indicates that we are outside of the
c       element
        if (iflag.ne.0) then
          Idec=0 ; return
        endif
        dxi(1:3)=0.d0
        do k=1,3
          dxi(1:3)=dxi(1:3)+dxidx(1:3,k)*(Xp(k)-x(k))
        enddo
c
c  .....update xi
        if (iprint.eq.1) then
          write(*,7002) eps,dxi(1:3)
 7002     format('invmap: eps, dxi = ',e12.5,4x,3(e12.5,2x))
        endif
        xi(1:3) = xi(1:3) + dxi(1:3)
c
c  .....check convergence
        dmax = max(abs(dxi(1)),abs(dxi(2)),abs(dxi(3)))
        if (dmax.lt.eps) goto 10
      enddo
      Idec=0; return
 10   continue
      if (iprint.eq.1) then
        write(*,*) 'invmap: NR LOOP FINISHED'
        write(*,*) '***************************'
        call pause
      endif

c
c  ...check if the point lies within the element
      Idec=1
      select case(Type)
      case('mdlp')
        if ( (xi(1).lt.-GEOM_TOL).or.
     .       (xi(2).lt.-GEOM_TOL).or.
     .       (xi(1)+xi(2).gt.1.d0+GEOM_TOL).or.
     .       (xi(3).lt.-GEOM_TOL).or.
     .       (xi(3).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      case('mdlb')
        if ((xi(1).lt.-GEOM_TOL).or.(xi(1).gt.1.d0+GEOM_TOL).or.
     .      (xi(2).lt.-GEOM_TOL).or.(xi(2).gt.1.d0+GEOM_TOL).or.
     .      (xi(3).lt.-GEOM_TOL).or.(xi(3).gt.1.d0+GEOM_TOL)) then
          Idec=0; 
        endif
      case('mdln')
        if ((xi(1).lt.-GEOM_TOL).or.(xi(2).lt.-GEOM_TOL).or.
     .      (xi(3).lt.-GEOM_TOL).or.
     .      ((xi(1)+xi(2)+xi(3)).gt.1.d0+GEOM_TOL)) then
          Idec=0; 
        endif
      case('mdld')
        if ((xi(1).lt.-GEOM_TOL).or.
     .      (xi(2).lt.-GEOM_TOL).or.
     .      (xi(3).lt.-GEOM_TOL).or.((xi(1)+xi(3)).gt.1.d0+GEOM_TOL).or.
     .      ((xi(2)+xi(3)).gt.1.d0+GEOM_TOL)) then
          Idec=0; 
        endif
      end select
c
      if(iprint.eq.1)then
        write(*,*)'Exiting invmap with Idec=',Idec
      endif
c     
      end
