c----------------------------------------------------------------------
c
c   routine name       - invmap
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - given physical coordinates Xp of a point in
c                        space and geometry dof for an element,
c                        routine determines whether the point lies within
c                        the element and returns its master element
c                        coordinates
c
c   arguments :
c     in:
c             Ntype    - element type
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
      subroutine invmap(Ntype,Xp,Norder,Nedge_orient,Nface_orient,
     .                  Xnod, Idec,Xi)
c
      use data_structure3D
      use element_data
      use control
#include "syscom.blk"
      common /cinvmap/ iprint_invmap
c
      integer, intent(in)  :: Ntype,Norder(19)
      integer, intent(in)  :: Nedge_orient(12),Nface_orient(6)
      real(8), intent(in)  :: Xnod(3,MAXbrickH)
      integer, intent(out) :: Idec
      real(8), intent(out) :: Xi(3)
c
      real(8) :: Xp(3)
c
c  ...shape functions and their derivatives wrt master coordinates
      real(8) :: shapH(MAXbrickH),gradH(3,MAXbrickH)
      integer :: nrdofH
c
c  ...geometry
      real(8) :: dxdxi(3,3),dxidx(3,3)
c
c  ...increment in xi
      real(8) :: dxi(3)
c
      real(8) :: x(3)
      integer :: i,iter,k
      integer :: iflag,iprint
c
      real(8), parameter :: eps = 1.d-8
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
      Xi(1:3)=0.d0
      do i=1,nvert(Ntype)
        select case(Ntype)
        case(MDLB,BRIC) ; Xi(1:3) = Xi(1:3) + BRICK_COORD(1:3,i)
        case(MDLN,TETR) ; Xi(1:3) = Xi(1:3) + TETRA_COORD(1:3,i)
        case(MDLP,PRIS) ; Xi(1:3) = Xi(1:3) + PRISM_COORD(1:3,i)
        case(MDLD,PYRA) ; Xi(1:3) = Xi(1:3) + PYRAM_COORD(1:3,i)
        endselect
      enddo
      Xi(1:3)=Xi(1:3)/nvert(Ntype)

      do iter=1,15
        call shape3DH(Ntype,Xi,Norder, Nedge_orient,Nface_orient,
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
          write(*,7001) iter,Xi(1:3),x(1:3)
 7001     format('invmap: iter,Xi,x = ',i4,3(e12.5,2x),2x,3(e12.5,2x))
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
        Xi(1:3) = Xi(1:3) + dxi(1:3)
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
      select case(Ntype)
      case(MDLP)
        if ( (Xi(1).lt.-GEOM_TOL).or.
     .       (Xi(2).lt.-GEOM_TOL).or.
     .       (Xi(1)+Xi(2).gt.1.d0+GEOM_TOL).or.
     .       (Xi(3).lt.-GEOM_TOL).or.
     .       (Xi(3).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      case(MDLB)
        if ((Xi(1).lt.-GEOM_TOL).or.(Xi(1).gt.1.d0+GEOM_TOL).or.
     .      (Xi(2).lt.-GEOM_TOL).or.(Xi(2).gt.1.d0+GEOM_TOL).or.
     .      (Xi(3).lt.-GEOM_TOL).or.(Xi(3).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      case(MDLN)
        if ((Xi(1).lt.-GEOM_TOL).or.(Xi(2).lt.-GEOM_TOL).or.
     .      (Xi(3).lt.-GEOM_TOL).or.
     .      ((Xi(1)+Xi(2)+Xi(3)).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      case(MDLD)
        if ((Xi(1).lt.-GEOM_TOL).or.
     .      (Xi(2).lt.-GEOM_TOL).or.
     .      (Xi(3).lt.-GEOM_TOL).or.((Xi(1)+Xi(3)).gt.1.d0+GEOM_TOL).or.
     .      ((Xi(2)+Xi(3)).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      end select
c
      if(iprint.eq.1)then
        write(*,*)'Exiting invmap with Idec=',Idec
      endif
c
      end subroutine invmap
