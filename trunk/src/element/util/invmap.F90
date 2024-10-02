!----------------------------------------------------------------------
!
!   routine name       - invmap
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - given physical coordinates Xp of a point in
!                        space and geometry dof for an element,
!                        routine determines whether the point lies within
!                        the element and returns its master element
!                        coordinates
!
!   arguments :
!     in:
!             Ntype        - element type
!             Xp           - physical coordinates of a point
!             Norder       - element order
!             Nedge_orient - edge orientations
!             Nface_orient - face orientations
!             Xnod         - element geometry dof
!     out:
!             Idec     = 1 if the point is within the element
!                        0 otherwise or NR iterations have not converged
!             Xi       - master element coordinates of Xp
!
!----------------------------------------------------------------------
!
   subroutine invmap(Ntype,Xp,Norder,Nedge_orient,Nface_orient, &
                     Xnod, Idec,Xi)
!
      use data_structure3D
      use element_data
      use control
!
      implicit none
!
      integer, intent(in)  :: Ntype,Norder(19)
      integer, intent(in)  :: Nedge_orient(12),Nface_orient(6)
      real(8), intent(in)  :: Xnod(3,MAXbrickH)
      integer, intent(out) :: Idec
      real(8), intent(out) :: Xi(3)
!
      real(8) :: Xp(3)
!
!  ...shape functions and their derivatives wrt master coordinates
      real(8) :: shapH(MAXbrickH),gradH(3,MAXbrickH)
      integer :: nrdofH
!
!  ...geometry
      real(8) :: dxdxi(3,3),dxidx(3,3)
!
!  ...increment in xi
      real(8) :: dxi(3)
!
      real(8) :: x(3),rjac,dmax
      integer :: i,iter,k
      integer :: iflag
!
#if HP3D_DEBUG
      integer :: iprint
#endif
!
      real(8), parameter :: eps = 1.d-8
!
!-----------------------------------------------------------------------
!
#if HP3D_DEBUG
      iprint=0
!
!  ...initiate the NR iterations
      if (iprint.eq.1) then
        write(*,*)'***************************'
        write(*,*)'invmap: NEWTON-RAPHSON LOOP'
      endif
#endif
!
!  ...use element centroid as initial guess
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
        call shape3DH(Ntype,Xi,Norder, Nedge_orient,Nface_orient, &
                      nrdofH,shapH,gradH)
!
!  .....determine physical coordinates and the derivatives of
!       the physical coordinates wrt master element coordinates
        x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0
        do k=1,nrdofH
          x(1:3) = x(1:3) + Xnod(1:3,k)*shapH(k)
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + Xnod(1:3,k)*gradH(i,k)
          enddo
        enddo
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,*)'----------------------------'
          write(*,7001) iter,Xi(1:3),x(1:3)
 7001     format('invmap: iter,Xi,x = ',i4,3(e12.5,2x),2x,3(e12.5,2x))
        endif
#endif
!
!  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag)
!
!  .....negative jacobian indicates that we are outside of the
!       element
        if (iflag.ne.0) then
          Idec=0 ; return
        endif
        dxi(1:3)=0.d0
        do k=1,3
          dxi(1:3)=dxi(1:3)+dxidx(1:3,k)*(Xp(k)-x(k))
        enddo
!
!  .....update xi
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7002) eps,dxi(1:3)
 7002     format('invmap: eps, dxi = ',e12.5,4x,3(e12.5,2x))
        endif
#endif
        Xi(1:3) = Xi(1:3) + dxi(1:3)
!
!  .....check convergence
        dmax = max(abs(dxi(1)),abs(dxi(2)),abs(dxi(3)))
        if (dmax.lt.eps) goto 10
      enddo
      Idec=0; return
 10   continue
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'invmap: NR LOOP FINISHED'
        write(*,*) '***************************'
        call pause
      endif
#endif
!
!  ...check if the point lies within the element
      Idec=1
      select case(Ntype)
      case(MDLP)
        if ( (Xi(1).lt.-GEOM_TOL).or.           &
             (Xi(2).lt.-GEOM_TOL).or.           &
             (Xi(1)+Xi(2).gt.1.d0+GEOM_TOL).or. &
             (Xi(3).lt.-GEOM_TOL).or.           &
             (Xi(3).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      case(MDLB)
        if ((Xi(1).lt.-GEOM_TOL).or.(Xi(1).gt.1.d0+GEOM_TOL).or.  &
            (Xi(2).lt.-GEOM_TOL).or.(Xi(2).gt.1.d0+GEOM_TOL).or.  &
            (Xi(3).lt.-GEOM_TOL).or.(Xi(3).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      case(MDLN)
        if ((Xi(1).lt.-GEOM_TOL).or.   &
            (Xi(2).lt.-GEOM_TOL).or.   &
            (Xi(3).lt.-GEOM_TOL).or.   &
            ((Xi(1)+Xi(2)+Xi(3)).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      case(MDLD)
        if ((Xi(1).lt.-GEOM_TOL).or.               &
            (Xi(2).lt.-GEOM_TOL).or.               &
            (Xi(3).lt.-GEOM_TOL).or.               &
            ((Xi(1)+Xi(3)).gt.1.d0+GEOM_TOL).or.   &
            ((Xi(2)+Xi(3)).gt.1.d0+GEOM_TOL)) then
          Idec=0;
        endif
      end select
!
#if HP3D_DEBUG
      if(iprint.eq.1)then
        write(*,*)'Exiting invmap with Idec=',Idec
      endif
#endif
!
   end subroutine invmap
