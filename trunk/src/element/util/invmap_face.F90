#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - invmap_face
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - given an element face and a point, both
!                        projected onto a plane, find the correspon-
!                        ding master element coordinates for a point
!                        on the element face such that the projection
!                        of the point coincides with the given point
!
!   arguments :
!     in:
!                Mdle  - An element (middle node)
!                Ifc   - face number
!                Xy_in - coordinates of a projected point
!     out
!                T     - master face coordinates of the point
!                X     - physical coordinates of the point (before
!                        the projection)
!                Ierr  = 0, if the procedure has converged
!                        1  otherwise
!
!----------------------------------------------------------------------
!
   subroutine invmap_face(Mdle,Ifc,Xy_in, T,X,Ierr)
!
      use graphmod
      use data_structure3D
!
      implicit none
!
      integer, intent(in)  :: Mdle,Ifc
      real(8), intent(in)  :: Xy_in(2)
      real(8), intent(out) :: T(2),X(3)
      integer, intent(out) :: Ierr
!
!  ...geometry dof
      real(8) :: xnod(3,MAXbrickH)
!
!  ...orientation for element nodes, order
      integer :: nedge_orient(12),nface_orient(6),norder(19)
!
!  ...geometry
      real(8) :: xi(3),dxdxi(3,3), &
                 dt(2),dXidt(3,2),dxdt(3,2),xp(3),dxpdt(3,2)
!
!  ...shape functions
      real(8) :: vshapH(MAXbrickH),dvshapH(3,MAXbrickH)
!
!  ...miscellanea
      real(8) :: s,det,det1,det2
      integer :: nrdofH
      integer :: iter,i,j,k
!
!  ...tolerance for Newton-Raphson iterations
      real(8), parameter :: epsilon = 1.d-4
!
#if HP3D_DEBUG
      integer :: ivar
      integer :: iprint
      iprint = 0
#endif
!
!  ...find order of approXimation for the element
      call find_order(Mdle, norder)
!
!  ...determine geometry dof
      call nodcor(Mdle, xnod)
!
!  ...determine node orientation
      call find_orient(Mdle, nedge_orient,nface_orient)
!
!  ...use origin for the initial guess
      T(1:2) = 0.d0
!
!  ...start NR iterations
      do iter=1,10
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7001) iter
 7001     format('invmap_face: iter = ',i2)
        endif
#endif
!
!  .....determine master element coordinates of the point and their
!       derivatives wrt master face coordinates
        call face_param(NODES(Mdle)%ntype,Ifc,t, Xi,dxidt)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7002) S_Type(NODES(Mdle)%ntype),Ifc,T
 7002     format('invmap_face: NODES(Mdle)%type,Ifc,T = ', &
                 a5,i2,2x,2f8.3)
          do ivar=1,3
            write(*,7003) ivar,Xi(ivar),dxidt(ivar,1:2)
 7003       format('ivar,Xi,dxidt = ',i2,f8.3,2x,2f8.3)
          enddo
        endif
#endif
!
!  .....evaluate shape functions at the point
        call shape3DH(NODES(Mdle)%ntype, &
                      Xi,norder,nedge_orient,nface_orient, &
                      nrdofH,vshapH,dvshapH)
!
!  .....compute physical coordinates and their derivatives wrt
!       master element coordinates
        X(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0
        do k=1,nrdofH
          X(1:3) = X(1:3) + xnod(1:3,k)*vshapH(k)
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dvshapH(i,k)
          enddo
        enddo
#if HP3D_DEBUG
        if (iprint.eq.1) then
          do ivar=1,3
            write(*,7004) ivar,X(ivar),dxdxi(ivar,1:3)
 7004       format('ivar,X,dxdxi = ',i2,f8.3,2x,3f8.3)
          enddo
        endif
#endif
!
!  .....compute derivatives of the physical coordinates wrt master face
!       coordinates
        dxdt(1:3,1:2) = 0.d0
        do i=1,2
          do j=1,3
            dxdt(1:3,i) = dxdt(1:3,i) + dxdxi(1:3,j)*dxidt(j,i)
          enddo
        enddo
#if HP3D_DEBUG
        if (iprint.eq.1) then
          do ivar=1,3
            write(*,7005) ivar,dxdt(ivar,1:2)
 7005       format('ivar,dxdt = ',i2,2f8.3)
          enddo
        endif
#endif
!
!  .....compute the coordinates of the projected point and their
!       derivatives wrt master face coordinates
        xp(1:3)=0.d0; dxpdt(1:3,1:2) = 0.d0
        do j=1,3
          xp(1:3) = xp(1:3)+RMTR(1:3,j)*X(j)
        enddo
        do i=1,2
          do j=1,3
            dxpdt(1:3,i) = dxpdt(1:3,i)+RMTR(1:3,j)*dxdt(j,i)
          enddo
        enddo
!
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,*) 'invmap_face: dxpdt,Xy_in-xp = '
          do i=1,2
            write(*,7006) dxpdt(i,1:2),Xy_in(i)-xp(i)
 7006       format(2e12.5, 3x,e12.5)
          enddo
        endif
#endif
!  .....solve for the next iterate
        det  = dxpdt(1,1)*dxpdt(2,2) - dxpdt(1,2)*dxpdt(2,1)
        if (abs(det).lt.1.0d-10) then
          Ierr=1; return
        endif
        det1 = (Xy_in(1)-xp(1))*dxpdt(2,2) - (Xy_in(2)-xp(2))*dxpdt(1,2)
        det2 = (Xy_in(2)-xp(2))*dxpdt(1,1) - (Xy_in(1)-xp(1))*dxpdt(2,1)
        dt(1) = det1/det; dt(2) = det2/det
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7007) dt
 7007     format('invmap_face: dt = ',2f8.3)
        endif
#endif
!
!  .....update
        T(1:2) = T(1:2) + dt(1:2)
        s = max(abs(dt(1)),abs(dt(2)))
        if (s.lt.epsilon) goto 20
      enddo
!
      Ierr=1; return
!
 20   Ierr=0; return
!
!
   end subroutine invmap_face

#endif
