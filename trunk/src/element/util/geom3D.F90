!----------------------------------------------------------------------
!
!   routine name       - geom3D
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - element geometry routine
!
!   arguments :
!     in:
!             Mdle     - element middle node
!             Xi       - master element coordinates of a point
!             Xnod     - element geometry dof
!             ShapH    - values of shape functions at a point
!             GradH    - derivatives of the shape functions wrt
!                        master coordinates at the point
!             NrdofH   - number of element H1 dof
!     out:
!             X        - physical coordinates of the point
!             Dxdxi    - Jacobian matrix
!             Dxidx    - inverse Jacobian matrix
!             Rjac     - Jacobian (determinant of the Jacobian matrix)
!             Iflag    = 0 OK
!                        1 negative Jacobian
!
!----------------------------------------------------------------------
!
   subroutine geom3D(Mdle,Xi,Xnod,ShapH,GradH,NrdofH, &
                     X,Dxdxi,Dxidx,Rjac,Iflag)
!
      use control
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Mdle,NrdofH
      real(8), intent(out) :: Rjac
      integer, intent(out) :: Iflag
!
      real(8) :: Xi(3),Xnod(3,MAXbrickH),             &
                 ShapH(MAXbrickH),GradH(3,MAXbrickH), &
                 X(3),Dxdxi(3,3),Dxidx(3,3)
!
      integer :: i,k
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) Mdle,Xi(1:3),EXGEOM
 7001   format('geom3D: Mdle,Xi = ',i6,3x,3f8.3,' EXGEOM = ',i3)
      endif
#endif
!
!  ...determine physical coordinates and the derivatives of
!     the physical coordinates wrt master element coordinates
      select case(EXGEOM)
!
!  ...isoparametric element
      case(0)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,*) 'geom3D: Xnod  = '
          do i=1,3
            write(*,7010) Xnod(i,1:NrdofH)
 7010       format(10(27f7.2,1x,/))
          enddo
          write(*,*) '        ShapH = '
          write(*,7010) ShapH(1:NrdofH)
          write(*,*) '        GradH = '
          do i=1,3
            write(*,7010) GradH(i,1:NrdofH)
          enddo
        endif
#endif
        X(1:3)=0.d0 ; Dxdxi(1:3,1:3)=0.d0
        do k=1,NrdofH
          X(1:3) = X(1:3) + Xnod(1:3,k)*ShapH(k)
          do i=1,3
            Dxdxi(1:3,i) = Dxdxi(1:3,i) + Xnod(1:3,k)*GradH(i,k)
          enddo
        enddo
!
!  ...exact geometry element
      case(1)
        call exact_geom(Mdle,Xi, X,Dxdxi)
      end select
!
!  ...evaluate the inverse Jacobian of the geometry map
      call geom(Dxdxi, Dxidx,Rjac,Iflag)
!
      if (Iflag .ne. 0) then
        !$OMP CRITICAL
        write(*,*) 'geom3D: negative Jacobian!'
        write(*,*) 'Mdle  = ', Mdle
        write(*,*) 'Rjac  = ', Rjac
        write(*,*) 'X     = '
        write(*,7002) X(1:3)
        write(*,*) 'Dxdxi = '
        do i=1,3
          write(*,7002) Dxdxi(i,1:3)
        enddo
        write(*,*) 'Dxidx = '
        do i=1,3
          write(*,7002) Dxidx(i,1:3)
        enddo
 7002   format(3f10.3)
        !$OMP END CRITICAL
      endif
!
   end subroutine geom3D
!
!----------------------------------------------------------------------
!
!   routine name       - bgeom3D
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - element boundary geometry routine
!
!   arguments :
!     in:
!             Mdle     - element middle node
!             Xi       - master element coordinates of a point
!             Xnod     - element geometry dof
!             ShapH    - values of shape functions at a point
!             GradH    - derivatives of the shape functions at the
!                        point
!             NrdofH   - number of element H1 dof
!             Dxidt    - derivatives of a point on the element boundary
!                        wrt face parameters
!             Nsign    - sign factor indicating the actual direction of the
!                        normal wrt to the local one (implied by the
!                        face parametrization)
!     out:
!             X        - physical coordinates of the point
!             Dxdxi    - Jacobian matrix
!             Dxidx    - inverse Jacobian matrix
!             Rjac     - Jacobian (determinant of the Jacobian matrix)
!             Dxdt     - derivatives wrt the face parameters
!             Rn       - the normal unit vector
!             Bjac     - boundary Jacobian (norm of Dxdt)
!
!----------------------------------------------------------------------
!
   subroutine bgeom3D(Mdle,Xi,Xnod,ShapH,GradH,NrdofH,Dxidt,Nsign, &
                      X,Dxdxi,Dxidx,Rjac,Dxdt,Rn,Bjac)
!
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Mdle,NrdofH,Nsign
      real(8), intent(out) :: Rjac,Bjac
!
      real(8) :: Xi(3),Xnod(3,MAXbrickH),             &
                 ShapH(MAXbrickH),GradH(3,MAXbrickH), &
                 Dxidt(3,2),Dxdxi(3,3),Dxidx(3,3),    &
                 X(3),Dxdt(3,2),Rn(3)
!
      integer :: i,j,iflag
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      call geom3D(Mdle,Xi,Xnod,ShapH,GradH,NrdofH, &
                  X,Dxdxi,Dxidx,Rjac,iflag)
      if (iflag.ne.0) then
        write(*,5999) Mdle,Rjac
 5999   format('bgeom3D: Negative Jacobian. Mdle,Rjac=',i10,',',e12.5)
        stop
      endif
!
!  ...derivatives of geometry map wrt face parameters
      Dxdt(1:3,1:2)=0.d0
      do i=1,2
        do j=1,3
          Dxdt(1:3,i) = Dxdt(1:3,i) + Dxdxi(1:3,j)*Dxidt(j,i)
        enddo
      enddo
!
!  ...surface normal
      call cross_product(Dxdt(1:3,1),Dxdt(1:3,2), Rn(1:3))
!
!  ...surface Jacobian
      call norm(Rn(1:3), Bjac)
!
!  ...unit normal vector
      Rn(1:3) = Rn(1:3)*Nsign/Bjac
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) X,Rn,Rjac,Bjac
 7001   format('bgeom3D: X,Rn,Rjac,Bjac = ',4(3f8.3,2x))
        call pause
      endif
#endif
!
   end subroutine bgeom3D
!
!
!----------------------------------------------------------------------
!
!   routine name       - refgeom3D
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - element reference geometry routine
!
!   arguments :
!     in:
!             Mdle     - element middle node
!             Xi       - master element coordinates of a point
!             Etav     - element reference vertex coordinates
!             ShapH    - values of shape functions at a point
!             GradH    - derivatives of the shape functions wrt
!                        master coordinates at the point
!             NrdofH   - number of (vertex) shape functions
!     out:
!             Eta      - reference coordinates of the point
!             Detadxi  - Jacobian matrix
!             Dxideta  - inverse Jacobian matrix
!             Rjac     - Jacobian (determinant of the Jacobian matrix)
!             Iflag    = 0 OK
!                        1 negative Jacobian
!
!----------------------------------------------------------------------
!
   subroutine refgeom3D(Mdle,Xi,Etav,ShapH,GradH,NrdofH, &
                        Eta,Detadxi,Dxideta,Rjac,Iflag)
!
      use control
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Mdle,NrdofH
      real(8), intent(out) :: Rjac
      integer, intent(out) :: Iflag
!
      real(8) :: Xi(3),Etav(3,8),      &
                 ShapH(8),GradH(3,8),  &
                 Eta(3),Detadxi(3,3),Dxideta(3,3)
!
      integer :: k,imas
!
#if HP3D_DEBUG
      integer :: i
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) Mdle,Xi(1:3)
 7001   format('refgeom3D: Mdle,Xi = ',i6,3x,3f8.3)
      endif
#endif
!
!  ...determine reference coordinates and the derivatives of
!     the reference coordinates wrt master element ccordinates
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'refgeom3D: Etav,ShapH = '
        do i=1,3
          write(*,7010) Etav(i,1:NrdofH)
 7010     format(10(27f7.2,1x,/))
        enddo
        write(*,*) '        GradH = '
        do i=1,3
          write(*,7010) GradH(i,1:NrdofH)
        enddo
      endif
#endif
      Eta(1:3) = 0.d0
      do k=1,NrdofH
        Eta(1:3) = Eta(1:3) + Etav(1:3,k)*ShapH(k)
      enddo
      do imas=1,3
        Detadxi(1:3,imas) = 0.d0
        do k=1,NrdofH
          Detadxi(1:3,imas) = &
          Detadxi(1:3,imas) + Etav(1:3,k)*GradH(imas,k)
        enddo
      enddo
!
!  ...evaluate the inverse Jacobian of the geometry map
      call geom(Detadxi, Dxideta,Rjac,Iflag)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'Eta,Detadxi = '
        do i=1,3
          write(*,7002) Eta(i),Detadxi(i,1:3)
 7002     format(f8.3,3x,3f8.3)
        enddo
        write(*,*) 'Rjac, Dxideta = ',Rjac
        do i=1,3
          write(*,7003) Dxideta(i,1:3)
 7003     format(11x,3f8.3)
        enddo
        call pause
      endif
#endif
!
   end subroutine refgeom3D
!
!
!----------------------------------------------------------------------
!
!   routine name       - brefgeom3D
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - element boundary reference geometry routine
!
!   arguments :
!     in:
!             Mdle     - element middle node
!             Xi       - master element coordinates of a point
!             Etav     - reference vertex coordinates
!             ShapH    - values of shape functions at a point
!             GradH    - derivatives of the shape functions at the
!                        point
!             NrdofH   - number of element H1 vertex dof
!             Dxidt    - derivatives of a point on the element boundary
!                        wrt face parameters
!             Nsign    - sign factor indicating the actual direction of the
!                        normal wrt to the local one (implied by the
!                        face parametrization)
!     out:
!             Eta      - reference coordinates of the point
!             Detadxi  - derivatives wrt master coordinates
!             Dxideta  - inverse Jacobian
!             Detadt   - derivatives wrt the face parameters
!             Rn       - the normal unit vector
!             Bjac     - boundary Jacobian (norm of Detadxi)
!
!----------------------------------------------------------------------
!
   subroutine brefgeom3D(Mdle,Xi,Etav,ShapH,GradH,NrdofH,Dxidt,Nsign, &
                         Eta,Detadxi,Dxideta,Rjac,Detadt,Rn,Bjac)
!
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Mdle,NrdofH,Nsign
      real(8), intent(out) :: Rjac,Bjac
!
      real(8) :: Xi(3),Etav(3,8),                  &
                 ShapH(8),GradH(3,8),Dxidt(3,2),   &
                 Detadxi(3,3),Dxideta(3,3),        &
                 Eta(3),Detadt(3,2),Rn(3)
!
      integer :: i,j,iflag
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      call refgeom3D(Mdle,Xi,Etav,ShapH,GradH,NrdofH, &
                     Eta,Detadxi,Dxideta,Rjac,iflag)
      if (iflag.ne.0) then
        write(*,*) 'brefgeom3D: ERROR'; stop 1
      endif
!
!  ...derivatives of geometry map wrt face parameters
      Detadt(1:3,1:2)=0.d0
      do i=1,2
        do j=1,3
          Detadt(1:3,i) = Detadt(1:3,i) + Detadxi(1:3,j)*Dxidt(j,i)
        enddo
      enddo
!
!  ...surface normal
      call cross_product(Detadt(1:3,1),Detadt(1:3,2), Rn(1:3))
!
!  ...surface Jacobian
      call norm(Rn(1:3), Bjac)
!
!  ...unit normal vector
      Rn(1:3) = Rn(1:3)*Nsign/Bjac
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Eta,Rn,Rjac,Bjac
 7001   format('brefgeom3D: Eta,Rn,Rjac,Bjac = ',4(3f8.3,2x))
        call pause
      endif
#endif
!
   end subroutine brefgeom3D
