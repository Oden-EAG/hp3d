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
c                        coordinates for the transfinite teterahedron
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
c
c-----------------------------------------------------------------------
c
      subroutine tetra_TraTet_OLD(No,Eta, X,Dxdeta)
c
c-----------------------------------------------------------------------
      use GMP
c-----------------------------------------------------------------------
      implicit none
      common /ctetra_TraTet/ iprint
      common /cface_bubble/ iprint_face_bubble
c-----------------------------------------------------------------------
c     EXTERNAL PROCEDURES
      double precision, external :: lambda3D
c-----------------------------------------------------------------------
c     DUMMY ARGUMENTS
      integer, intent(in)                           :: No                ! tetra number
      double precision, dimension(3), intent(in)    :: Eta               ! master tetrahedron coordinates
      double precision, dimension(3), intent(out)   :: X                 ! physical coordinates
      double precision, dimension(3,3), intent(out) :: Dxdeta            ! d(x1,x2,x3)/d(eta1,eta2,eta3)
c-----------------------------------------------------------------------
c     LOCAL VARIABLES
      double precision, dimension(3)   :: x_aux
      double precision, dimension(3,3) :: dx_aux
      double precision, dimension(3,4) :: xv(3,4)                        ! vertex coordinates
      integer                          :: np                             ! vertex number
      integer                          :: nc                             ! curve number
      integer                          :: nick                           ! triangle nickname (triangle number & orientation)
      integer                          :: iprint, i, j, ie, ivar,
     .                                    iprint_face_bubble
      double precision                 :: jac
c-----------------------------------------------------------------------
c
ccc  ...linear shape functions for a tetrahedron
cc      dimension norder(19),vshap(4),dvshap(3,4)
c
c      iprint = 0
c
      iprint_face_bubble = iprint
      if (iprint .eq. 1) then
        write(*,7001) No, Eta, TETRAS(No)%Type
 7001   format(1x,'tetra_TraTet -
     . TETRA: No, Eta, Type = ',i8,2x,3f8.3,2x,a10)
      endif
c
c  ...get the vertex coordinates
      do i = 1, 4
        np = TETRAS(No)%VertNo(i)
        call pointr(np, xv(1:3,i))
      enddo
c
c
c PG, Sep 08
ccc  ...temporary............................
cc      norder(1:11)=1
cc      call shape3t(Eta,norder, nrdof,vshap,dvshap)
cc
cc      X(1:3)=0.d0; Dxdeta(1:3,1:3) = 0.d0
cc      do k=1,4
cc        X(1:3) = X(1:3) + xv(1:3,k)*vshap(k)
cc        do i=1,3
cc          Dxdeta(1:3,i) = Dxdeta(1:3,i) + xv(1:3,k)*dvshap(i,k)
cc        enddo
cc      enddo
c
c
      X(1:3) = 0.d0                                                      ! compute linear interpolant
      do i = 1, 4
        X(1:3) = X(1:3) + xv(1:3,i)*lambda3D(i,Eta)
      end do
      Dxdeta(1:3,1) = xv(1:3,2) - xv(1:3,1)                              ! compute derivative of linear interpolant
      Dxdeta(1:3,2) = xv(1:3,3) - xv(1:3,1)
      Dxdeta(1:3,3) = xv(1:3,4) - xv(1:3,1)
c
      do i = 1, 6                                                        ! add edge bubbles
        nc = TETRAS(No)%EdgeNo(i)                                        !   get curve number
        call edge_bubble(nc,i,Eta, x_aux,dx_aux)                         !   compute edge bubble
c
        if (iprint .eq. 2) then
          if (i .eq. 1) then
           write(*,*)'tetra_TraTet: derivatives of bubble associated...'
            write(*,*)'...to edge 1 (dxdeta2,dxdeta3 should be 0)'
            write(*,*)'dxdeta1 = ',dx_aux(1:3,1)
            write(*,*)'dxdeta2 = ',dx_aux(1:3,2)
            write(*,*)'dxdeta3 = ',dx_aux(1:3,3)
            call pause
          end if
c
          if (i .eq. 3) then
           write(*,*)'tetra_TraTet: derivatives of bubble associated...'
            write(*,*)'...to edge 3 (dxdeta1,dxdeta3 should be 0)'
            write(*,*)'dxdeta2 = ',dx_aux(1:3,2)
            write(*,*)'dxdeta1 = ',dx_aux(1:3,1)
            write(*,*)'dxdeta3 = ',dx_aux(1:3,3)
            call pause
          end if
c
          if (i .eq. 4) then
           write(*,*)'tetra_TraTet: derivatives of bubble associated...'
            write(*,*)'...to edge 4 (dxdeta1,dxdeta2 should be 0)'
            write(*,*)'dxdeta3 = ',dx_aux(1:3,3)
            write(*,*)'dxdeta1 = ',dx_aux(1:3,1)
            write(*,*)'dxdeta2 = ',dx_aux(1:3,2)
            call pause
          end if
        end if
c
        X(1:3) = X(1:3) + x_aux(1:3)                                     !   accumulate coordinates
        Dxdeta(1:3,1:3) = Dxdeta(1:3,1:3) + dx_aux(1:3,1:3)              !   accumulate derivatives
        if (iprint .eq. 1) then
          write(*,*)'tetra_TraTet: added edge bubble'
          write(*,*)'X = ', X
          do ivar = 1, 3
            write(*,7002) ivar,Dxdeta(ivar,1:3)
 7002       format('               ivar, Dxdeta = ',i1,2x,3e12.5)
          enddo
          jac = Dxdeta(1,1)*Dxdeta(2,2)*Dxdeta(3,3) +
     .          Dxdeta(2,1)*Dxdeta(3,2)*Dxdeta(1,3) +
     .          Dxdeta(1,2)*Dxdeta(2,3)*Dxdeta(3,1) -
     .          Dxdeta(3,1)*Dxdeta(2,2)*Dxdeta(1,3) -
     .          Dxdeta(2,1)*Dxdeta(1,2)*Dxdeta(3,3) -
     .          Dxdeta(1,1)*Dxdeta(3,2)*Dxdeta(2,3)
          write(*,*) '             jac = ',jac
          call pause
        end if
      end do
c
      do i = 1, 4                                                        ! add face bubbles
        call face_bubble(No,i,Eta, x_aux,dx_aux)
        X(1:3) = X(1:3) + x_aux(1:3)                                     !   accumulate coordinates
        Dxdeta(1:3,1:3) = Dxdeta(1:3,1:3) + dx_aux(1:3,1:3)              !   accumulate derivatives
        if (iprint .eq. 1) then
          write(*,*)'tetra_TraTet: added face bubble'
          write(*,*)'X = ', X
          do ivar = 1, 3
            write(*,7002) ivar,Dxdeta(ivar,1:3)
          enddo
          jac = Dxdeta(1,1)*Dxdeta(2,2)*Dxdeta(3,3) +
     .          Dxdeta(2,1)*Dxdeta(3,2)*Dxdeta(1,3) +
     .          Dxdeta(1,2)*Dxdeta(2,3)*Dxdeta(3,1) -
     .          Dxdeta(3,1)*Dxdeta(2,2)*Dxdeta(1,3) -
     .          Dxdeta(2,1)*Dxdeta(1,2)*Dxdeta(3,3) -
     .          Dxdeta(1,1)*Dxdeta(3,2)*Dxdeta(2,3)
          write(*,*) '             jac = ',jac
          call pause
        end if
      end do
c
      end subroutine tetra_TraTet_OLD
c
c
c
c
c-----------------------------------------------------------------------
c
      subroutine edge_bubble(Nc,Ie,Eta, X,Dxdeta)
c
c-----------------------------------------------------------------------
c     LATEST REVISION: Dec 08
c
c     PURPOSE: routine evaluates edge bubble function for a tetrahedron
c
c     IN:      Nc - curve number
c              Ie - local number for edge
c              Eta - master tetrahedron coordinate
c
c     OUT:     X - physical coordinates
c              Dxdeta - derivative
c
c     REMARKS: reconstruction becomes undefined at verices of edge
c-----------------------------------------------------------------------
      use GMP
      use element_data
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c     EXTERNAL PROCEDURES
      double precision, external :: lambda3D
      integer, external          :: mod3
c-----------------------------------------------------------------------
c     DUMMY ARGUMENTS
      integer,                          intent(in)  :: Nc                ! curve number
      integer,                          intent(in)  :: Ie                ! local number for edge
      double precision, dimension(3),   intent(in)  :: Eta               ! master tetrahedron coordinate
      double precision, dimension(3),   intent(out) :: X                 ! physical coordinates
      double precision, dimension(3,3), intent(out) :: Dxdeta            ! derivatives
c-----------------------------------------------------------------------
c     LOCAL VARIABLES
      double precision                 :: zeta                           ! curvilinear coordinate
      double precision, dimension(3)   :: dzetadeta                      ! derivative
      double precision, dimension(2)   :: lambda                         ! affine coordinates wrt start and end point
      double precision, dimension(3,2) :: dl                             ! dlambda_i/deta_j
      double precision                 :: blend                          ! blending function
      double precision, dimension(3)   :: dblend                         ! derivatives of blending function
      integer                          :: nv1, nv2                       ! start and end point number
      double precision, dimension(3)   :: dxdzeta                        ! d(x1,x2,x3)/dzeta
      double precision, dimension(3,2) :: v                              ! vertices: i-th coordinate of j-th vertex
      integer                          :: i,j,iprint
      double precision                 :: temp
c-----------------------------------------------------------------------
c     PARAMETERS
      double precision, parameter, dimension(4,3) ::                     ! derivatives of 3D affine coordinates
     .  dlambda = reshape((/-1.d0, 1.d0, 0.d0, 0.d0,                     ! dlambda_i/deta_j
     .                      -1.d0, 0.d0, 1.d0, 0.d0,
     .                      -1.d0, 0.d0, 0.d0, 1.d0/),
     .                                             (/4,3/))
      double precision, parameter :: eps = 1.d-13                        ! geometric tolerance
c-----------------------------------------------------------------------
c
      iprint = 0
c
      if (iprint .eq. 1) then
        write(*,7002) Nc,Ie, Eta
 7002   format(1x,'edge_bubble: Nc, Ie, Eta = ',i8,2x,i3,2x,3f8.3)
      endif
c
      lambda(1) = lambda3D(TETRA_EDGE_TO_VERT(1,Ie), Eta)                ! edge baricentrical coordinates for point
      lambda(2) = lambda3D(TETRA_EDGE_TO_VERT(2,Ie), Eta)
      blend = lambda(1)*lambda(2)                                        ! blending function
c
      dl(1,1:3) = dlambda(TETRA_EDGE_TO_VERT(1,Ie),1:3)                  ! derivatives of edge baricentrical coordinates
      dl(2,1:3) = dlambda(TETRA_EDGE_TO_VERT(2,Ie),1:3)
      do j = 1, 3                                                        ! dblend/deta
        dblend(j) = dl(1,j)*lambda(2) + lambda(1)*dl(2,j)
      end do
c
      zeta = 0.5d0*(lambda(2) - lambda(1) + 1.d0)                        ! curvilinear coordinate wrt to LOCAL (tetra induced!) orientation for edge...
      dzetadeta(1:3) = 0.5d0*(dl(2,1:3) - dl(1,1:3))                     ! ... zeta = (l_end - l_start + 1)/2; dzeta/d(eta1,eta2,eta3)
c
      if (iprint .eq. 1) then
        write(*,*)'---------------------------------------------------'
        write(*,*)'edge_bubble...'
        write(*,*)'zeta = ',zeta
        write(*,*)'dzetadeta = ',dzetadeta(1:3)
        write(*,*)'---------------------------------------------------'
        call pause
      end if
c
      if (Nc .gt. 0) then                                                ! get start and end point (wrt LOCAL orientation!) of curve
        nv1 = CURVES(abs(Nc))%EndPoNo(1)
        nv2 = CURVES(abs(Nc))%EndPoNo(2)
      else
        nv1 = CURVES(abs(Nc))%EndPoNo(2)
        nv2 = CURVES(abs(Nc))%EndPoNo(1)
      end if
      v(1:3,1) = POINTS(nv1)%Rdata(1:3)
      v(1:3,2) = POINTS(nv2)%Rdata(1:3)
c
      if (abs(zeta) .le. eps) then                                       ! IF start point
        X(1:3) = 0.d0
        Dxdeta(1:3,1:3) = 0.d0
c
      else if (abs(zeta-1.d0) .le. eps) then                             ! ELSE IF end point
        X(1:3) = 0.d0
        Dxdeta(1:3,1:3) = 0.d0
c
      else                                                               ! ELSE neither start nor end point
        if (Nc .gt. 0) then                                              !   IF global and local orientations agree
          call curve(abs(Nc),zeta, X,dxdzeta)
        else                                                             !   ELSE global and local orientations do no agree
          call curve(abs(Nc),(1.d0-zeta), X,dxdzeta)
          dxdzeta(1:3) = -dxdzeta(1:3)
        end if                                                           !   END IF
c
        do i = 1, 3
          X(i) = X(i) - (1.d0 - zeta)*v(i,1) - zeta*v(i,2)               !   subtract linear interpolant
        end do
c
        do i = 1, 3                                                      !   compute dx/deta
          do j = 1, 3
            Dxdeta(i,j) = (dxdzeta(i) + v(i,1) - v(i,2))*dzetadeta(j)    !   rank-1 update
          end do
        end do
c
        if (iprint .eq. 1) then
          write(*,*)'-------------------------------------------------'
          write(*,*)'edge_bubble: subtraction linear interpolant      '
          write(*,*)'X = ', X(1:3)
          write(*,*)'Dxdeta1 = ', Dxdeta(1:3,1)
          write(*,*)'Dxdeta2 = ', Dxdeta(1:3,2)
          write(*,*)'Dxdeta3 = ', Dxdeta(1:3,3)
          write(*,*)'-------------------------------------------------'
          call pause
        end if
c
        do i = 1, 3                                                      !   compute derivatives of kernel function
          do j = 1, 3
            Dxdeta(i,j) = (Dxdeta(i,j)*zeta*(1.d0 - zeta) - X(i)*        !   rank-1 update
     .      (1.d0 - 2.d0*zeta)*dzetadeta(j))/(zeta**2*(1.d0 - zeta)**2)
          end do
        end do
c
        do i = 1, 3
          X(i) = X(i)/((1.d0 - zeta)*zeta)                               !   compute kernel function
        end do
c
        if (iprint .eq. 1) then
          write(*,*)'-------------------------------------------------'
          write(*,*)'edge_bubble: computing kernel function           '
          write(*,*)'X = ', X(1:3)
          write(*,*)'Dxdeta1 = ', Dxdeta(1:3,1)
          write(*,*)'Dxdeta2 = ', Dxdeta(1:3,2)
          write(*,*)'Dxdeta3 = ', Dxdeta(1:3,3)
          write(*,*)'-------------------------------------------------'
          call pause
        end if
c
        do i = 1, 3                                                      !   compute derivative of bubble function
          do j = 1, 3
            Dxdeta(i,j) = Dxdeta(i,j)*blend + X(i)*dblend(j)             !   rank-1 update
          end do
        end do
c
        do i = 1, 3
          X(i) = X(i)*blend                                              !   compute bubble function
        end do
c
        if (iprint .eq. 1) then
          write(*,*)'-------------------------------------------------'
          write(*,*)'edge_bubble: computing bubble function           '
          write(*,*)'X = ', X(1:3)
          write(*,*)'Dxdeta1 = ', Dxdeta(1:3,1)
          write(*,*)'Dxdeta2 = ', Dxdeta(1:3,2)
          write(*,*)'Dxdeta3 = ', Dxdeta(1:3,3)
          write(*,*)'-------------------------------------------------'
          call pause
        end if
c
      end if                                                             ! END IF
c
      end subroutine edge_bubble
c-----------------------------------------------------------------------
c
c
c
c
c
c
c
c-----------------------------------------------------------------------
c
      subroutine face_bubble(No,Iface,Eta, X,Dxdeta)
c
c-----------------------------------------------------------------------
c     LATEST REVISION: Dec 08
c
c     PURPOSE: routine evaluates face bubble function for a tetrahedron
c
c     IN:      No - tet number
c              Iface - face local number
c              Eta -  master tetrahedron coordinates
c
c     OUT:     X - physical coordinates
c              Dxdeta - derivatives
c
c     REMARKS: reconstruction becomes undefined at edges of face
c-----------------------------------------------------------------------
      use GMP
      use element_data
c-----------------------------------------------------------------------
      implicit none
      common /cface_bubble/ iprint
c-----------------------------------------------------------------------
c     EXTERNAL PROCEDURES
      double precision, external :: lambda3D
c-----------------------------------------------------------------------
c     DUMMY ARGUMENTS
      integer,                          intent(in)  :: No                ! triangle nickname (Nt,orientation)
      integer,                          intent(in)  :: Iface             ! local number for face
      double precision, dimension(3),   intent(in)  :: Eta               ! master triangle coordinates
      double precision, dimension(3),   intent(out) :: X                 ! physical coordinates
      double precision, dimension(3,3), intent(out) :: Dxdeta            ! derivatives
c-----------------------------------------------------------------------
c     LOCAL VARIABLES
      integer                          :: nt                             ! triangle number
      integer                          :: lab                            ! triangle orientation
      integer                          :: nc                             ! curve number
      double precision, dimension(3)   :: lambda2D                       ! 2D baricentrical coordinates
      double precision, dimension(4)   :: lambda                         ! 3D baricentrical coordinates
      double precision                 :: blend                          ! blending function
      double precision, dimension(3)   :: dblend                         ! dblend/deta
      integer                          :: nv                             ! vertex number
      double precision, dimension(3,3) :: v                              ! face vertices stored column-wise
      double precision, dimension(3,2) :: dxdeta_F
      double precision, dimension(3,2) :: dxdzeta                        ! d(x1,x2,x3)/d(zeta1,zeta2)
      double precision, dimension(2)   :: zeta                           ! face coordinate of projection (NOT accounting for orientation!)
      double precision, dimension(2)   :: eta_F                          ! face coordinate of projection (accounting for orientation!)
      double precision, dimension(3)   :: x_aux
      double precision, dimension(3)   :: temp
      double precision                 :: den
      double precision, dimension(3,3) :: dx_aux
      double precision, dimension(2,3) :: dzetadeta                      ! d(zeta1,zeta2)/d(eta1,eta2,eta3)
      double precision, dimension(3,3) :: dl                             ! d(lambda1,lambda2,lambda3)/d(eta1,eta2,eta3) derivatives of FACE coordinates
      integer                          :: i,j,k
      integer                          :: iprint
c-----------------------------------------------------------------------
c     PARAMETERS
      integer, parameter, dimension(4) ::                                ! tetra_face_to_opposite_vertex_connectivity
     .  face2opp_vert = (/4, 3, 1, 2/)
      double precision, parameter, dimension(4,3) ::                     ! derivatives of 3D baricentrical coordinates
     .  dlambda = reshape((/-1.d0, 1.d0, 0.d0, 0.d0,                     ! dlambda_i/deta_j
     .                      -1.d0, 0.d0, 1.d0, 0.d0,
     .                      -1.d0, 0.d0, 0.d0, 1.d0/),
     .                                             (/4,3/))
      double precision, parameter, dimension(2,3,4) ::
     .  dzeta = reshape((/      1.d0,       0.d0,                        ! FACE 1: dzeta_i/deta_j
     .                          0.d0,       1.d0,
     .                     1.d0/3.d0,  1.d0/3.d0,
     .
     .                          1.d0,       0.d0,                        ! FACE 2: dzeta_i/deta_j
     .                     1.d0/3.d0,  1.d0/3.d0,
     .                          0.d0,       1.d0,
     .
     .                    -1.d0/3.d0, -1.d0/3.d0,                        ! FACE 3: dzeta_i/deta_j
     .                     2.d0/3.d0, -1.d0/3.d0,
     .                    -1.d0/3.d0,  2.d0/3.d0,
     .
     .                     1.d0/3.d0,  1.d0/3.d0,                        ! FACE 4: dzeta_i/deta_j
     .                          1.d0,       0.d0,
     .                          0.d0,       1.d0/),
     .                                             (/2,3,4/))
      double precision, parameter, dimension(2,2,6) ::
     .  deta_Fdzeta = reshape((/1.d0,  0.d0,  0.d0,  1.d0,               ! ORIENT 0: deta_F/dzeta
     .                          0.d0, -1.d0,  1.d0, -1.d0,               ! ORIENT 1: deta_F/dzeta
     .                         -1.d0,  1.d0, -1.d0,  0.d0,               ! and so on...
     .                          0.d0,  1.d0,  1.d0,  0.d0,
     .                         -1.d0,  0.d0, -1.d0,  1.d0,
     .                          1.d0, -1.d0,  0.d0, -1.d0/),
     .                                                  (/2,2,6/))
      double precision, parameter :: eps = 1.d-10                        ! geometric tolerance
c-----------------------------------------------------------------------
c
cccc      iprint = 1
c
      if (iprint .eq. 1) then
        write(*,7001) No,Iface,Eta(1:3)
 7001   format('face_bubble: No,Iface,Eta = ',i5,2x,i1,2x,3e12.5)
      endif
c
      do i = 1, 3                                                        ! affine coordinates wrt LOCAL orientation of face
        lambda(i) = lambda3D(TETRA_FACE_TO_VERT(i,Iface), Eta)
      end do
      lambda(4) = lambda3D(face2opp_vert(Iface), Eta)                    ! affine coordinate wrt opposite vertex
      blend = lambda(1)*lambda(2)*lambda(3)                              ! compute blending function
c
      do i = 1, 3                                                        ! retrive dlambda_i/d(eta1,eta2,eta3)
        dl(i,1:3) = dlambda(TETRA_FACE_TO_VERT(i,Iface),1:3)
      end do
      do j = 1, 3                                                        ! compute derivative of blending function
        dblend(j) = dl(1,j)*lambda(2)*lambda(3) +
     .              lambda(1)*dl(2,j)*lambda(3) +
     .              lambda(1)*lambda(2)*dl(3,j)
      end do
c
      call decode(TETRAS(No)%FigNo(Iface), nt,lab)                       ! get triangle number and orientation
c
      do i = 1, 3                                                        ! get vertices
        j = TETRA_FACE_TO_VERT(i,Iface)
        nv = TETRAS(No)%VertNo(j)
        v(1:3,i) = POINTS(nv)%Rdata(1:3)
      end do
c
      lambda2D(1) = lambda(1) + lambda(4)/3.d0                           ! 2D affine coordinates of projection onto face
      lambda2D(2) = lambda(2) + lambda(4)/3.d0
      lambda2D(3) = lambda(3) + lambda(4)/3.d0
c
      zeta(1) = lambda2D(2)                                              ! local cartesian coordinates of projection
      zeta(2) = lambda2D(3)
c
      do i = 1, 2                                                        ! retrive d(zeta1,zeta2)/d(eta1,eta2,eta3)
        do j = 1, 3
          dzetadeta(i,j) = dzeta(i,j,Iface)
        end do
      end do
      if (iprint.eq.1) then
        do i=1,2
          write(*,7002) i,zeta(i),dzetadeta(i,1:3)
 7002     format('face_bubble: ivar,zeta(ivar),dzetadeta(ivar,*) = ',
     .                         i1,2x,e12.5,2x,3e12.5)
        enddo
      endif
c
      if (abs(zeta(2)) .lt. eps) then                                    ! IF first edge
        X(1:3) = 0.d0
        Dxdeta(1:3,1:3) = 0.d0
c
      else if (abs(1.d0-zeta(1)-zeta(2)) .lt. eps) then                  ! ELSEIF second edge
        X(1:3) = 0.d0
        Dxdeta(1:3,1:3) = 0.d0
c
      else if (abs(zeta(1)) .lt. eps) then                               ! ELSEIF third edge
        X(1:3) = 0.d0
        Dxdeta(1:3,1:3) = 0.d0
c
      else                                                               ! ELSE not an edge
        select case (lab)                                                !   account for face orientation
        case (0)
          eta_F(1) = lambda2D(2); eta_F(2) = lambda2D(3)
        case (1)
          eta_F(1) = lambda2D(3); eta_F(2) = lambda2D(1)
        case (2)
          eta_F(1) = lambda2D(1); eta_F(2) = lambda2D(2)
        case (3)
          eta_F(1) = lambda2D(3); eta_F(2) = lambda2D(2)
        case (4)
          eta_F(1) = lambda2D(1); eta_F(2) = lambda2D(3)
        case (5)
          eta_F(1) = lambda2D(2); eta_F(2) = lambda2D(1)
        case default
          write(*,*)'face_bubble: unknown face orientation.'
          stop
        end select
        if (iprint.eq.1) then
          write(*,7003) lab,eta_F(1:2)
 7003     format('face_bubble: orientation, eta_F = ',i1,2x,2e12.5)
        endif
c
        call trian(nt,eta_F, X,dxdeta_F)                                 !   compute coordinates in physical space and derivatives
        if (iprint.eq.1) then
          do i=1,3
            write(*,7004) i,X(i),Dxdeta(i,1:3)
 7004       format('face_bubble: ivar,X,Dxeta from trian = ',
     .                           i1,2x,e12.5,2x,3e12.5)
          enddo
        endif
c
        dxdzeta = MATMUL(dxdeta_F(1:3,1:2), deta_Fdzeta(1:2,1:2,lab+1))  !   chain rule
        Dxdeta  = MATMUL(dxdzeta(1:3,1:2),  dzetadeta(1:2,1:3))          !   chain rule
c
        do i = 1, 3
          X(i) = X(i) - (1.d0 - zeta(1) - zeta(2))*v(i,1) -              !   subtract linear interpolant
     .           zeta(1)*v(i,2) - zeta(2)*v(i,3)
        end do
c
        do i = 1, 3                                                      !   compute derivative
          do j = 1, 3
            Dxdeta(i,j) = Dxdeta(i,j) +                                  !   rank-1 update
     .                    v(i,1)*(dzetadeta(1,j) + dzetadeta(2,j)) -
     .                    v(i,2)*dzetadeta(1,j) -
     .                    v(i,3)*dzetadeta(2,j)
          end do
        end do
        if (iprint.eq.1) then
          do i=1,3
            write(*,7005) i,X(i),Dxdeta(i,1:3)
 7005       format('face_bubble: ivar,X,Dxeta after LI = ',
     .                           i1,2x,e12.5,2x,3e12.5)
          enddo
        endif
c
        do i = 1, 3                                                      !   subtract edge bubbles
          j = TETRA_FACE_TO_EDGE(i,Iface)
          nc = TETRAS(No)%EdgeNo(j)
          call edge_bubble(nc,TETRA_FACE_TO_EDGE(i,Iface),Eta,
     .                                                  x_aux,dx_aux)
          X(1:3) = X(1:3) - x_aux(1:3)
          Dxdeta(1:3,1:3) = Dxdeta(1:3,1:3) - dx_aux(1:3,1:3)
        end do
        if (iprint.eq.1) then
          do i=1,3
            write(*,7006) i,X(i),Dxdeta(i,1:3)
 7006       format('face_bubble: ivar,X,Dxeta of face bubble = ',
     .                           i1,2x,e12.5,2x,3e12.5)
          enddo
        endif

c
        den = (1 - zeta(1) - zeta(2))*zeta(1)*zeta(2)                    !   denominator of kernel function
        do j = 1, 3                                                      !   derivative of denominator of kernel function
          temp(j) = -(dzetadeta(1,j) + dzetadeta(2,j))*zeta(1)*zeta(2)+
     .              (1.d0 - zeta(1) - zeta(2))*dzetadeta(1,j)*zeta(2) +
     .              (1.d0 - zeta(1) - zeta(2))*zeta(1)*dzetadeta(2,j)
        end do
c
        do i = 1, 3                                                      !   derivative of kernel function
          do j = 1, 3
            Dxdeta(i,j) = (Dxdeta(i,j)*den - X(i)*temp(j))/(den**2)      !   rank-1 update
          end do
        end do
c
        do i = 1, 3                                                      !   compute kernel function
          X(i) = X(i)/den
        end do
        if (iprint.eq.1) then
          do i=1,3
            write(*,7007) i,X(i),Dxdeta(i,1:3)
 7007       format('face_bubble: ivar,X,Dxeta of face kernel = ',
     .                           i1,2x,e12.5,2x,3e12.5)
          enddo
        endif
c
        do i = 1, 3                                                      !   compute derivative of bubble function
          do j = 1, 3
            Dxdeta(i,j) = Dxdeta(i,j)*blend + X(i)*dblend(j)             !   rank-1 update
          end do
        end do
c
        do i = 1, 3                                                      !   compute bubble function
          X(i) =  X(i)*blend
        end do
c
      end if                                                             ! END IF
c
      end subroutine face_bubble
c-----------------------------------------------------------------------
c
c
c
c-----------------------------------------------------------------------
c
      double precision function lambda3D(I,Eta)
c
c-----------------------------------------------------------------------
c     LATEST REVISION: Sep 08
c
c     PURPOSE: function evaluates 3D baricentrical coordinates with
c              respect to master tetrahedron
c
c     IN:      I - baricentrical coordinate number
c                  1 -> (0,0,0)
c                  2 -> (1,0,0)
c                  3 -> (0,1,0)
c                  4 -> (0,0,1)
c
c              Eta - cartesian coordinates of point
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c     DUMMY ARGUMENTS
      integer, intent(in) :: I
      double precision, dimension(3), intent(in) :: Eta
c-----------------------------------------------------------------------
c
      select case (I)
      case (1)
        lambda3D = 1.d0 - Eta(1) - Eta(2) - Eta(3)
      case (2)
        lambda3D = Eta(1)
      case (3)
        lambda3D = Eta(2)
      case (4)
        lambda3D = Eta(3)
      case default
        lambda3D = 0.d0
        write(*,*)'lambda3D: unknown baricentrical coordinate'
        stop
      end select
c
      end
