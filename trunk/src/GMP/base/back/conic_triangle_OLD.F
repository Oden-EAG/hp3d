!*********************************************************************************************************************************************
      SUBROUTINE conic_triangle_OLD(nt, eta, x, Dx)
!*********************************************************************************************************************************************
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             routine parametrizes a conic triangle using transfinite interpolation in the parameters space.
      !
      !*** REQUIRED PROCEDURES: scalar_product, coneANDplane, rotation, cone2polar, Dpolar2cone, Dcartesian2cone.
      !
      !*** REMARKS:             when we unfold a cone onto a plane we obtain a circular sector in the (x,y) coordinates system; in such
      !                         system the cone geodesics are straight lines. When mapping such circular sector back into a cone, however,
      !                         we want to use polar coordinates (rho,theta) in place of (x,y).
!*********************************************************************************************************************************************
      use GMP
!*********************************************************************************************************************************************
      IMPLICIT NONE
!*********************************************************************************************************************************************
      !*** EXTERNAL PROCEDURES
      double precision, external :: lambda
      integer, external :: mod3
      !
      !*** DUMMY ARGUMENTS
      integer, intent(in) :: nt                                          !*** triangle number
      double precision, dimension(2), intent(in) :: eta                  !*** master element coordinates
      double precision, dimension(3) :: x                                !*** coordinates in physical space (x,y,z)
      double precision, dimension(3,2) :: Dx                             !*** d(x,y,z)/d(eta1,eta2)
      !
      !*** LOCAL VARIABLES
      double precision, dimension(2) :: r                                !*** cartesian coordinates (x,y) in parameter space
      double precision, dimension(2,2) :: Dr                             !*** d(x,y)/d(eta1,eta2)
      double precision, dimension(2) :: r_aux                            !*** polar coordinates (rho,theta) in parameter space
      double precision, dimension(2,2) :: Dr_aux                         !*** d(rho,theta)/d(eta1,eta2)
      double precision, dimension(3,3) :: RBM                            !*** rigid body motion
      double precision, dimension(3,3) :: RBM_aux                        !*** auxiliary matrix
      double precision, dimension(3) :: v_aux                            !*** auxiliary vector
      double precision, dimension(3) :: vertex                           !*** cone vertex
      double precision, dimension(3) :: axis                             !*** cone axis
      double precision :: alpha                                          !*** cone half-aperture
      integer :: ns                                                      !*** surface number
      integer :: nc                                                      !*** curve number
      integer, dimension(3,3) :: icurve                                  !*** icurve(i,j): i = local number;
                                                                         !    j: | global number | orientation | curve type |
      double precision, dimension(3) :: zeta                             !*** curvilinear coordinates in parameter space (one per edge)
      double precision, dimension(3) :: zeta_aux                         !*** auxiliary variables
      double precision, dimension(2,3) :: Dzeta                          !*** derivatives of curvilinear coordinates (one per edge)
      double precision, dimension(3) :: l                                !*** baricentrical coordinate (one per edge)
      double precision, dimension(2,3) :: Dl                             !*** derivatives of baricentrical coordinates (one per edge)
      integer :: nv                                                      !*** vertex number
      double precision, dimension(3,3) :: V                              !*** vertices in physical space, stored column-wise
      double precision, dimension(2,3) :: A                              !*** vertices in parameter space, stored column-wise
      double precision, dimension(3,3) :: x_curve                        !*** curve parameterization in physical space (one per edge)
      double precision, dimension(3,3) :: Dx_curve                       !*** derivative
      double precision, dimension(2,3) :: r_curve                        !*** curve parameterization in parameter space (one per edge)
      double precision, dimension(2,3) :: Dr_curve                       !*** derivative
      double precision, dimension(2,3) :: r_kernel                       !*** kernel functions in parameter space (one per edge)
      double precision, dimension(2,3) :: Dr_kernel                      !*** derivative
      double precision, dimension(3) :: x_surface                        !*** explicite surface parameterization in cartesian coordinates
      double precision, dimension(3,2) :: Dx_surface                     !*** derivative of explicit surface parameterization
      double precision :: rho, theta                                     !*** polar coordinates
      double precision, dimension(3,2) :: LHS                            !*** left hand side of linear system; also stores the solution
      double precision, dimension(3) :: RHS                              !*** right hand side of linear system; also stores the solution
      integer, dimension(4) :: WORK                                      !*** workspace array
      integer :: INFO                                                    !*** error flag for routine DGELS
      double precision :: pi                                             !*** pi
      double precision :: prod                                           !*** auxiliary variable
      double precision, dimension(2,2) :: M                              !*** auxiliary matrix
      integer :: i, j, k
      !
      !*** PARAMETERS
      double precision, parameter :: eps = 1.d-13                        !*** geometric tolerance
      double precision, dimension(3), parameter :: e1 =
     .                               (/1.d0, 0.d0, 0.d0/)                !*** vector (1,0,0)
      double precision, dimension(3), parameter :: e3 =
     .                               (/0.d0, 0.d0, 1.d0/)                !*** vector (0,0,1)
c
!*********************************************************************************************************************************************
!     STEP 0: get cone data and check compatibility
!---------------------------------------------------------------------------------------------------------------------------------------------
      ns = TRIANGLES(nt)%Idata(1)                                        !*** get surface number
      !
      if (SURFACES(ns)%Type .ne. 'Cone') then                            !*** check surface type
        write(*,*) 'conic_triangle: INCONSISTENT SURFACE TYPE'
        stop 1
      end if
      !
      vertex(1:3) = SURFACES(ns)%Rdata(1:3)                              !*** get cone data
      axis(1:3) = SURFACES(ns)%Rdata(4:6)
      alpha = SURFACES(ns)%Rdata(7)
      !
      icurve(1:3,2) = 1                                                  !*** initialize orientations to 1
      icurve(1:3,3) = 1                                                  !*** initialize curve type to straight segment
      !
      do i = 1, 3
        nv = TRIANGLES(nt)%VertNo(i)                                     !*** get vertex number
        V(1:3,i) = POINTS(nv)%Rdata(1:3) - vertex(1:3)                   !*** store vertex coordinates wrt vertex
        icurve(i,1) = abs(TRIANGLES(nt)%EdgeNo(i))                       !*** get curve number
        if (nv .ne. CURVES(icurve(i,1))%EndPoNo(1)) icurve(i,2) = 0      !*** update orientation if needed
      end do
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!     STEP 1: distinguish between vertices and non-vertex
!---------------------------------------------------------------------------------------------------------------------------------------------
      if ((abs(eta(1)) .le. eps) .and. (abs(eta(2)) .le. eps)) then      !*** 1st VERTEX
        nc = icurve(1,1)                                                 !*** compute d(x,y,z)/deta1
        if (icurve(1,2) .eq. 1) then
          call curve(nc, 0.d0, x, Dx(1:3,1))
        else
          call curve(nc, 1.d0, x, Dx(1:3,1))
          Dx(1:3,1) = -Dx(1:3,1)
        end if
      !
        nc = icurve(3,1)                                                 !*** compute d(x,y,z)/deta2
        if (icurve(3,2) .eq. 1) then
          call curve(nc, 1.d0, x, Dx(1:3,2))
          Dx(1:3,2) = -Dx(1:3,2)
        else
          call curve(nc, 0.d0, x, Dx(1:3,2))
        end if
      !
        x(1:3) = V(1:3,1) + vertex(1:3)                                  !*** compute (x,y,z)
!---------------------------------------------------------------------------------------------------------------------------------------------
      else if ((abs(eta(1) - 1.d0) .le. eps) .and.                       !*** 2nd VERTEX
     .                                (abs(eta(2)) .le. eps)) then
        nc = icurve(1,1)                                                 !*** compute d(x,y,z)/deta1
        if (icurve(1,2) .eq. 1) then
          call curve(nc, 1.d0, x, Dx(1:3,1))
        else
          call curve(nc, 0.d0, x, Dx(1:3,1))
          Dx(1:3,1) = -Dx(1:3,1)
        end if
      !
        nc = icurve(2,1)                                                 !*** compute d(x,y,z)/deta2
        if (icurve(2,2) .eq. 1) then
          call curve(nc, 0.d0, x, Dx(1:3,2))
        else
          call curve(nc, 1.d0, x, Dx(1:3,2))
          Dx(1:3,2) = -Dx(1:3,2)
        end if
        do i = 1, 3
          Dx(i,2) = Dx(i,1) + sqrt(2.d0)*Dx(i,2)
        end do
      !
        x(1:3) = V(1:3,2) + vertex(1:3)                                  !*** compute (x,y,z)
!---------------------------------------------------------------------------------------------------------------------------------------------
      else if ((abs(eta(1)) .le. eps) .and.                              !*** 3rd VERTEX
     .                          (abs(eta(2) - 1.d0) .le. eps)) then
        nc = icurve(3,1)                                                 !*** compute d(x,y,z)/deta2
        if (icurve(3,2) .eq. 1) then
          call curve(nc, 0.d0, x, Dx(1:3,2))
          Dx(1:3,2) = -Dx(1:3,2)
        else
          call curve(nc, 1.d0, x, Dx(1:3,2))
        end if
      !
        nc = icurve(2,1)                                                 !*** compute d(x,y,z)/deta1
        if (icurve(2,2) .eq. 1) then
          call curve(nc, 1.d0, x, Dx(1:3,1))
        else
          call curve(nc, 0.d0, x, Dx(1:3,1))
          Dx(1:3,1) = -Dx(1:3,1)
        end if
        do i = 1, 3
          Dx(i,1) = -sqrt(2.d0)*Dx(i,1) - Dx(i,2)
        end do
      !
        x(1:3) = V(1:3,3) + vertex(1:3)                                  !*** compute (x,y,z)
      !
      else                                                               !*** NOT A VERTEX
!---------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.1: get l, Dl, Dzeta
!---------------------------------------------------------------------------------------------------------------------------------------------
        do i = 1, 3                                                      !*** get baricentrical coordinates...
          l(i) = lambda(i, eta)
        end do
      !
        Dl(1:2,1) = (/-1.d0, -1.d0/)                                     !*** ...and their derivatives
        Dl(1:2,2) = (/ 1.d0,  0.d0/)
        Dl(1:2,3) = (/ 0.d0,  1.d0/)
      !
        Dzeta(1:2,1) = (/ 1.0d0,  0.5d0/)                                !*** derivatives of curvilinear coordinate zeta
        Dzeta(1:2,2) = (/-0.5d0,  0.5d0/)
        Dzeta(1:2,3) = (/-0.5d0, -1.0d0/)
      !
        do i = 1, 3                                                      !*** account for orientation of axis
          call scalar_product(axis, V(1:3,i), prod)
          if (prod .lt. 0.d0) then
            pi = acos(-1.d0)
            alpha = pi - alpha
            exit
          end if
        end do
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.2: reconstruct triangle's edges in physical space
!---------------------------------------------------------------------------------------------------------------------------------------------
        do i = 1, 3
          nc = icurve(i,1)
!          if (CURVES(nc)%Type .ne. '1SursCur') then                      !*** IF curve is not a cylinder geodesic
            icurve(i,3) = 0
            zeta(i) = (l(mod3(i+1)) - l(mod3(i)) + 1.d0)/2.d0            !*** get curvilinear coordinate zeta
            if (icurve(i,2) .eq. 1) then
              zeta_aux(i) = zeta(i)
              call curve(nc, zeta_aux(i), x_curve(1:3,i),
     .                                                Dx_curve(1:3,i))
            else
              zeta_aux(i) = 1.d0 - zeta(i)                               !*** account for orientation
              call curve(nc, zeta_aux(i), x_curve(1:3,i),
     .                                                Dx_curve(1:3,i))
              Dx_curve(1:3,i) = -Dx_curve(1:3,i)
            end if
!          end if
        end do
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.3: motion
!---------------------------------------------------------------------------------------------------------------------------------------------
        call rotation(axis, e3, RBM)                                     !*** align cone axis with z-axis
      !
        v_aux(1:3) =  MATMUL(RBM(1:3,1:3), V(1:3,1))                     !*** apply rotation to 1st vertex
        v_aux(3) = 0.d0                                                  !*** project onto xy-plane
        call rotation(v_aux, e1, RBM_aux)                                !*** rotation that brings 1st vertex onto xz-plane
      !
        RBM(1:3,1:3) = MATMUL(RBM_aux(1:3,1:3), RBM(1:3,1:3))            !*** determine rigid body motion
      !
        V(1:3,1:3) = MATMUL(RBM(1:3,1:3), V(1:3,1:3))                    !*** apply rotation to vertices
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.4: determine vertices and edges in parameter space
!---------------------------------------------------------------------------------------------------------------------------------------------
        do i = 1, 3                                                      !*** get vertices in parameter space (x,y)
          call cone2polar(alpha, V(1:3,i), rho, theta)
          A(1,i) = rho*cos(theta)
          A(2,i) = rho*sin(theta)
        end do
      !
        do i = 1, 3
!          if (icurve(i,3) .ne. 1) then                                   !*** IF curve is not a cone geodesic
            x_curve(1:3,i) = x_curve(1:3,i) - vertex(1:3)
            x_curve(1:3,i) = MATMUL(RBM(1:3,1:3), x_curve(1:3,i))        !*** apply RBM to x_curve
            Dx_curve(1:3,i) = MATMUL(RBM(1:3,1:3), Dx_curve(1:3,i))      !*** apply RBM to Dx_curve
      !
            call cone2polar(alpha, x_curve(1:3,i), rho, theta)           !*** get edge in parmeter space (x(zeta),y(zeta))
            r_curve(1,i) = rho*cos(theta)
            r_curve(2,i) = rho*sin(theta)
            M(1,1) = cos(theta); M(1,2) = -rho*sin(theta)                !*** compute d(x,y)/d(rho,theta)
            M(2,1) = sin(theta); M(2,2) = rho*cos(theta)
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!      STEP 1.5: determine derivatives of edges in parameter space
!---------------------------------------------------------------------------------------------------------------------------------------------
            call Dpolar2cone(alpha, rho, theta, x_surface, Dx_surface)   !*** compute x_surface and Dx_surface
      !
            LHS(1:3,1:2) = Dx_surface(1:3,1:2)
            RHS(1:3) = Dx_curve(1:3,i)
      !
            call DGELS('N', 3, 2, 1, LHS, 3, RHS, 3, WORK, 4, INFO)      !*** compute d(rho,theta)/dzeta using least squares
      !
            if (INFO .ne. 0) then
              write (*,*) 'conic_triangle: linear system not solved.'
              stop
            end if
      !
            Dr_curve(1:2,i) = MATMUL(M(1:2,1:2), RHS(1:2))               !*** chain rule: d(x,y)/dzeta = d(x,y)/d(rho,theta) d(rho,theta)/dzeta
!          end if                                                         !*** END IF
        end do
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.6: compute r_kernel, Dr_kernel
!---------------------------------------------------------------------------------------------------------------------------------------------
        curves_1: do j = 1, 3                                            !*** compute r_kernel = (x(zeta),y(zeta))
          if (icurve(j,3) .eq. 1) cycle                                  !*** IF curve is a straight segment cycle
          components_1: do i = 1, 2
            r_kernel(i,j) = (r_curve(i,j) -
     .	          (1.d0 - zeta(j))*A(i,j) - zeta(j)*A(i,mod3(j+1)))
     .                        /(zeta(j)*(1.d0 - zeta(j)))
          end do components_1
        end do curves_1
      !
        curves_2: do j = 1, 3                                            !*** compute Dr_kernel = d(x,y)/dzeta
          if (icurve(j,3) .eq. 1) cycle                                  !*** IF curve is a straight segment cycle
          components_2: do i = 1, 3
            Dr_kernel(i,j) = ((Dr_curve(i,j) + A(i,j) - A(i,mod3(j+1)))*
     .                                     zeta(j)*(1.d0 - zeta(j)) -
     .(r_curve(i,j) - (1.d0 - zeta(j))*A(i,j) - zeta(j)*A(i,mod3(j+1)))*
     .                                     (1.d0 - 2.d0*zeta(j)))/
     .         (zeta(j)*(1.d0 - zeta(j)))**2
          end do components_2
        end do curves_2
      !
!--------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.7: compute r, Dr
!--------------------------------------------------------------------------------------------------------------------------------------------
        r(1:2) = 0.d0                                                    !*** compute r = (x(eta1,eta2),y(eta1,eta2))
      !
        curves_3: do j = 1, 3                                            !*** linear contribution
	  components_3: do i = 1, 2
	    r(i) = r(i) + l(j)*A(i,j)
          end do components_3
        end do curves_3
      !
        curves_4: do j = 1, 3                                            !*** kernel contribution
          if (icurve(j,3) .eq. 1) cycle
          components_4: do i = 1, 2
            r(i) = r(i) + l(mod3(j))*l(mod3(j+1))*r_kernel(i,j)
	  end do components_4
        end do curves_4
      !
      !                                                                  !*** compute Dr = d(x,y)/d(eta1,eta2)
        Dr(1:2,1) = - A(1:2,1) + A(1:2,2)                                !*** derivative of linear contribution
        Dr(1:2,2) = - A(1:2,1) + A(1:2,3)
      !
        curves_5: do j = 1, 3                                            !*** derivative of kernel contribution
          if (icurve(j,3) .eq. 1) cycle
          components_5: do i = 1, 2
            derivatives_5: do k = 1, 2
              Dr(i,k) = Dr(i,k) +
     .   (Dl(k,mod3(j))*l(mod3(j+1)) + l(mod3(j))*Dl(k,mod3(j+1)))*
     .                                                 r_kernel(i,j) +
     .   l(mod3(j))*l(mod3(j+1))*Dr_kernel(i,j)*Dzeta(k,j)
            end do derivatives_5
          end do components_5
        end do curves_5
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.7: mapping between parameter space and physical space
!---------------------------------------------------------------------------------------------------------------------------------------------
        call Dcartesian2polar(r, r_aux, Dr_aux)                          !*** get (rho(eta1,eta2),theta(eta1,eta2)) and d(rho,theta)/d(x,y)
        call Dpolar2cone(alpha, r_aux(1), r_aux(2), x, Dx)               !*** get (x,y,z) and d(x,y,z)/d(rho,theta)
      !
        Dr(1:2,1:2) = MATMUL(Dr_aux(1:2,1:2), Dr(1:2,1:2))               !*** d(rho,theta)/d(eta1,eta2) = d(rho,theta)/d(x,y) d(x,y)/d(eta1,eta2)
      !
        RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                           !*** invert rigid body motion
      !
        x(1:3) = MATMUL(RBM(1:3,1:3), x(1:3)) + vertex(1:3)              !*** apply inverse rigid body motion
        Dx(1:3,1:2) = MATMUL(Dx(1:3,1:2), Dr(1:2,1:2))                   !*** d(x,y,z)/d(eta1,eta2) = d(x,y,z)/d(rho,theta) d(rho,theta)/d(eta1,eta2)
        Dx(1:3,1:2) = MATMUL(RBM(1:3,1:3), Dx(1:3,1:2))                  !*** apply inverse rigid body motion
      end if
      !
      END SUBROUTINE conic_triangle_OLD
