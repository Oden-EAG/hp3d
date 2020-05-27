!**************************************************************************************************************************
      SUBROUTINE cylindrical_triangle_OLD(nt, eta, x, Dx)
!**************************************************************************************************************************
      !*** LEATEST REVISION:    Aug 08
      !
      !*** PURPOSE:             routine parametrizes a cylindrical triangle using transfinite interpolatin in the
      !                         parameters space (theta,z)
      !
      !*** REQUIRED PROCEDURES: lambda, mod3, cylinderANDsphere, cylinderANDplane,  rotation, catesian2cylidrical,
      !                         Dcartesian2cylindrical, Dcylindrical2cartesian
!**************************************************************************************************************************
      use GMP
!**************************************************************************************************************************
      IMPLICIT NONE
!**************************************************************************************************************************
      !*** EXTERNAL PROCEDURES
      double precision, external :: lambda                               !*** baricentrical coordinates
      integer, external :: mod3
      !
      !*** DUMMY ARGUMENTS
      integer, intent(in) :: nt                                          !*** triangle number
      double precision, dimension(2), intent(in) :: eta                  !*** master element coordinates
      double precision, dimension(3), intent(out) :: x                   !*** coordinates in physical space
      double precision, dimension(3,2), intent(out) :: Dx                !*** d(x,y,z)/d(eta1,eta2)
      !
      !*** LOCAL VARIABLES
      integer :: ns                                                      !*** cylinder surface number
      integer :: nc                                                      !*** curve number
      double precision, dimension(3) :: center                           !*** cylinder center
      double precision, dimension(3) :: axis                             !*** cylinder axis
      double precision :: radius                                         !*** cylinder radius
      double precision, dimension(3) :: l                                !*** baricentrical coordinates
      double precision, dimension(2,3) :: Dl                             !*** derivatives
      integer :: nv                                                      !*** vertex number
      double precision, dimension(2) :: r                                !*** coordinates in parameter space (theta,z)
      double precision, dimension(2,2) :: Dr                             !*** d(theta,z)/d(eta1,eta2)
      double precision, dimension(3,3) :: V                              !*** vertices in physical space, stored column-wise
      double precision, dimension(2,3) :: A                              !*** vertices in parameter space, stored column-wise
      double precision, dimension(3,3) :: A_aux                          !*** auxiliary matrix
      double precision, dimension(3,3) :: Dx_aux                         !*** auxiliary matrix
      integer, dimension(3,3) :: icurve                                  !*** curve(i,j): i = local number;
                                                                         !    j: | global number | orientation | curve type |
      double precision, dimension(3) :: zeta                             !*** curvilinear coordinates in parameter space
      double precision, dimension(3) :: zeta_aux                         !*** auxiliary variable
      double precision, dimension(2,3) :: Dzeta                          !*** derivatives
      double precision, dimension(3,3) :: x_curve                        !*** curve parameterization in physical space
      double precision, dimension(3,3) :: Dx_curve                       !*** derivatives
      double precision, dimension(2,3) :: r_curve                        !*** curve parameterization in parameter space
      double precision, dimension(2,3) :: Dr_curve                       !*** derivatives
      double precision, dimension(2,3) :: r_kernel                       !*** kernel functions in parameter space
      double precision, dimension(2,3) :: Dr_kernel                      !*** derivatives
      double precision, dimension(3,3) :: RBM                            !*** rigid body motion
      double precision, dimension(3,3) :: RBM_aux                        !*** auxiliary matrix
      double precision, dimension(3) :: v_aux                            !*** auxiliary vector
      double precision, dimension(3) :: x_surface                        !*** surface parameterization
      double precision, dimension(3,3) :: Dx_surface                     !*** derivatives
      double precision :: rho, theta, z                                  !*** cylindrical coordinates
      double precision, dimension(3,2) :: LHS                            !*** left hand side
      double precision, dimension(3) :: RHS                              !*** right hand side; also stores the solution
      integer, dimension(4) :: WORK                                      !*** workspace array
      integer :: INFO                                                    !*** error flag for routine DGELS
      integer :: i, j, k
      !
      !*** PARAMETERS
      double precision, dimension(3), parameter :: e1 =
     .                                     (/1.d0, 0.d0, 0.d0/)          !*** vector (1,0,0)
      double precision, dimension(3), parameter :: e3 =
     .                                     (/0.d0, 0.d0, 1.d0/)          !*** vector (0,0,1)
      double precision, parameter :: eps = 1.d-13                        !*** geometric tolerance
!**************************************************************************************************************************
!     STEP 0: get cylinder data and check compatibility
!--------------------------------------------------------------------------------------------------------------------------
      ns = TRIANGLES(nt)%Idata(1)                                        !*** get surface number
      !
      if (SURFACES(ns)%Type .ne. 'Cylinder') then                        !*** check surface type
        write(*,*) 'cylindrical_triangle:
     . inconsistent surface type'
	stop
      end if
      !
      center(1:3) = SURFACES(ns)%Rdata(1:3)                              !*** get surface data
      axis(1:3) = SURFACES(ns)%Rdata(4:6)
      radius = SURFACES(ns)%Rdata(7)
      !
      icurve(1:3,2) = 1                                                  !*** initialize orientation to 1
      icurve(1:3,3) = 1                                                  !*** initialize curve type to straight segment
      !
      do i = 1, 3
        nv = TRIANGLES(nt)%VertNo(i)                                     !*** get vertex number
        V(1:3,i) = POINTS(nv)%Rdata(1:3) - center(1:3)                   !*** store vertex coordinates wrt center
        icurve(i,1) = abs(TRIANGLES(nt)%EdgeNo(i))                       !*** get curve number
        if (nv .ne. CURVES(icurve(i,1))%EndPoNo(1)) icurve(i,2) = 0      !*** update orientation if needed
      end do
      !
!---------------------------------------------------------------------------------------------------------------------------
!     STEP 1: distinguish between vertices and non-vertex
!---------------------------------------------------------------------------------------------------------------------------
      if ((abs(eta(1)) .le. eps) .and.
     .                          (abs(eta(2)) .le. eps)) then             !*** 1st VERTEX
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
        x(1:3) = V(1:3,1) + center(1:3)                                  !*** compute (x,y,z)
!-----------------------------------------------------------------------------------------------------------------------------
      else if ((abs(eta(1) - 1.d0) .le. eps) .and.
     .                          (abs(eta(2)) .le. eps)) then             !*** 2nd VERTEX
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
        x(1:3) = V(1:3,2) + center(1:3)                                  !*** compute (x,y,z)
!------------------------------------------------------------------------------------------------------------------------------
      else if ((abs(eta(1)) .le. eps) .and.
     .                          (abs(eta(2) - 1.d0) .le. eps)) then      !*** 3rd VERTEX
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
        x(1:3) = V(1:3,3) + center(1:3)                                  !*** compute (x,y,z)
      !
      else                                                               !*** NOT a vertex
      !
!---------------------------------------------------------------------------------------------------------------------------
!       STEP 1.1: get preprocessing data
!---------------------------------------------------------------------------------------------------------------------------
        do i = 1, 3                                                      !*** get baricentrical coordinates...
          l(i) = lambda(i, eta)
        end do
      !
        Dl(1:2,1) = (/-1.d0, -1.d0/)                                     !*** ...and their derivatives
        Dl(1:2,2) = (/ 1.d0,  0.d0/)
        Dl(1:2,3) = (/ 0.d0,  1.d0/)
      !
        Dzeta(1:2,1) = (/ 1.0d0,  0.5d0/)                                !*** get derivatives of curvilinear coordinates
        Dzeta(1:2,2) = (/-0.5d0,  0.5d0/)
        Dzeta(1:2,3) = (/-0.5d0, -1.0d0/)
     !
!----------------------------------------------------------------------------------------------------------------------------
!       STEP 1.2: reconstruct triangle's edges in physical space
!----------------------------------------------------------------------------------------------------------------------------
        do i = 1, 3
          nc = icurve(i,1)
          if (CURVES(nc)%Type .ne. '1SurfsCur') then                     !*** IF curve is not a straight segment
            icurve(i,3) = 0
            zeta(i) = (l(mod3(i+1)) - l(mod3(i)) + 1.d0)/2.d0            !*** get curvilinear coordinate zeta
            if (icurve(i,2) .eq. 1) then
              zeta_aux(i) = zeta(i)
              call curve(nc, zeta_aux(i), x_curve(1:3,i),
     .                                                Dx_curve(1:3,i))
            else
              zeta_aux(i) = 1.d0 - zeta(i)
              call curve(nc, zeta_aux(i), x_curve(1:3,i),
     .                                                Dx_curve(1:3,i))
              Dx_curve(1:3,i) = -Dx_curve(1:3,i)
            end if
          end if
        end do
      !
!-------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.3: motion
!-------------------------------------------------------------------------------------------------------------------------------------
        call rotation(axis, e3, RBM)                                     !*** determine rotation that maps cylinder axis on z-axis
        v_aux(1:3) = MATMUL(RBM(1:3,1:3), V(1:3,1))                      !*** determine image of 1st vector
        v_aux(3) = 0.d0                                                  !*** project onto xy-plane
      !
        call rotation(v_aux, e1, RBM_aux)
      !
        RBM(1:3,1:3) = MATMUL(RBM_aux(1:3,1:3), RBM(1:3,1:3))            !*** determine rigid body motion
        V(1:3,1:3) = MATMUL(RBM(1:3,1:3), V(1:3,1:3))                    !*** rotate triangle in physical space
      !
!--------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.4: determine vertices and edges in pameter space
!--------------------------------------------------------------------------------------------------------------------------------------
        do i = 1, 3                                                      !*** determine vertices in parmeter space
          call cartesian2cylindrical(V(1:3,i), A_aux(1,i), A_aux(2,i),
     .                                         A_aux(3,i))
        end do
      !
        A(1:2,1:3) = A_aux(2:3,1:3)                                      !*** store vertices
      !
        do i = 1, 3
          if (icurve(i,3) .ne. 1) then                                   !*** IF curve is not a cylinder geodesic
            x_curve(1:3,i) = x_curve(1:3,i) - center(1:3)
            x_curve(1:3,i) = MATMUL(RBM(1:3,1:3), x_curve(1:3,i))        !*** apply RBM to r_curve
            Dx_curve(1:3,i) = MATMUL(RBM(1:3,1:3), Dx_curve(1:3,i))      !*** apply RBM to Dr_curve
      !
            call cartesian2cylindrical(x_curve(1:3,i), rho, theta, z)    !*** get coordinates in parameter space
            r_curve(1,i) = theta                                         !*** store theta
            r_curve(2,i) = z                                             !*** store z
      !
            call Dcylindrical2cartesian(radius, theta, z, x_surface,
     .                                                      Dx_surface)  !*** compute r_surface and Dr_surface
      !
            LHS(1:3,1:2) = Dx_surface(1:3,2:3)
            RHS(1:3) = Dx_curve(1:3,i)
            call DGELS('N', 3, 2, 1, LHS, 3, RHS, 3, WORK, 4, INFO)      !*** compute least square solution
      !
            if (INFO .ne. 0) then                                        !*** check error flag
              write (*,*) 'cylindrical_triangle:
     .  linear system not solved.'
              stop
            end if
      !
            Dr_curve(1,i) = RHS(1)                                       !*** store derivatives
            Dr_curve(2,i) = RHS(2)
          end if                                                         !*** END IF
        end do
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.5: compute r_kernel, Dr_kernel
!---------------------------------------------------------------------------------------------------------------------------------------------
        curves_1: do j = 1, 3
          if (icurve(j,3) .eq. 1) cycle                                  !*** IF curve is a straight segment cycle
          components_1: do i = 1, 2
            r_kernel(i,j) = (r_curve(i,j) -
     .	          (1.d0 - zeta(j))*A(i,j) - zeta(j)*A(i,mod3(j+1)))
     .                        /(zeta(j)*(1.d0 - zeta(j)))
          end do components_1
        end do curves_1
      !
        curves_2: do j = 1, 3                                            !*** derivative of kernel function wrt zeta
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
!       STEP 1.6: compute r, Dr
!--------------------------------------------------------------------------------------------------------------------------------------------
        r(1:2) = 0.d0
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
!---------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.7: mapping between parameter space and physical space
!---------------------------------------------------------------------------------------------------------------------------------------
        call Dcylindrical2cartesian(radius, r(1), r(2), x, Dx_aux)       !*** get (x,y,z) and d(x,y,z)/d(theta,z)
        Dx(1:3,1:2) = Dx_aux(1:3,2:3)
        !
        RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                           !*** invert rigid body motion
        !
        x(1:3) = MATMUL(RBM(1:3,1:3), x(1:3)) + center(1:3)              !*** apply inverse rigid body motion
        !
        Dx(1:3,1:2) = MATMUL(Dx(1:3,1:2), Dr(1:2,1:2))                   !*** d(x,y,z)/d(eta1,eta2) = d(x,y,z)/d(theta,z) d(theta,z)/d(eta1,eta2)
        Dx(1:3,1:2) = MATMUL(RBM(1:3,1:3), Dx(1:3,1:2))                  !*** apply inverse rigid body motion
      end if
      !
      END SUBROUTINE cylindrical_triangle_OLD
