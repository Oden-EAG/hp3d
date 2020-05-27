!**************************************************************************************************************************
      SUBROUTINE cylinderANDsphere(nc, ncylinder, nsphere, eta, r, Dr)
!**************************************************************************************************************************
      !*** LATEST REVISION:       Aug 08
      !
      !*** PURPOSE:               routine parametrizes a curve resulting from the intersection of a cylinder and a sphere.
      !
      !*** REQUIRED SUBROUTINES:  planeANDline, rotation, cartesian2cylindrical, scalar_product.
!**************************************************************************************************************************
      use GMP
!**************************************************************************************************************************
      IMPLICIT NONE
!**************************************************************************************************************************
      !*** DUMMY ARGUMENTS
      integer, intent(in) :: nc                                          !*** curve number
      integer, intent(in) :: ncylinder                                   !*** cylinder number
      integer, intent(in) :: nsphere                                     !*** sphere number
      double precision, intent(in) :: eta                                !*** curvilinear coordinate
      double precision, dimension(3), intent(out) :: r                   !*** cartesian coordinates (x,y,z)
      double precision, dimension(3), intent(out) :: Dr                  !*** d(x,y,z)/deta
      !
      !*** LOCAL VARIABLES
      integer :: nv                                                      !*** vertex number
      double precision, dimension(3,2) :: V                              !*** curve vertices stored column-wise
      double precision, dimension(3) :: cylinder_center                  !*** cylinder center
      double precision, dimension(3) :: cylinder_axis                    !*** cylinder axis
      double precision, dimension(3) :: v_aux                            !*** auxiliary variable
      double precision :: cylinder_radius                                !*** cylinder radious
      double precision, dimension(3) :: sphere_center                    !*** sphere center
      double precision :: sphere_radius                                  !*** sphere radious
      double precision, dimension(3,3) :: RBM                            !*** rigid body motion
      double precision, dimension(3,3) :: RBM_aux                        !*** auxiliary matrix
      double precision :: rho, theta, z                                  !*** cylindrical coordinates
      double precision :: theta_aux                                      !*** auxiliary variable
      double precision :: Delta
      double precision :: prod
      integer :: i
      !
      !*** PARAMETERS
      double precision, dimension(3), parameter :: e1
     .                                = (/1.d0, 0.d0, 0.d0/)             !*** vector (1,0,0)
      double precision, dimension(3), parameter :: e3
     .                                = (/0.d0, 0.d0, 1.d0/)             !*** vector (0,0,1)
!*************************************************************************************************************************
      !
      do i = 1, 2                                                        !*** get vertices
        nv = CURVES(nc)%EndPoNo(i)
        V(1:3,i) = POINTS(nv)%Rdata(1:3)
      end do
      !
      cylinder_center(1:3) = SURFACES(ncylinder)%Rdata(1:3)              !*** get cylinder data
      cylinder_axis(1:3) = SURFACES(ncylinder)%Rdata(4:6)
      cylinder_radius = SURFACES(ncylinder)%Rdata(7)
      !
      sphere_center(1:3) = SURFACES(nsphere)%Rdata(1:3)                  !*** get sphere data
      sphere_radius = SURFACES(nsphere)%Rdata(4)
      !
      call planeANDline(cylinder_axis, sphere_center, cylinder_axis,     !*** redefine cylinder center
     .                             cylinder_center, v_aux)
      cylinder_center(1:3) = v_aux(1:3)
      !
      sphere_center(1:3) = sphere_center(1:3) - cylinder_center(1:3)     !*** redefine quantities wrt new cylinder center
      V(1:3,1) = V(1:3,1) - cylinder_center(1:3)
      V(1:3,2) = V(1:3,2) - cylinder_center(1:3)
      !
      call rotation(cylinder_axis, e3, RBM_aux)                          !*** get rotation that maps cylinder axis onto z-axis
      v_aux(1:3) = MATMUL(RBM_aux(1:3,1:3), V(1:3,1))                    !*** find image of 1st vertex
      v_aux(3) = 0.d0                                                    !*** project onto xy-plane
      call rotation(v_aux, e1, RBM)                                      !*** get rotation that maps projection onto x-axis
      RBM(1:3,1:3) = MATMUL(RBM(1:3,1:3), RBM_aux(1:3,1:3))              !*** determine rigid body motion
      !
      sphere_center(1:3) = MATMUL(RBM(1:3,1:3), sphere_center(1:3))      !*** apply rigid body motion
      V(1:3,1:2) = MATMUL(RBM(1:3,1:3), V(1:3,1:2))
      !
      call cartesian2cylindrical(V(1:3,2), rho, theta_aux, z)            !*** get theta coordinate for end point
      !
      theta = eta*theta_aux
      !
      Delta = sphere_radius**2 - cylinder_radius**2 -                    !*** evaluate Delta = z^2
     .   sphere_center(1)**2 - sphere_center(2)**2 +
     .   2.d0*sphere_center(1)*cylinder_radius*cos(theta) +
     .   2.d0*sphere_center(2)*cylinder_radius*sin(theta)
      !
      if (Delta .le. 0.d0) then
        write(*,*) 'cylinderANDsphere:
     . cylinder and sphere do not intersect.'
        stop
      end if
      !
      r(1) = cylinder_radius*cos(theta)                                  !*** evaluate r
      r(2) = cylinder_radius*sin(theta)
      r(3) = sqrt(Delta)
      !
      call scalar_product(e3, V(1:3,1), prod)
      if (prod .lt. 0.d0) r(3) = -r(3)                                   !*** invert sign of z-component if needed
      !
      Dr(1) = - r(2)                                                     !*** evaluate d(x,y,z)/dtheta; Dr(3) is well defined since r(3) in not zero
      Dr(2) = r(1)
      Dr(3) =  (r(1)*sphere_center(2) - r(2)*sphere_center(1))/r(3)
      !
      Dr(1) = Dr(1)*theta_aux                                            !*** chain rule: d(x,y,z)/deta = d(x,y,z)/dtheta dtheta/deta
      Dr(2) = Dr(2)*theta_aux
      Dr(3) = Dr(3)*theta_aux
      !
      RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                             !*** invert rigid body motion
      !
      r(1:3) = MATMUL(RBM(1:3,1:3), r(1:3)) + cylinder_center(1:3)       !*** apply inverse rigid body motion
      Dr(1:3) = MATMUL(RBM(1:3,1:3), Dr(1:3))
      !
      END SUBROUTINE cylinderANDsphere

