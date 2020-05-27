!************************************************************************************************************************
      SUBROUTINE sphereANDplane(nc, nsphere, nplane, eta, r, Dr)
!************************************************************************************************************************
      !*** LATEST REVISION:      Aug 08
      !
      !*** PURPOSE:              routine parametrizes a curve resulting from the intersection of a sphere and a plane.
      !
      !*** REQUIRED PROCEDURES:  cross_product, rotation, cartesian2cylindrical, Dcylindrical2cartesian.
!************************************************************************************************************************
      use GMP
!************************************************************************************************************************
      IMPLICIT NONE
!************************************************************************************************************************
      !*** DUMMY ARGUMENTS
      integer, intent(in) :: nc                                          !*** curve number
      integer, intent(in) :: nsphere                                     !*** cylinder surface number
      integer, intent(in) :: nplane                                      !*** plane surface number
      double precision, intent(in) :: eta                                !*** curvilinear coordinate
      double precision, dimension(3), intent(out) :: r, Dr               !*** cartesian coordinates and derivatives
      !
      !*** LOCAL VARIABLES
      double precision, dimension(3) :: center                           !*** sphere center
      double precision :: radius                                         !*** sphere radius
      double precision, dimension(3) :: r_plane, n_plane                 !*** data for plane of type 'VecPt'
      double precision, dimension(3) :: p1, p2                           !*** auxiliary variables
      double precision, dimension(3,3) :: RBM                            !*** rigid body motion
      double precision, dimension(3,3) :: RBM_aux                        !*** auxiliary matrix
      double precision, dimension(3) :: v_aux                            !*** auxiliary vector
      double precision, dimension(3,2) :: V                              !*** vertices stored columns-wise
      double precision :: rho, theta, z                                  !*** cylindrical coordinates
      double precision :: theta_aux                                      !*** auxiliary variable
      double precision, dimension(3,3) :: Dr_aux                         !*** auxiliary variable
      integer :: nv                                                      !*** vertex number
      integer :: i
      !
      !*** PARAMETERS
      double precision, dimension(3), parameter :: e1 =                  !*** vector (1,0,0)
     .                                         (/1.d0, 0.d0, 0.d0/)
      double precision, dimension(3), parameter :: e3 =                  !*** vector (0,0,1)
     .                                         (/0.d0, 0.d0, 1.d0/)
      double precision, parameter :: eps = 1.d-13                        !*** geometric tolerance
!************************************************************************************************************************
      !
      center(1:3) = SURFACES(nsphere)%Rdata(1:3)
      radius = SURFACES(nsphere)%Rdata(4)
      !
      r_plane(1:3) = SURFACES(nplane)%Rdata(1:3)
      !
      if (SURFACES(nplane)%Type .eq. 'VecPt') then
        n_plane(1:3) = SURFACES(nplane)%Rdata(4:6)
      else
        p1(1:3) = SURFACES(nplane)%Rdata(4:6) - r_plane(1:3)
        p2(1:3) = SURFACES(nplane)%Rdata(7:9) - r_plane(1:3)
        call cross_product(p1, p2, n_plane)
      end if
      !
      do i = 1, 2                                                        !*** get endpoints wrt center
        nv = CURVES(nc)%EndPoNo(i)
        V(1:3,i) = POINTS(nv)%Rdata(1:3) - center(1:3)
      end do
      !
      call rotation(n_plane, e3, RBM_aux)                                !*** rotation that maps cylinder axis to z-axis
      !
      v_aux(1:3) = MATMUL(RBM_aux(1:3,1:3), V(1:3,1))                    !*** determine images of 1st vertex
      !
      v_aux(3) = 0.d0                                                    !*** project onto xy-plane
      !
      call rotation(v_aux, e1, RBM)                                      !*** rotation that maps projection onto x-axis
      !
      RBM(1:3,1:3) = MATMUL(RBM(1:3,1:3), RBM_aux(1:3,1:3))              !*** determine rigid body motion
      !
      V(1:3,2) = MATMUL(RBM(1:3,1:3), V(1:3,2))                          !*** determine theta
      call cartesian2cylindrical(V(1:3,2), rho, theta_aux, z)
      theta = eta*theta_aux
      !
      if (z .gt. radius) then
        write(*,*)'sphereANDplane: sphere and plane do not intersect.'
        stop
      end if
      !
      V(1:3,2) = MATMUL(RBM(1:3,1:3), V(1:3,2))                           !*** determine theta
      call cartesian2cylindrical(V(1:3,2), rho, theta_aux, z)
      theta = eta*theta_aux
      !
      call Dcylindrical2cartesian(rho, theta, z, r, Dr_aux)               !*** get (x,y,z) and d(x,y,z)/dtheta
      Dr(1:3) = Dr_aux(1:3,2)
      !
      Dr(1) = Dr(1)*theta_aux                                             !*** chain rule: d(x,y,z)/deta = d(x,y,z)/dtheta dtheta/deta
      Dr(2) = Dr(2)*theta_aux
      Dr(3) = Dr(3)*theta_aux
      !
      RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                              !*** invert rigid body motion
      !
      r(1:3) = MATMUL(RBM(1:3,1:3), r(1:3)) + center(1:3)
      Dr(1:3) = MATMUL(RBM(1:3,1:3), Dr(1:3)) + center(1:3)
      !
      END SUBROUTINE sphereANDplane
