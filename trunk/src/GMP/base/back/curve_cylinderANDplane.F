!************************************************************************************************************************
      SUBROUTINE cylinderANDplane(nc, ncylinder, nplane, eta, r, Dr)
!************************************************************************************************************************
      !*** LATEST REVISION:      Aug 08
      !
      !*** PURPOSE:              routine parametrizes a curve resulting from the intersection of a cylinder and a plane.
      !
      !*** REQUIRED PROCEDURES:  cross_product, rotation, scalar_product, cartesian2cylindrical.
!************************************************************************************************************************
      use GMP
!************************************************************************************************************************
      IMPLICIT NONE   
!************************************************************************************************************************
      !*** DUMMY ARGUMENTS   
      integer, intent(in) :: nc                                          !*** curve number
      integer, intent(in) :: ncylinder                                   !*** cylinder surface number
      integer, intent(in) :: nplane                                      !*** plane surface number
      double precision, intent(in) :: eta                                !*** curvilinear coordinate
      double precision, dimension(3), intent(out) :: r, Dr               !*** cartesian coordinates and derivatives
      !
      !*** LOCAL VARIABLES
      double precision, dimension(3) :: center, axis                     !*** center and axis for cylinder
      double precision :: radius                                         !*** cylinder radius
      double precision, dimension(3) :: r_plane, n_plane                 !*** data for plane of type 'VecPt'
      double precision, dimension(3) :: p1, p2                           !*** auxiliary variables
      double precision :: prod                                           !*** auxiliary variable
      double precision, dimension(3,3) :: RBM                            !*** rigid body motion                        
      double precision, dimension(3,3) :: RBM_aux                        !*** auxiliary matrix
      double precision, dimension(3,2) :: V                              !*** vertices stored columns-wise
      double precision, dimension(3) :: v_aux                            !*** auxiliary matrix
      double precision :: rho, theta, z                                  !*** cylindrical coordinates
      double precision :: theta_aux                                      !*** auxiliary variable
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
      center(1:3) = SURFACES(ncylinder)%Rdata(1:3)                       !*** get cylinder data
      axis(1:3) = SURFACES(ncylinder)%Rdata(4:6)
      radius = SURFACES(ncylinder)%Rdata(7)
      !
      r_plane(1:3) = SURFACES(nplane)%Rdata(1:3)                         !*** get plane data
      !
      if (SURFACES(nplane)%Type .eq. 'VecPt') then      
        n_plane(1:3) = SURFACES(nplane)%Rdata(4:6)
      else
        p1(1:3) = SURFACES(nplane)%Rdata(4:6) - r_plane(1:3)
        p2(1:3) = SURFACES(nplane)%Rdata(7:9) - r_plane(1:3)
        call cross_product(p1, p2, n_plane)
      end if
      !
      do i = 1, 2
        nv = CURVES(nc)%EndPoNo(i)
        V(1:3,i) = POINTS(nv)%Rdata(1:3) - center(1:3)
      end do
      !
      r_plane(1:3) = r_plane(1:3) - center(1:3)
      !
      call rotation(axis, e3, RBM)                                       !*** rotation that maps cylinder axis to z-axis
      !
      r_plane(1:3) = MATMUL(RBM(1:3,1:3), r_plane(1:3))                  !*** apply rotation to r_plane
      !
      call scalar_product(axis, n_plane, prod)
      !
      if (abs(prod) .le. eps) then                                       !*** IF plane is parallel to cylinder axis
        if (sqrt(r_plane(1)**2 + r_plane(2)**2) .le. radius) then        !***** IF intersection is not empty
      !
          do i = 1, 3                                                    !***** determine r
            r(i) = (1.d0 - eta)*V(i,1) + eta*V(i,2)                     
          end do
      !
          Dr(1:3) = V(1:3,2) - V(1:3,1)                                  !***** determine Dr
      !
          RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                         !***** invert rigid body motion
      !     
         r(1:3) = MATMUL(RBM(1:3,1:3), r(1:3)) + center(1:3)             !***** apply inverse rigid body motion
         Dr(1:3) = MATMUL(RBM(1:3,1:3), Dr(1:3))
      !
        else                                                             !***** ELSE intersection is empty
          write(*,*)'-----------------------------'
          write(*,*)'cylinderANDplane: '
          write(*,*)'ERROR: intersection is empty.'
          write(*,*)'-----------------------------'
          stop
      !
        end if                                                           !***** END IF
      !
      else                                                               !*** ELSE plane is not parallel to cylinder axis
        V(1:3,1:2) = MATMUL(RBM(1:3,1:3), V(1:3,1:2))                    !*** rotate vertices
      !         
        v_aux(1:3) = V(1:3,1)                                            !*** project 1st vertex onto xy-plane
        v_aux(3) = 0.d0
      !
        call rotation(v_aux, e1, RBM_aux)                                !*** rotation that maps 1st vertex to xz-plane
      !
        V(1:3,2) = MATMUL(RBM_aux(1:3,1:3), V(1:3,2))                    !*** determine theta
        call cartesian2cylindrical(V(1:3,2), rho, theta_aux, z)
        theta = eta*theta_aux
      !
        RBM(1:3,1:3) = MATMUL(RBM_aux(1:3,1:3), RBM(1:3,1:3))            !*** determine rigid body motion
      !
        r_plane(1:3) = MATMUL(RBM_aux(1:3,1:3), r_plane(1:3))            !*** apply rigid body motion to plane   
        n_plane(1:3) = MATMUL(RBM(1:3,1:3), n_plane(1:3)) 
      ! 
        call scalar_product(r_plane, n_plane, prod)                      !*** determine r
        r(1) = radius*cos(theta)
        r(2) = radius*sin(theta)
        r(3) = (prod - n_plane(1)*r(1) - n_plane(2)*r(2))/n_plane(3)
      !
        Dr(1) = -r(2)                                                    !*** d(x,y,z)/dtheta
        Dr(2) = r(1)
        Dr(3) = (n_plane(2)*r(1) - n_plane(1)*r(2))/n_plane(3)
      !
        Dr(1) = Dr(1)*theta_aux                                          !*** chain rule: d(x,y,z)/dtheta dtheta/deta
        Dr(2) = Dr(2)*theta_aux
        Dr(3) = Dr(3)*theta_aux
      !
        RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                           !*** invert rigid body motion
      !
        r(1:3) = MATMUL(RBM(1:3,1:3), r(1:3)) + center(1:3)              !*** apply inverse rigid body motion
        Dr(1:3) = MATMUL(RBM(1:3,1:3), Dr(1:3))
      !
      end if                                                             !*** END IF
      !
      END SUBROUTINE cylinderANDplane
