!*************************************************************************************************************************
      SUBROUTINE coneANDplane(nc, ncone, nplane, eta, r, Dr)
!*************************************************************************************************************************    
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             routine parametrizes a curve that results from the intersection between a 
      !                         cone and a plane.
      !
      !*** REQUIRED PROCEDURES: rotation, cross_product, scalar_product, normalize, cartesian2cylindrical.
!*************************************************************************************************************************
      use GMP
!*************************************************************************************************************************
      IMPLICIT NONE
!*************************************************************************************************************************
      !*** DUMMY ARGUMENTS
      integer, intent(in) :: nc                                         !*** curve number
      integer, intent(in) :: ncone                                      !*** surface number for cone
      integer, intent(in) :: nplane                                     !*** surface number for plane
      double precision, intent(in) :: eta                               !*** curvilinear coordinate
      double precision, dimension(3), intent(out) :: r                  !*** cartesian coordinates
      double precision, dimension(3), intent(out) :: Dr                 !*** derivative wrt eta
      !
      !*** LOCAL VARIABLES
      double precision,  dimension(3) :: vertex, axis                   !*** vertex and axis for cylinder
      double precision :: alpha                                         !*** half-aperture for cylinder
      double precision, dimension(3) :: r_plane, n_plane                !*** point and normal for plane
      double precision, dimension(3) :: p1, p2                          !*** auxiliary variables
      integer :: nv                                                     !*** vertex number
      double precision, dimension(3,2) :: V                             !*** start and end point for curve
      double precision, dimension(3) :: v_aux                           !*** auxiliary variable
      double precision, dimension(3,3) :: RBM                           !*** rigid body motion
      double precision, dimension(3,3) :: A, B                          !*** auxiliary matrices
      double precision :: rho, theta, z                                 !*** cylindrical coordinates
      double precision :: theta_aux                                     !*** auxiliary variable
      double precision :: prod                                          !*** auxiliary variable
      double precision :: beta                                          !*** auxiliary variable
      double precision :: cosine                                        !*** auxiliary variable
      double precision :: pi_halves                                     !*** auxiliary variable
      integer :: i
      !
      !*** PARAMETERS
      double precision, dimension(3), parameter :: e1 = 
     .                               (/1.d0, 0.d0, 0.d0/)               !*** vector (1,0,0)
      double precision, dimension(3), parameter :: e3 = 
     .                               (/0.d0, 0.d0, 1.d0/)               !*** vector (0,0,1)
      double precision, parameter :: eps = 1.d-13                       !*** geometric tolerance
!*************************************************************************************************************************
      !  
      if (SURFACES(ncone)%Type .ne. 'Cone') then                        !*** check surface type
        write(*,*)'------------------------------------------'
        write(*,*)'coneANDplane:                             '
        write(*,*)'ERROR: inconsistent surface type for cone.'
        write(*,*)'------------------------------------------'
        stop
      end if
      ! 
      vertex = SURFACES(ncone)%Rdata(1:3)                               !*** get cone data
      axis = SURFACES(ncone)%Rdata(4:6)              
      alpha = SURFACES(ncone)%Rdata(7)
      !
      r_plane(1:3) = SURFACES(nplane)%Rdata(1:3)                        !*** get plane data
      if (SURFACES(nplane)%Type .eq. 'VecPt') then      
        n_plane(1:3) = SURFACES(nplane)%Rdata(4:6)
      else if (SURFACES(nplane)%Type .eq. 'ThreePt') then
        p1(1:3) = SURFACES(nplane)%Rdata(4:6) - r_plane(1:3)
        p2(1:3) = SURFACES(nplane)%Rdata(7:9) - r_plane(1:3)
        call cross_product(p1, p2, n_plane)
      else
        write(*,*)'-------------------------------------------'
        write(*,*)'coneANDplane:                              '
        write(*,*)'ERROR: inconsistent surface type for plane.'
        write(*,*)'-------------------------------------------'
        stop
      end if
      !
      do i = 1, 2                                                       !*** get vertices
        nv = CURVES(nc)%EndPoNo(i)
        V(1:3,i) = POINTS(nv)%Rdata(1:3) - vertex(1:3)
      end do
      !
      do i = 1, 2                                                       !*** flip axis orientation if needed
        call scalar_product(axis, V(1:3,i), prod)
        if (prod .lt. 0.d0) then                                        
          axis(1:3) = -axis(1:3)
          exit
        end if
      end do
      !
      v_aux(1:3) = vertex(1:3) - r_plane(1:3)
      call scalar_product(v_aux, n_plane, prod)
      !
      if (abs(prod) .lt. eps) then                                      !*** IF cone vertex lays on plane
        call normalize(axis)
        call normalize(n_plane)
        call scalar_product(axis, n_plane, cosine)
        beta = acos(cosine)
        pi_halves = acos(-1.d0)/2.d0
      !  
        if (abs(beta - pi_halves) .le. alpha) then                      !***** IF intersection does not reduce to vertex only
          do i = 1, 3                                                   !***** intersection is a straight line
            r(i) = (1.d0 - eta)*V(i,1) + eta*V(i,2) + vertex(i)
          end do
          Dr(1:3) = V(1:3,2) - V(1:3,1)
      !    
        else                                                            !***** ELSE intersection reduces to vertex only
          r(1:3) = vertex(1:3)
          Dr(1:3) = 0.d0
          write(*,*)'---------------------------------------------'
          write(*,*)'coneANDplane:                                '
          write(*,*)'WARNING: intersection reduces to origin only;'
          write(*,*)'WARNING: derivative is undefined at origin.  '
          write(*,*)'---------------------------------------------'
        end if                                                          !***** END IF
      !
      else                                                              !*** ELSE cone vertex does not lay on plane
        call rotation(axis, e3, A)                                      !*** determine rotation that aligns cone axis to z-axis      
        V(1:3,1:2) = MATMUL(A(1:3,1:3), V(1:3,1:2))                     !*** apply A to vertices
        v_aux(1:3) = V(1:3,1)                                           !*** project 1st vertex onto xy plane
        v_aux(3) = 0.d0
        call rotation(v_aux, e1, B)                                     !*** determine rotation that align projection to x-axis
        V(1:3,2) = MATMUL(B(1:3,1:3), V(1:3,2))                         !*** apply B to vertices
        call cartesian2cylindrical(V(1:3,2), rho, theta_aux, z)         !*** determine theta of end point
        theta = eta*theta_aux
      !
        RBM(1:3,1:3) = MATMUL(B(1:3,1:3), A(1:3,1:3))                   !*** determine rigid body motion
      !
        r_plane(1:3) = r_plane(1:3) - vertex(1:3)                       !*** apply rigid body motion to plane
        r_plane(1:3) = MATMUL(RBM(1:3,1:3), r_plane(1:3))
        n_plane(1:3) = MATMUL(RBM(1:3,1:3), n_plane(1:3))
      !
        call scalar_product(r_plane, n_plane, prod)                     !*** prod is needed for computing r(3)
      !
        r(3) = prod/(n_plane(1)*tan(alpha)*cos(theta) +                 !*** evaluate r
     .                n_plane(2)*tan(alpha)*sin(theta) + n_plane(3))
        r(1) = tan(alpha)*cos(theta)*r(3)
        r(2) = tan(alpha)*sin(theta)*r(3)
      !  
        Dr(3) = - (n_plane(2)*tan(alpha)*cos(theta) -                   !*** evaluate d(x,y,z)/dtheta
     .                      n_plane(1)*tan(alpha)*sin(theta))/
     .                              (n_plane(1)*tan(alpha)*cos(theta) +
     .                n_plane(2)*tan(alpha)*sin(theta) + n_plane(3))**2
        Dr(1) = tan(alpha)*(cos(theta)*Dr(3) - sin(theta)*r(3))
        Dr(2) = tan(alpha)*(sin(theta)*Dr(3) + cos(theta)*r(3))  
      !  
        Dr(1) = Dr(1)*theta_aux                                         !*** chain rule: d(x,y,z)/deta = d(x,y,z)/dtheta dtheta/deta
        Dr(2) = Dr(2)*theta_aux 
        Dr(3) = Dr(3)*theta_aux
      !
        RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                          !*** determine inverse rigid body motion
      !
        r(1:3) = MATMUL(RBM(1:3,1:3), r(1:3)) + vertex(1:3)             !*** apply inverse rigid body motion
        Dr(1:3) = MATMUL(RBM(1:3,1:3), Dr(1:3))
      end if                                                            !*** END IF
      !
      END SUBROUTINE coneANDplane
