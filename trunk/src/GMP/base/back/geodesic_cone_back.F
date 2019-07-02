!**************************************************************************************************************
      SUBROUTINE cone_geodesic_old(nc, eta, r, Dr)
!**************************************************************************************************************
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             routine parametrizes a cone geodesic.
      !
      !*** REQUIRED PROCEDURES: scalar_product, rotation, cone2polar, Dpolar2cone, Dcartesian2polar.
      !
!**************************************************************************************************************      
      use GMP
!**************************************************************************************************************      
      IMPLICIT NONE
!**************************************************************************************************************
      !*** DUMMY ARGUMENTS      
      integer, intent(in) :: nc                                          !*** curve number
      double precision, intent(in) :: eta                                !*** curvilinear coordinate
      double precision, dimension(3), intent(out) :: r                   !*** cartesian coordinates (x,y,z)
      double precision, dimension(3), intent(out) :: Dr                  !*** d(x,y,z)/deta
      !
      !*** LOCAL VARIABLES
      integer :: ns                                                      !*** surface number
      double precision, dimension(3) :: vertex                           !*** cone vertex
      double precision, dimension(3) :: axis                             !*** cone axis
      double precision :: alpha                                          !*** cone half-aperture
      integer :: nv                                                      !*** end point number
      double precision, dimension(3,2) :: V                              !*** curve end points
      double precision :: prod                                           !*** auxiliary variable
      double precision :: pi                                             !*** pi
      double precision, dimension(3) :: v_aux                            !*** auxiliary vector
      double precision, dimension(3,3) :: RBM                            !*** rigid body motion
      double precision, dimension(3,3) :: RBM_aux                        !*** auxiliary matrix
      double precision :: rho, theta                                     !*** polar coordinates     
      double precision, dimension(2,2) :: A                              !*** end points in parameter space
      double precision, dimension(2) :: r_par                            !*** coordinates (x,y) in parameter space
      double precision, dimension(2) :: Dr_par                           !*** derivatives d(x,y)/d(eta) in parameter space 
      double precision, dimension(3,2) :: Dr_aux                         !*** auxiliary variables
      double precision, dimension(2) :: polar                            !*** coordinates (rho,theta) in parameter space
      double precision, dimension(2,2) :: Dpolar                         !*** derivatives d(rho,theta)/d(eta1,eta2) in parameter space
      integer :: i
      !
      !*** PARAMETERS
      double precision, dimension(3), parameter :: e1 
     .                                = (/1.d0, 0.d0, 0.d0/)             !*** vector (1,0,0)
      double precision, dimension(3), parameter :: e3 
     .                                = (/0.d0, 0.d0, 1.d0/)             !*** vector (0,0,1)
!**************************************************************************************************************  
      !   
      ns = CURVES(nc)%Idata(1)
      !
      if (SURFACES(ns)%Type .ne. 'Cone') then
        write(*,*)'cone_geodesic: inconsistent surface type.'
        stop
      end if
      !
      vertex(1:3) = SURFACES(ns)%Rdata(1:3)                              !*** get cone data
      axis(1:3) = SURFACES(ns)%Rdata(4:6)
      alpha = SURFACES(ns)%Rdata(7)
      !                                                                   
      do i = 1, 2                                                        !*** get endpoints wrt center
        nv = CURVES(nc)%EndPoNo(i)
        V(1:3,i) = POINTS(nv)%Rdata(1:3) - vertex(1:3)
      end do
      !
      do i = 1, 2
        call scalar_product(axis, V(1:3,i), prod)                        !*** account for orientation of axis
        if (prod .lt. 0.d0) then   
          pi = acos(-1.d0)                                     
          alpha = pi - alpha
          exit
        end if
      end do
      !
      call rotation(axis, e3, RBM_aux)
      !
      v_aux = MATMUL(RBM_aux, V(1:3,1))                                  !*** get image of start point
      v_aux(3) = 0.d0                                                    !*** project onto xy-plane 
      call rotation(v_aux, e1, RBM)
      !
      RBM(1:3,1:3) = MATMUL(RBM(1:3,1:3), RBM_aux(1:3,1:3))              !*** rigid body motion
      !
      V(1:3,1:2) = MATMUL(RBM(1:3,1:3), V(1:3,1:2))                      !*** apply rigid body motion to end points
      !
      do i = 1, 2                                                        !*** get end points in parameter space
        call cone2polar(alpha, V(1:3,i), rho, theta)
        A(1,i) = rho*cos(theta)
        A(2,i) = rho*sin(theta)
      end do
      !
      do i = 1, 2                                                        !*** get cartesian coordinates in parameter space
        r_par(i) = (1.d0 - eta)*A(i,1) + eta*A(i,2)
      end do
      !
      Dr_par = A(1:2,2) - A(1:2,1)                                       !*** d(x,y)/deta
      !
      call Dcartesian2polar(r_par, polar, Dpolar)                        !*** get (rho,theta) and d(rho,theta)/d(x,y)
      !
      Dr_par(1:2) = MATMUL(Dpolar(1:2,1:2), Dr_par(1:2))                 !*** chain rule --> d(rho,theta)/deta = d(rho,theta)/d(x,y) d(x,y)/deta
      !
      call Dpolar2cone(alpha, polar(1), polar(2), r, Dr_aux)             !*** get (x,y,z) and d(x,y,z)/d(rho,theta)
      !
      Dr(1:3) = MATMUL(Dr_aux(1:3,1:2), Dr_par(1:2))                     !*** chain rule --> d(x,y,z)/deta = d(x,y,z)/d(rho,theta) d(rho,theta)/deta
      !
      RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                             !*** invert rigid body motion
      !
      r(1:3) = MATMUL(RBM(1:3,1:3), r(1:3)) + vertex(1:3)                !*** apply inverse rigid body motion
      Dr(1:3) = MATMUL(RBM(1:3,1:3), Dr(1:3))
      !  
      END SUBROUTINE cone_geodesic
