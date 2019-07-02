!*******************************************************************************************************************************************     
      SUBROUTINE spherical_triangle_OLD(nt,eta, x,Dx)     
!*******************************************************************************************************************************************
      !*** LATEST REVISION::    Dec 08
      !
      !*** PURPOSE:             routine parametrizes a spherical triagle using transfinite interpolation in the parameters
      !                         space (theta,phi)
      !
      !*** REQUIRED PROCEDURES: circular_segment, cylinderANDsphere, rotation, cartesian2spherical, mod3, Dspherical2cartesian
!*******************************************************************************************************************************************
      use GMP
!*******************************************************************************************************************************************
      IMPLICIT NONE
!*******************************************************************************************************************************************
      !*** EXTERNAL PROCEDURES
c      double precision, external :: lambda 
      integer, external :: mod3  
!-------------------------------------------------------------------------------------------------------------------------------------------
      !*** DUMMY ARGUMENTS
      integer, intent(in)                           :: nt                !*** triangle number
      double precision, dimension(2), intent(in)    :: eta               !*** master element coordinates
      double precision, dimension(3), intent(out)   :: x                 !*** coordinates in physical space (x,y,z)
      double precision, dimension(3,2), intent(out) :: Dx                !*** d(x,y,z)/d(eta1,eta2)
!--------------------------------------------------------------------------------------------------------------------------------------------
      !*** LOCAL VARIABLES
      integer                            :: ns                           !*** surface number
      integer                            :: nv                           !*** vertex number
      integer                            :: nc                           !*** curve number
      double precision, dimension(2)     :: r                            !*** coordinates in parameter space (theta,phi)
      double precision, dimension(2,2)   :: Dr                           !*** d(theta,phi)/d(eta1,eta2)
      double precision, dimension(3,3)   :: Dx_aux                       !*** auxiliary variable
      double precision, dimension(3,3)   :: RBM                          !*** rigid body motion
      double precision, dimension(3,3)   :: RBM_aux                      !*** auxiliary matrix
      double precision, dimension(3)     :: v_aux                        !*** auxiliary vector
      double precision, dimension(3)     :: center                       !*** sphere center
      double precision                   :: radius                       !*** sphere radius
      integer, dimension(3,2)            :: icurve                       !*** icurve(i,j): i = curve local number;
                                                                         !    j: | curve global number | orientation |
      double precision, dimension(3)     :: l                            !*** baricentrical coordinates
      double precision, dimension(3)     :: zeta                         !*** curvilinear coordinates for edges (one per edge)
      double precision, dimension(3)     :: zeta_aux                     !*** auxiliary variable
      double precision, dimension(3,3)   :: x_curve                      !*** curve parameterization in physical space
      double precision, dimension(3,3)   :: Dx_curve                     !*** derivative
      double precision, dimension(3,3)   :: V                            !*** vertices in physical space stored column-wise
      double precision, dimension(2,3)   :: A                            !*** vertices in parameter space stored column-wise   
      double precision, dimension(3,3)   :: A_aux                        !*** auxiliary matrix
      double precision, dimension(2,3)   :: r_curve                      !*** curve parameterization in parameter space
      double precision, dimension(2,3)   :: Dr_curve                     !*** d(phi,theta)/dzeta_j
      double precision, dimension(2,3)   :: r_ker                        !*** kernel functions in parameter space (one per edge)   
      double precision, dimension(2,3)   :: dr_KERdzeta                  !*** derivative of kernel functions
      double precision, dimension(2,2,3) :: dr_KERdeta                                
      double precision, dimension(3)     :: x_surface                    !*** explicite surface parameterization in cartesian coordinates
      double precision, dimension(3,3)   :: Dx_surface                   !*** derivative of explicit surface parameterization
      double precision                   :: rho, theta, phi              !*** spherical coordinates
      double precision, dimension(3,2)   :: LHS                          !*** left hand side of linear system; also stores the solution
      double precision, dimension(3)     :: RHS                          !*** right hand side of linear system; also stores the solution
      integer, dimension(4)              :: WORK                         !*** workspace array
      integer                            :: INFO                         !*** error flag for routine DGELS             
      integer                            :: i, j, k
      double precision, dimension(2,2)   :: temp1,temp2  
!-------------------------------------------------------------------------------------------------------------------------------------------
      !*** PARAMETERS
      double precision, parameter :: eps = 1.d-13                        !*** geometric tolerance
      double precision, dimension(3), parameter :: e1 = 
     .                                       (/1.d0, 0.d0, 0.d0/)        !*** vector (1,0,0)
      double precision, dimension(3), parameter :: e2 = 
     .                                       (/0.d0, 1.d0, 0.d0/)        !*** vector (0,1,0)
      double precision, dimension(3,2), parameter ::                     !*** dl_i/deta_j
     .  dldeta = reshape((/-1.d0, 1.d0, 0.d0, -1.d0, 0.d0, 1.d0/),
     .                                                        (/3,2/))
      double precision, dimension(3,2), parameter ::                     !*** dzeta_i/deta_j
     .  dzetadeta = reshape((/ 1.d0, -0.5d0, -0.5d0,
     .                               0.5d0, 0.5d0, -1.d0/), (/3,2/))
!*******************************************************************************************************************************************
!     STEP 0: get sphere data and check compatibility
!-------------------------------------------------------------------------------------------------------------------------------------------
      ns = TRIANGLES(nt)%Idata(1)                                        !*** get surface number
      !
      if (SURFACES(ns)%Type .ne. 'Sphere') then                          !*** check surface type
        write(*,*) 'spherical_triangle: inconsistent surface type.'
	    stop
      end if
      !
      center(1:3) = SURFACES(ns)%Rdata(1:3)                              !*** get sphere data
      radius = SURFACES(ns)%Rdata(4)    
      !
      icurve(1:3,2) = 1                                                  !*** initialize orientation to 1
      !
      do i = 1, 3                                                        
        nv = TRIANGLES(nt)%VertNo(i)                                     !*** get vertex number
        V(1:3,i) = POINTS(nv)%Rdata(1:3) - center(1:3)                   !*** store vertex coordinates wrt center
        icurve(i,1) = abs(TRIANGLES(nt)%EdgeNo(i))                       !*** get curve number
        if (nv .ne. CURVES(icurve(i,1))%EndPoNo(1)) icurve(i,2) = 0      !*** update orientation if needed
      end do                      
      !	  
!-------------------------------------------------------------------------------------------------------------------------------------------
!     STEP 1: distinguish between vertices and non-vertex
!-------------------------------------------------------------------------------------------------------------------------------------------
      if ((abs(eta(1)) .le. eps) .and.                                   !*** 1st VERTEX
     .                  (abs(eta(2)) .le. eps)) then                  
         nc = icurve(1,1)                                                !*** compute d(x,y,z)/deta1
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
!--------------------------------------------------------------------------------------------------------------------------------------------
      else if ((abs(eta(1) - 1.d0) .le. eps) .and.                       !*** 2nd VERTEX
     .                       (abs(eta(2)) .le. eps)) then       
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
!-------------------------------------------------------------------------------------------------------------------------------------------
      else if ((abs(eta(1)) .le. eps) .and.                              !*** 3rd VERTEX
     .                (abs(eta(2) - 1.d0) .le. eps)) then     
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
      else                                                               !*** NOT A VERTEX
      !
!-------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.1: get preprocessing data
!-------------------------------------------------------------------------------------------------------------------------------------------
        l(1) = 1.d0 - eta(1) - eta(2)                                    !*** get baricentrical coordinates... 
        l(2) = eta(1)
        l(3) = eta(2)
      !
        do i = 1, 3
          nc = icurve(i,1)
          zeta(i) = (l(mod3(i+1)) - l(mod3(i)) + 1.d0)/2.d0              !*** get curvilinear coordinate zeta
          if (icurve(i,2) .eq. 1) then                                   !*** account for orientation
            zeta_aux(i) = zeta(i)
            call curve(nc,zeta_aux(i), x_curve(1:3,i),Dx_curve(1:3,i))
          else
            zeta_aux(i) = 1.d0 - zeta(i)        
            call curve(nc,zeta_aux(i), x_curve(1:3,i),Dx_curve(1:3,i))
            Dx_curve(1:3,i) = -Dx_curve(1:3,i)
          end if
        end do
      !
!-------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.3: motion
!-------------------------------------------------------------------------------------------------------------------------------------------       
        call rotation(V(1:3,1),e1, RBM)                                  !*** determine rotation that maps 1st vertex on x-axis and...
      !  
        v_aux(1:3) = MATMUL(RBM(1:3,1:3), V(1:3,2))                      !*** ...find image of 2nd vertex under such rotation
        v_aux(1) = 0.d0                                                  !*** project onto yz-plane
        call rotation(v_aux,e2, RBM_aux)
      !  
        RBM(1:3,1:3) = MATMUL(RBM_aux(1:3,1:3), RBM(1:3,1:3))            !*** determine rigid body motion
      !
        V(1:3,1:3) = MATMUL(RBM(1:3,1:3), V(1:3,1:3))                    !*** rotate triangle in physical space
      !
!-------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.4: determine vertices and edges in parameter space
!-------------------------------------------------------------------------------------------------------------------------------------------
        do i = 1, 3
          call cartesian2spherical(V(1:3,i), A_aux(1,i),A_aux(2,i),      !*** determine spherical coordinates of vertices
     .                                                    A_aux(3,i))
        end do
        A(1:2,1:3) = A_aux(2:3,1:3)                                      !*** coordinates of vertices in parameter space (theta,phi)
      !
        do i = 1, 3                        
          x_curve(1:3,i)  = x_curve(1:3,i) - center(1:3) 
          x_curve(1:3,i)  = MATMUL(RBM(1:3,1:3), x_curve(1:3,i))         !*** apply RBM to x_curve 
          Dx_curve(1:3,i) = MATMUL(RBM(1:3,1:3), Dx_curve(1:3,i))        !*** apply RBM to Dx_curve
      !    
          call cartesian2spherical(x_curve(1:3,i), rho,theta,phi)        !*** get spherical coordinates
          r_curve(1,i) = theta                                           !*** store theta
          r_curve(2,i) = phi                                             !*** store phi
      !
!-------------------------------------------------------------------------------------------------------------------------------------------
!      STEP 1.5: determine derivatives of edges in parameter space
!-------------------------------------------------------------------------------------------------------------------------------------------
          call Dspherical2cartesian(radius,theta,phi, x_surface,         !*** compute x_surface and Dx_surface   
     .                                                      Dx_surface)  
      !
          LHS(1:3,1:2) = Dx_surface(1:3,2:3)
          RHS(1:3) = Dx_curve(1:3,i)   
      !
          call DGELS('N', 3, 2, 1, LHS, 3, RHS, 3, WORK, 4, INFO)        !*** compute least square solution
      !
          if (INFO .ne. 0) then
            write (*,*) 'spherical_triangle: linear system not solved.'
            stop
          end if
      ! 
          Dr_curve(1,i) = RHS(1)                                         !*** store d(theta,phi)/dzeta
          Dr_curve(2,i) = RHS(2)   
        end do
      !
!---------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.6: compute r_ker, Dr_ker
!---------------------------------------------------------------------------------------------------------------------------------------------
        do j = 1, 3                                                      !*** compute r_ker -- (one per edge)       
           r_ker(:,j) = (r_curve(:,j) - 
     .           (1.d0 - zeta(j))*A(:,j) - zeta(j)*A(:,mod3(j+1))) 
     .                        /(zeta(j)*(1.d0 - zeta(j)))
        end do
      !
        do j = 1, 3                                                      !*** dr_KER/dzeta -- (one per edge) 
           dr_KERdzeta(:,j) = ((Dr_curve(:,j) + A(:,j) 
     .              - A(:,mod3(j+1)))*zeta(j)*(1.d0 - zeta(j)) -
     .           (r_curve(:,j) - (1.d0 - zeta(j))*A(:,j) - zeta(j)*
     .            A(:,mod3(j+1)))*(1.d0 - 2.d0*zeta(j)))/
     .                            (zeta(j)*(1.d0 - zeta(j)))**2  
        end do
      !    
        do k = 1, 3                                                      !*** dr_KER/d(eta1,eta2) -- (one per edge) 
          do i = 1, 2
            do j = 1, 2
              dr_KERdeta(i,j,l) = dr_KERdzeta(i,k)*dzetadeta(k,j)
            end do    
          end do
        end do        
      !
!--------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1.7: compute r, Dr
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
          components_4: do i = 1, 2
            r(i) = r(i) + l(j)*l(mod3(j+1))*r_ker(i,j) 
	  end do components_4
        end do curves_4
      !
        Dr(1:2,1) = - A(1:2,1) + A(1:2,2)                                !*** derivative of linear contribution
        Dr(1:2,2) = - A(1:2,1) + A(1:2,3)
      !
c
c----------------------------------------------------------------------
c RECALL: let F be a vector and u be scalar, then Leibnitz rule
c         immediately applies: d(Fu) = dF u + F du where the last 
c         operation is a rank-1 update. If we consider the product 
c         uF instead of Fu, the dimensions of the jacobian matrices
c         do not conform any longer!!!
c----------------------------------------------------------------------
c
        do k = 1, 3                                                      !*** dr/deta  
           do i = 1, 2
             do j = 1, 2
               temp1(i,j) = r_ker(i,k)*dldeta(k,j)
               temp2(i,j) = r_ker(i,k)*dldeta(mod3(k+1),j)
             end do
           end do
           Dr = Dr + l(k)*l(mod3(k+1))*dr_KERdeta(:,:,k) +
     .          + l(mod3(k+1))*temp1(:,:) + l(k)*temp2(:,:)
        end do
      !
!-------------------------------------------------------------------------------------------------------------------------------------------
!       STEP 1,8: mapping between parameter space and physical space
!-------------------------------------------------------------------------------------------------------------------------------------------
        call Dspherical2cartesian(radius,r(1),r(2), x,Dx_aux)            !*** get (x,y,z) and d(x,y,z)/d(theta,phi)
        Dx(1:3,1:2) = Dx_aux(1:3,2:3)         
      !
        RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                           !*** invert rigid body motion
      !
        x(1:3) = MATMUL(RBM(1:3,1:3), x(1:3)) + center(1:3)              !*** apply inverse rigid body motion
      !
        Dx(1:3,1:2) = MATMUL(Dx(1:3,1:2), Dr(1:2,1:2))                   !*** d(x,y,z)/d(eta1,eta2) = d(x,y,z)/d(theta,phi) d(theta,phi)/d(eta1,eta2)
        Dx(1:3,1:2) = MATMUL(RBM(1:3,1:3), Dx(1:3,1:2))                  !*** apply inverse rigid body motion
      end if
      !
      END SUBROUTINE spherical_triangle_OLD
