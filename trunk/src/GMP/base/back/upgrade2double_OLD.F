!****************************************************************************************************************************
!*** WARNING: this routine is specifically tailored for the geometry of the HEAD PROBLEM. There is NO guarantee
!***          whatsoever it will work for ANY other geometry!
!****************************************************************************************************************************
      SUBROUTINE upgrade2double_OLD
!****************************************************************************************************************************
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             routines updates single precision data provided by NETGEN to double precision data.
      !
      !*** REQUIRED PROCEDURES: rotation, planeANDline, cartesian2cylindrical, cylindrical2cartesian
!****************************************************************************************************************************
      use GMP
!****************************************************************************************************************************
      IMPLICIT NONE
!****************************************************************************************************************************
      !*** LOCAL VARIABLES
      double precision :: eps                                               !*** geometric tolerance
      integer :: npoint                                                     !*** vertex number
      integer :: nt                                                         !*** triangle number
      integer :: ns                                                         !*** surface number
      integer :: nv                                                         !*** vertex number
      double precision, dimension(3) :: v                                   !*** triangle vertex
      double precision, dimension(3) :: v_aux                               !*** auxiliary variable
      double precision, dimension(3,3) :: RBM                               !*** rigid body motion
      double precision, dimension(3,3) :: RBM_aux                           !*** auxiliary matrix
      integer :: status                                                     !*** integer flag
      logical :: B_status                                                   !*** boolean flag
      integer :: i
      integer :: iprint
      !
      !*** PARAMETERS
      double precision, dimension(3), parameter :: e1 =
     .                                              (/1.d0, 0.d0, 0.d0/)    !*** unit vector (1, 0, 0)
      double precision, dimension(3), parameter :: e3 =
     .                                              (/0.d0, 0.d0, 1.d0/)    !*** unit vector (0, 0, 1)
!****************************************************************************************************************************
      !
      write(*,*)'upgrade2double: upgrading to double
     . precision geometry...'
      !
      iprint = 0
      !
      !*** STEP 1: determine geometric tolerance for initial mesh. Geometric tolerance is an upper bound for the distance
      !***         between points and all of the algebraic surfaces they belong to.
      !
      !
      if (iprint .eq. 1) then
        call geometric_tolerance(eps, npoint)
        write(*,*) '-------------------------------------------------'
        write(*,*) 'upgrade2double:                                  '
        write(*,*) 'initial tolerance, npoint = ', eps, npoint
        write(*,*) '-------------------------------------------------'
        call print_GMP
      end if
      !
      !----------------------------------------------------------------------------------------------------------------------
      !*** STEP 2: project onto algebraic surfaces one at a time; the order of surfaces is - hopefully! - cleverly
      !***         picked according to the geometry of the problem.
      !
      call projectONTOplanes
      !
      if (iprint .eq. 1) then
        call geometric_tolerance(eps, npoint)
        write(*,*) '-------------------------------------------------'
        write(*,*) 'upgrade2double:                                  '
        write(*,*) 'tolerance after planes, npoint = ', eps, npoint
        write(*,*) '-------------------------------------------------'
        call print_GMP
      end if
      !
      !-----------------------------------------------------------------------------------------------------------------------
      call projectONTOcones
      !
      if (iprint .eq. 1) then
        call geometric_tolerance(eps, npoint)
        write(*,*) '-------------------------------------------------'
        write(*,*) 'upgrade2double:                                  '
        write(*,*) 'tolerance after cones, npoint = ', eps, npoint
        write(*,*) '-------------------------------------------------'
        call print_GMP
      end if
      !
      !-----------------------------------------------------------------------------------------------------------------------
      call projectONTOspheres
      !
      if (iprint .eq. 1) then
        call geometric_tolerance(eps, npoint)
        write(*,*) '-------------------------------------------------'
        write(*,*) 'upgrade2double:                                  '
        write(*,*) 'tolerance after spheres, npoint = ', eps, npoint
        write(*,*) '-------------------------------------------------'
        call print_GMP
      end if
      !
      !-----------------------------------------------------------------------------------------------------------------------
      call projectONTOcylinders
      !
      if (iprint .eq. 1) then
        call geometric_tolerance(eps, npoint)
        write(*,*) '-------------------------------------------------'
        write(*,*) 'upgrade2double:                                  '
        write(*,*) 'tolerance after cylinders, npoint = ', eps, npoint
        write(*,*) '-------------------------------------------------'
        call print_GMP
      end if
      !
      CONTAINS
      !
      !
      !**********************************************************************************************************************
      SUBROUTINE geometric_tolerance(eps, npoint)
      !**********************************************************************************************************************
      !*** DUMMY ARGUMENTS
      double precision, intent(out) :: eps                                  !*** geometric tolerance
      integer, intent(out) :: npoint                                        !*** vertex number
      !
      !*** LOCAL VARIABLES
      double precision, dimension(3) :: vertex                              !*** cone vertex
      double precision, dimension(3) :: axis                                !*** cone axis, cylinder axis
      double precision :: alpha                                             !*** cone half-aperture
      double precision, dimension(3) :: center                              !*** sphere center, cylinder center
      double precision :: radius                                            !*** sphere radius, cylinder radius
      double precision, dimension(3) :: r_line, n_line                      !*** line data
      double precision, dimension(3) :: r_plane, n_plane                    !*** plane data
      double precision :: dist                                              !*** distance between point and algebraic surface
      double precision, dimension(3) :: p1, p2                              !*** auxiliary variables
      integer :: status                                                     !*** flag
      !**********************************************************************************************************************
      !
      eps = 0.d0
      !
      triangles5: do nt = 1, NRTRIAN                                        !*** loop over triangles
      !
        if (TRIANGLES(nt)%Type .ne. 'PTITri') cycle                         !*** IF triangle does not belong to an algebraic surface cycle
      !
        ns = TRIANGLES(nt)%Idata(1)                                         !*** get surface number
      !
        select case (SURFACES(ns)%Type)                                     !*** select case for surface type
      !
        case ('Sphere')
          center(1:3) = SURFACES(ns)%Rdata(1:3)                             !*** get sphere data
          radius = SURFACES(ns)%Rdata(4)
      !
!          v_aux(1) = 0.d0; v_aux(2) = 0.d0; v_aux(3) = radius               !*** set up vector v_aux
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)


              if (nv .eq. 17) cycle


            v(1:3) = POINTS(nv)%Rdata(1:3) - center(1:3)                    !*** get vertex coordinates wrt center
      !
            dist = abs(sqrt(v(1)**2 + v(2)**2 + v(3)**2) - radius)          !*** compute distance to algebraic surface
      !
            if (dist .gt. eps) then
              eps = dist
              npoint = nv
            end if
      !
          end do
      !
        case('Cylinder')
          center(1:3) = SURFACES(ns)%Rdata(1:3)                             !*** get cylinder data
          axis(1:3) = SURFACES(ns)%Rdata(4:6)
          radius = SURFACES(ns)%Rdata(7)
      !
          call rotation(axis, e3, RBM)                                      !*** align axis to z-axis
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)

              if (nv .eq. 17) cycle

            v(1:3) = POINTS(nv)%Rdata(1:3) - center(1:3)                    !*** get vertex coordinates wrt center
      !
            v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3))                           !*** apply rotation to vertex
      !
            dist = abs(sqrt(v(1)**2 + v(2)**2) - radius)                    !*** compute distance to algebraic surface
      !
            if (dist .gt. eps) then
              eps = dist
              npoint = nv
            end if
      !
          end do
      !
        case('Cone')
          vertex(1:3) = SURFACES(ns)%Rdata(1:3)                             !*** get cone data
          axis(1:3) = SURFACES(ns)%Rdata(4:6)
          alpha = SURFACES(ns)%Rdata(7)
      !
          call rotation(axis, e3, RBM)                                      !*** align axis to z-axis
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)

              if (nv .eq. 17) cycle

            v(1:3) = POINTS(nv)%Rdata(1:3) - vertex(1:3)                    !*** get vertex coordinates wrt vertex
      !
            v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3))                           !*** apply rotation to vertex
      !
            dist = abs(sqrt(v(1)**2 + v(2)**2) - tan(alpha)*abs(v(3)))      !*** compute distance to algebraic surface
!            write(*,*)'cone dist = ', dist
!            if (nv .eq. 136) call pause
      !
            if (dist .gt. eps) then
              eps = dist
              npoint = nv
            end if
      !
          end do
      !
        case('VecPt')
          r_plane(1:3) = SURFACES(ns)%Rdata(1:3)
          n_plane(1:3) = SURFACES(ns)%Rdata(4:6)
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)


              if (nv .eq. 17) cycle

            v(1:3) = POINTS(nv)%Rdata(1:3)                                  !*** get vertex coordinates
            call planeANDline(n_plane, r_plane, n_plane, v, v_aux)          !*** project prependicularly onto plane
            v_aux(1:3) = v(1:3) - v_aux(1:3)
      !
            dist = sqrt(v_aux(1)**2 + v_aux(2)**2 + v_aux(3)**2)            !*** compute distance to algebraic surface
      !
            if (dist .gt. eps) then
              eps = dist
              npoint = nv
            end if
      !
          end do
      !
        case('ThreePt')
          r_plane(1:3) = SURFACES(ns)%Rdata(1:3)
          p1(1:3) = SURFACES(ns)%Rdata(4:6) - r_plane(1:3)
          p2(1:3) = SURFACES(ns)%Rdata(7:9) - r_plane(1:3)
          call cross_product(p1, p2, n_plane)
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)

              if (nv .eq. 17) cycle

            v(1:3) = POINTS(nv)%Rdata(1:3)                                  !*** get vertex coordinates
            call planeANDline(n_plane, r_plane, n_plane, v, v_aux)          !*** project prependicularly onto plane
            v_aux(1:3) = v(1:3) - v_aux(1:3)
      !
            dist = sqrt(v_aux(1)**2 + v_aux(2)**2 + v_aux(3)**2)            !*** compute distance to algebraic surface
      !
            if (dist .gt. eps) then
              eps = dist
              npoint = nv
            end if
      !
          end do
      !
        end select
      !
      end do triangles5
      !
      END SUBROUTINE geometric_tolerance
      !
      !
      !**********************************************************************************************************************
      SUBROUTINE projectONTOplanes
      !**********************************************************************************************************************
      !*** LOCAL VARIABLES
      double precision, dimension(3) :: r_line, n_line                      !*** line data
      double precision, dimension(3) :: r_plane, n_plane                    !*** plane data
      double precision, dimension(3) :: p1, p2                              !*** auxiliary variables
      !**********************************************************************************************************************
      !
      triangles1: do nt = 1, NRTRIAN                                        !*** project onto PLANES
      !
        if (TRIANGLES(nt)%Type .ne. 'PTITri') cycle                         !*** IF triangle does not belong to an algebraic surface cycle
      !
        ns = TRIANGLES(nt)%Idata(1)                                         !*** get surface number
      !
        if (SURFACES(ns)%Type .eq. 'VecPt') then
          r_plane(1:3) = SURFACES(ns)%Rdata(1:3)                            !*** get plane data
          n_plane(1:3) = SURFACES(ns)%Rdata(4:6)
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)                                    !*** get vertex coordinates
            v(1:3) = POINTS(nv)%Rdata(1:3)
      !
            call planeANDline(n_plane, r_plane, n_plane, v, v_aux)          !*** project prependicularly onto plane
            v(1:3) = v_aux(1:3)                                             !*** store result into v
      !
            call store_vertex(nv, v, ns)
          end do
      !
        else if (SURFACES(ns)%Type .eq. 'ThreePt') then
          r_plane(1:3) = SURFACES(ns)%Rdata(1:3)                            !*** get plane data
          p1(1:3) = SURFACES(ns)%Rdata(4:6) - r_plane(1:3)
          p2(1:3) = SURFACES(ns)%Rdata(7:9) - r_plane(1:3)
          call cross_product(p1, p2, n_plane)
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)                                    !*** get vertex coordinates
            v(1:3) = POINTS(nv)%Rdata(1:3)
      !
            call planeANDline(n_plane, r_plane, n_plane, v, v_aux)          !*** project prependicularly onto plane
            v(1:3) = v_aux(1:3)                                             !*** store result into v
      !
            call store_vertex(nv, v, ns)
          end do
        end if
      !
      end do triangles1
      !
      END SUBROUTINE projectONTOplanes
      !
      !
      !**********************************************************************************************************************
      SUBROUTINE projectONTOcones
      !**********************************************************************************************************************
      !*** LOCAL VARIABLES
      double precision, dimension(3) :: vertex                              !*** cone vertex
      double precision, dimension(3) :: axis                                !*** cone axis
      double precision :: alpha                                             !*** cone half-aperture
      !**********************************************************************************************************************
      !
      triangles2: do nt = 1, NRTRIAN                                        !*** project onto CONES
      !
        if (TRIANGLES(nt)%Type .ne. 'PTITri') cycle                         !*** IF triangle does not belong to an algebraic surface cycle
      !
        ns = TRIANGLES(nt)%Idata(1)                                         !*** get surface number
      !
        if (SURFACES(ns)%Type .eq. 'Cone') then
      !
          vertex(1:3) = SURFACES(ns)%Rdata(1:3)                             !*** get cone data
          axis(1:3) = SURFACES(ns)%Rdata(4:6)
          alpha = SURFACES(ns)%Rdata(7)
      !
          call rotation(axis, e3, RBM_aux)                                  !*** align axis to z-axis
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)
            v(1:3) = POINTS(nv)%Rdata(1:3) - vertex(1:3)                    !*** get vertex coordinates wrt vertex
      !
            v(1:3) = MATMUL(RBM_aux(1:3,1:3), v(1:3))                       !*** apply rotation to vertex
      !
            v_aux(1) = v(1); v_aux(2) = v(2); v_aux(3) = 0.d0               !*** project v onto xy-plane
      !
            call rotation(v_aux, e1, RBM)                                   !*** rotation that moves vertex onto xz-plane
      !
            v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3))                           !*** apply rotation to vertex
      !
            v(1) = tan(alpha)*abs(v(3))                                     !*** project perpendicularly to cone axis
      !
            RBM(1:3,1:3) = MATMUL(RBM(1:3,1:3), RBM_aux(1:3,1:3))           !*** determine rigid body motion
            RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                          !*** invert rigid body motion
      !
            v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3)) + vertex(1:3)             !*** apply inverse rigid body motion
      !
            call store_vertex(nv, v, ns)
          end do
      !
        end if
      !
      end do triangles2
      !
      END SUBROUTINE projectONTOcones
      !
      !
      !**********************************************************************************************************************
      SUBROUTINE projectONTOspheres
      !**********************************************************************************************************************
      !*** LOCAL VARIABLES
      double precision, dimension(3) :: center                              !*** sphere center
      double precision :: radius                                            !*** sphere radius
      !**********************************************************************************************************************
      !
      triangles3: do nt = 1, NRTRIAN                                        !*** project onto SPHERES
      !
        if (TRIANGLES(nt)%Type .ne. 'PTITri') cycle                         !*** IF triangle does not belong to an algebraic surface cycle
      !
        ns = TRIANGLES(nt)%Idata(1)                                         !*** get surface number
      !
        if (SURFACES(ns)%Type .eq. 'Sphere') then
      !
          center(1:3) = SURFACES(ns)%Rdata(1:3)                             !*** get sphere data
          radius = SURFACES(ns)%Rdata(4)
      !
          v_aux(1) = 0.d0; v_aux(2) = 0.d0; v_aux(3) = radius               !*** set up vector v_aux
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)
            v(1:3) = POINTS(nv)%Rdata(1:3) - center(1:3)                    !*** get vertex coordinates wrt center
      !
            call rotation(v_aux, v, RBM)                                    !*** determine rigid body motion
      !
            v(1:3) = MATMUL(RBM(1:3,1:3), v_aux(1:3)) + center(1:3)         !*** apply rigid body motion
      !
            call store_vertex(nv, v, ns)
          end do
      !
        end if
      !
      end do triangles3
      !
      END SUBROUTINE projectONTOspheres
      !
      !
      !**********************************************************************************************************************
      SUBROUTINE projectONTOcylinders
      !**********************************************************************************************************************
      !*** LOCAL VARIABLES
      double precision, dimension(3) :: center                              !*** cylinder center
      double precision, dimension(3) :: axis                                !*** cylinder axis
      double precision :: radius                                            !*** cylinder radius
      double precision :: rho, theta, z                                     !*** cylindrical coordinates
      double precision :: sign
      integer :: nsphere                                                    !*** sphere surface number
      double precision :: radius_aux                                        !*** auxiliary variable
      !**********************************************************************************************************************
      !
      triangles4: do nt = 1, NRTRIAN                                        !*** project onto CYLINDERS
      !
        if (TRIANGLES(nt)%Type .ne. 'PTITri') cycle                         !*** IF triangle does not belong to an algebraic surface cycle
      !
        ns = TRIANGLES(nt)%Idata(1)                                         !*** get surface number
      !
        if (SURFACES(ns)%Type .eq. 'Cylinder') then
      !
          center(1:3) = SURFACES(ns)%Rdata(1:3)                             !*** get cylinder data
          axis(1:3) = SURFACES(ns)%Rdata(4:6)
          radius = SURFACES(ns)%Rdata(7)
      !
          call rotation(axis, e3, RBM_aux)                                  !*** align axis to z-axis
      !
          do i = 1, 3                                                       !*** loop over triangle vertices
            nv = TRIANGLES(nt)%VertNo(i)
            v(1:3) = POINTS(nv)%Rdata(1:3) - center(1:3)                    !*** get vertex coordinates wrt center
      !
            B_status = associated(POINTS(nv)%Idata)                         !*** check whether Idata is associated with a target
      !
            if (.not. B_status) then                                        !*** IF first visit, i.e. Idata is not associated
      !
              v(1:3) = MATMUL(RBM_aux(1:3,1:3), v(1:3))                     !*** apply rotation to vertex
      !
              v_aux(1) = v(1); v_aux(2) = v(2); v_aux(3) = 0.d0             !*** project v onto xy-plane
      !
              call rotation(v_aux, e1, RBM)                                 !*** rotation that moves vertex onto xz-plane
      !
              v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3))                         !*** apply rotation to vertex
      !
              v(1) = radius                                                 !*** set x-component equal to cylinder radius
      !
              RBM(1:3,1:3) = MATMUL(RBM(1:3,1:3), RBM_aux(1:3,1:3))         !*** determine rigid body motion
              RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                        !*** invert rigid body motion
      !
              v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3)) + center(1:3)           !*** apply inverse rigid body motion
      !
              call store_vertex(nv, v, ns)
      !
            else                                                            !*** ELSE vertex has already been visited
              if (POINTS(nv)%Idata(1) .ne. 0) then                          !***** IF previously visited surface was a sphere
                nsphere = POINTS(nv)%Idata(1)                               !***** get sphere number
                radius_aux = SURFACES(nsphere)%Rdata(4)                     !***** get sphere radius
      !
                v(1:3) = MATMUL(RBM_aux(1:3,1:3), v(1:3))                   !***** apply rotation to vertex
      !
                call cartesian2cylindrical(v, rho, theta, z)                !***** get cylindrical coordinates (theta is the only quantity of interest)
      !
                if (z.ge. 0.d0) then                                        !***** determine sign(z)
                  sign = 1.d0
                else
                  sign = -1.d0
                end if
      !
                rho = radius                                                !***** update rho and z
                z = sign*sqrt(radius_aux**2 - radius**2)
      !
                call cylindrical2cartesian(rho, theta, z, v)                !***** determine new coordinates
      !
                RBM(1:3,1:3) = TRANSPOSE(RBM_aux(1:3,1:3))                  !***** invert rigid body motion
                v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3)) + center(1:3)         !***** apply inverse rigid body motion
      !
                call store_vertex(nv, v, ns)
      !
              else                                                          !***** ELSE previously visited surface was NOT a sphere
                v(1:3) = MATMUL(RBM_aux(1:3,1:3), v(1:3))                   !***** apply rotation to vertex
      !
                v_aux(1) = v(1); v_aux(2) = v(2); v_aux(3) = 0.d0           !***** project v onto xy-plane
      !
                call rotation(v_aux, e1, RBM)                               !***** rotation that moves vertex onto xz-plane
      !
                v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3))                       !***** apply rotation to vertex
      !
                v(1) = radius                                               !***** set x-component equal to cylinder radius
      !
                RBM(1:3,1:3) = MATMUL(RBM(1:3,1:3), RBM_aux(1:3,1:3))       !***** determine rigid body motion
                RBM(1:3,1:3) = TRANSPOSE(RBM(1:3,1:3))                      !***** invert rigid body motion
      !
                v(1:3) = MATMUL(RBM(1:3,1:3), v(1:3)) + center(1:3)         !***** apply inverse rigid body motion
      !
                call store_vertex(nv, v, ns)
      !
              end if                                                        !***** END IF
            end if                                                          !*** END IF
          end do
      !
        end if
      !
      end do triangles4
      !
      END SUBROUTINE projectONTOcylinders
      !**********************************************************************************************************************
      !
      END SUBROUTINE upgrade2double_OLD
      !
      !
      !
      !
      !
      !**********************************************************************************************************************
      SUBROUTINE store_vertex(nv, v, ns)
      !**********************************************************************************************************************
      !*** LATEST REVISION:     Aug 08
      !
      !*** PURPOSE:             routine stores vertex coodinates into into the GMP data structure. If vertex belongs to a
      !***                      spherical surface it stores the sphere number into Idata(1), otherewise it stores 0.
      !
      !*** REQUIRED PROCEDURES: none
      !**********************************************************************************************************************
      use GMP
      !**********************************************************************************************************************
      IMPLICIT NONE
      !**********************************************************************************************************************
      !*** DUMMY ARGUMENTS
      integer, intent(in) :: nv                                           !*** vertex number
      double precision, dimension(3), intent(in) :: v                     !*** vertex coordinates
      integer :: ns                                                       !*** sphere number
      !
      !*** LOCAL VARIABLES
      integer :: status
      logical :: B_status                                                 !*** boolean flag
      !**********************************************************************************************************************
      !
      POINTS(nv)%Rdata(1:3) = v(1:3)                                      !*** store new coordinates into Rdata
      !
      B_status = associated(POINTS(nv)%Idata)                             !*** check whether Idata is associated with a target
      !
      if (.not. B_status) then                                            !*** IF first visit, i.e. Idata is not associated
        allocate (POINTS(nv)%Idata(1), STAT = status)                     !*** allocate Idata
      !
        if (status .ne. 0) then
          write(*,*)'---------------------------------'
          write(*,*)'store_vertex:                    '
          write(*,*)'ERROR: Idata array not allocated.'
          write(*,*)'---------------------------------'
          stop
        end if
      !
      end if                                                              !*** END IF
      !
      POINTS(nv)%Idata(1) = 0                                             !*** set Idata equal to 0
      !
      if (SURFACES(ns)%Type .eq. 'Sphere') POINTS(nv)%Idata(1) = ns       !*** IF on sphere set Idata equal to surface number
      !
      END SUBROUTINE store_vertex
