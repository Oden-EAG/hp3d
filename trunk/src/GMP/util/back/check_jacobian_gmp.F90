!-----------------------------------------------------------------------
!> Purpose : check of jacobian for GMP blocks      
!!
!! @revision Apr 11      
!-----------------------------------------------------------------------
!
subroutine check_jacobian_gmp
!
      use control
      use GMP
!-----------------------------------------------------------------------
      implicit none
      common /ctetra_TraTet/ iprint_tetra_TraTet
      common /cprism/ iprint_prism
      common /ctrian_PTITri/ iprint_trian_PTITri
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
      double precision, dimension(2,5000) :: x_tri                        ! 
      double precision, dimension(3,5000) :: x_tet                        ! 
      double precision, dimension(3,5000) :: x_pri                        ! 
      double precision, dimension(3,5000) :: x_pyr                        !  
      integer                            :: n_tri,n_tet,n_pri,n_pyr
      double precision, dimension(3)     :: x, prod_new,prod_old
      double precision, dimension(3,3)   :: dx,dxdeta
      double precision, dimension(3)     :: r
!      double precision, dimension(3,2)   :: dr
      double precision                   :: jac, prod
      integer                            :: i,j,k,m,n
      integer                            :: neg_tri,neg_tet
      integer                            :: neg_pri,neg_pyr
      integer, dimension(NRSURFS)        :: n_surface
      integer, dimension(NRDOMAIN)       :: n_domain
      integer                            :: nd, ns, np, np0, iv
      integer                            :: iflag, ierror
      integer :: iprint_tetra_TraTet, iprint_prism,iprint_trian_PTITri
      double precision        :: tet_min_jac, pri_min_jac, pyr_min_jac
      double precision        :: dr(3,3),s1,s2
!-----------------------------------------------------------------------
!     PARAMETERS
      double precision, parameter      :: eps = 1.d-13                   ! geometric tolerance 
      integer, parameter :: nsub = 5                                     ! number of subdivisions
!-----------------------------------------------------------------------
!     generate grid on triangle
!
      n=0
      do j=0,nsub
        do i=0,(nsub-j)
          if ((i.eq.0)    .AND. (j.eq.0))    cycle
          if ((i.eq.nsub) .AND. (j.eq.0))    cycle 
          if ((i.eq.0)    .AND. (j.eq.nsub)) cycle
          n=n+1
          x_tri(1,n) = dble(i)/dble(nsub)
          x_tri(2,n) = dble(j)/dble(nsub)
          if ((x_tri(1,n)+x_tri(2,n)).gt.1.d0) then
            x_tri(1:2,n) = (1.d0-GEOM_TOL)*x_tri(1:2,n)
          endif 
        enddo
      enddo
      n_tri=n
!-----------------------------------------------------------------------
!     generate grid on tet
!
      n = 0
      do k = 0, nsub
        do j = 0, (nsub-k)
          do i = 0, (nsub-k-j)
            if ((i.eq.0)    .AND. (j.eq.0)    .AND. (k.eq.0))    cycle   ! 1st vertex
            if ((i.eq.nsub) .AND. (j.eq.0)    .AND. (k.eq.0))    cycle   ! 2nd vertex
            if ((i.eq.0)    .AND. (j.eq.nsub) .AND. (k.eq.0))    cycle   ! 3rd vertex
            if ((i.eq.0)    .AND. (j.eq.0)    .AND. (k.eq.nsub)) cycle   ! 4th vertex
            n = n + 1
            x_tet(1,n) = dble(i)/dble(nsub)
            x_tet(2,n) = dble(j)/dble(nsub)
            x_tet(3,n) = dble(k)/dble(nsub)
            if ((x_tet(1,n) + x_tet(2,n) + x_tet(3,n)) .gt. 1.d0) then
              x_tet(1:3,n) = (1.d0 - GEOM_TOL)*x_tet(1:3,n)
            end if 
          end do
        end do
      end do
      n_tet = n
!-----------------------------------------------------------------------
!     generate grid on prism
!
      n = 0
      do k = 0, nsub
        do j = 0, nsub
          do i = 0, (nsub-j)
            if ((i.eq.0)    .AND. (j.eq.0))       cycle   ! 1st vertical edge
            if ((i.eq.nsub) .AND. (j.eq.0))       cycle   ! 2nd vertical edge
            if ((i.eq.0)    .AND. (j.eq.nsub))    cycle   ! 3rd vertical edge
            n = n + 1
            x_pri(1,n) = dble(i)/dble(nsub)
            x_pri(2,n) = dble(j)/dble(nsub)
            x_pri(3,n) = dble(k)/dble(nsub) 
            if ((x_pri(1,n) + x_pri(2,n)) .gt. 1.d0) then
              x_pri(1:2,n) = (1.d0 - GEOM_TOL)*x_pri(1:2,n)
            end if 
            if (x_pri(3,n) .gt. 1.d0) then
              x_pri(3,n) = (1.d0 - GEOM_TOL)*x_pri(3,n)
            end if
          end do
        end do
      end do
      n_pri = n
!
!-----------------------------------------------------------------------
!     generate grid on pyramid
      n = 0
      do k = 0, nsub
        do j = 0, (nsub-k)
          do i = 0, (nsub-k)
            if ((i.eq.0)    .AND. (j.eq.0)    .AND. (k.eq.0))    cycle   ! 1st vertex
            if ((i.eq.nsub) .AND. (j.eq.0)    .AND. (k.eq.0))    cycle   ! 2nd vertex
            if ((i.eq.nsub) .AND. (j.eq.nsub) .AND. (k.eq.0))    cycle   ! 3rd vertex
            if ((i.eq.0)    .AND. (j.eq.nsub) .AND. (k.eq.0))    cycle   ! 4th vertex
            if ((i.eq.0)    .AND. (j.eq.0)    .AND. (k.eq.nsub)) cycle   ! 5th vertex
            n = n + 1
            x_pyr(1,n) = dble(i)/dble(nsub)
            x_pyr(2,n) = dble(j)/dble(nsub)
            x_pyr(3,n) = dble(k)/dble(nsub) 
          end do
        end do
      end do
      n_pyr = n
!
!-----------------------------------------------------------------------
!
      n_surface = 0
!
      neg_tri = 0
      prod_new = 0.d0
      prod_old = 0.d0

!  ...T R I A N G L E S................................................      
      neg_tri=0 
      ! loop over triangles
      do k=1,NRTRIAN
        if (TRIANGLES(k)%Type.eq.'PlaneTri') cycle    

        ! loop over triangle grid
        iflag=0
        do i=1,n_tri
          call trian(k,x_tri(1:2,n), x,dxdeta(1:3,1:2))
          ! surface jacobian
          call cross_product(dxdeta(1:3,1),dxdeta(1:3,2), prod_old)
          call norm(prod_old, s1)
          if (s1.le.GEOM_TOL) then
            iflag=1
            write(*,7001) k,TRIANGLES(k)%Type
 7001       format('check_jacobian_gmp: DEGENERATE POINT IN TRIANGLE =', &
                    i5,' Type = ',a10)
          endif
        ! loop over triangle grid
        enddo
       
        ! check error flag
        if (iflag.eq.1) then
          neg_tri=neg_tri+1
          if (TRIANGLES(k)%Type.eq.'PTITri') then
            ns=TRIANGLES(k)%Idata(1) ; n_surface(ns)=n_surface(ns)+1  
          endif
        endif

      ! loop over triangles
      enddo



!!!        n = 0
!!!        do j = 0, nsub
!!!          do i=0,nsub-j
!!! 100        continue
!!!            if ((i.eq.0) .AND. (j.eq.0)) then                            ! skip 1st vertex
!!!              iflag = 0
!!!              cycle
!!!            end if
!!!            if ((i.eq.nsub) .AND. (j.eq.0)) then                         ! skip 2nd vertex
!!!              iflag = 0
!!!              cycle
!!!            end if
!!!            if ((i.eq.0) .AND. (j.eq.nsub)) then                         ! skip 3rd vertex
!!!              exit middle
!!!            end if
!!!            n = n + 1
!!!            call trian(k,x_tri(1:2,n), r,dr)
!!!            if (iflag.eq.0) then                                         ! IF starting a new row
!!!              call cross_product(dr(1:3,1),dr(1:3,2), prod_old)
!!!              call norm(prod_old, s1)
!!!              if (s1.le.GEOM_TOL) then
!!!                write(*,7001) k,TRIANGLES(k)%Type
!!! 7001      format('check_jacobian_gmp: DEGENERATE POINT IN TRIANGLE =',
!!!     .                  i5,' Type = ',a10)
!!!                call pause
!!!                iprint_trian_PTITri=1
!!!                go to 100
!!!              endif
!!!            else                                                         ! ELSE row already started
!!!              call cross_product(dr(1:3,1),dr(1:3,2), prod_new)
!!!              call norm(prod_new, s2)
!!!              if (s2.le.GEOM_TOL) then
!!!                write(*,7001) k,TRIANGLES(k)%Type
!!!                call pause
!!!                iprint_trian_PTITri=1
!!!                go to 100
!!!              endif
!!!              call scalar_product(prod_old,prod_new, prod)
!!!              prod=prod/s1/s2
!!!              prod_old = prod_new; s1=s2
!!!              if (prod .le. 0.5d0) then                                   ! IF angle is greater than pi/2 ???
!!!                neg_tri = neg_tri + 1
!!!                if (TRIANGLES(k)%Type.eq.'PTITri') then
!!!                  ns = TRIANGLES(k)%Idata(1)
!!!                  n_surface(ns) = n_surface(ns) + 1  
!!!                endif
!!!                exit middle
!!!              end if
!!!            end if 
!!!            iflag = 1                                                    ! raise iflag
!!!            if ((i+j) .eq. nsub) iflag = 0                               ! IF last point in row 
!!!          end do
!!!        end do
!!!      end do 
!
!-----------------------------------------------------------------------
!      do i = 1, n_tri
!        write (*,*)'tri = ', x_tri(1:2,i)
!      end do
!
!      do i = 1, n_tet
!        write (*,*)'tet = ', x_tet(1:3,i)
!      end do
!      call pause
!
!       do i = 1, n_pri
!        write (*,*)'pri = ', x_pri(1:3,i)
!      end do
!
!       do i = 1, n_pyr
!        write (*,*)'pyr = ', x_pyr(1:3,i)
!      end do
!-----------------------------------------------------------------------
!
      n_domain = 0
!
!  ...T E T R A H E D R A.................................................
      tet_min_jac = 1.d30
      neg_tet = 0
      do i=1,NRTETRA
        if (TETRAS(i)%Type .eq. 'Linear') cycle
        do j=1,n_tet
 121      continue
!cc          write(*,*) 'i,x_tet(1:3,j) = ',i,x_tet(1:3,j)
          call tetra(i,x_tet(1:3,j), x,dx)
          jac = dx(1,1)*dx(2,2)*dx(3,3)  &
              + dx(2,1)*dx(3,2)*dx(1,3)  &
              + dx(1,2)*dx(2,3)*dx(3,1)  &
              - dx(3,1)*dx(2,2)*dx(1,3)  &
              - dx(2,1)*dx(1,2)*dx(3,3)  &
              - dx(1,1)*dx(3,2)*dx(2,3)
          tet_min_jac = min(tet_min_jac,jac)
          if (jac .lt. 0.d0) then
            write(*,*) 'check_jacobian_gmp: jacobian = ',jac
            nd = TETRAS(i)%Domain
            n_domain(nd) = n_domain(nd) + 1  
            neg_tet = neg_tet + 1
            write(*,7024) x_tet(1:3,j)
 7024       format('check_jacobian_gmp: x_tet = ',3f8.3)
            write(*,*)'check_jacobian_gmp: neg jac for tet ',i
!cc            call pause
!
!  .........sanity check at vertices
            np0 = TETRAS(i)%VertNo(1)
            do iv=1,3
              np = TETRAS(i)%VertNo(1+iv)
              dr(1:3,iv) = POINTS(np)%Rdata(1:3)-POINTS(np0)%Rdata(1:3)
            enddo
            call mixed_product(dr(1:3,1),dr(1:3,2),dr(1:3,3), jac)
            write(*,*) 'check_jacobian_gmp: vertex jacobian = ',jac
!cc            call graphg
            call print_GMP
            iprint_tetra_TraTet=1
            go to 121
            exit
!            write(*,*)'negative jacobian for tetra = ',i
          end if
        end do
      end do
!
      pri_min_jac = 1.d30
      neg_pri = 0   
      do i = 1, NRPRISM
        if (PRISMS(i)%Type .eq. 'Linear') cycle
        do j = 1, n_pri
 151      continue
          call prism(i,x_pri(1:3,j), x,dx)
          jac = dx(1,1)*dx(2,2)*dx(3,3) + &
                dx(2,1)*dx(3,2)*dx(1,3) + &
                dx(1,2)*dx(2,3)*dx(3,1) - &
                dx(3,1)*dx(2,2)*dx(1,3) - &
                dx(2,1)*dx(1,2)*dx(3,3) - &
                dx(1,1)*dx(3,2)*dx(2,3)
          pri_min_jac = min(pri_min_jac,jac)
          if (jac .lt. 0.d0) then
            iprint_prism=1
            write(*,*)'negative jacobian for prism = ',i
            write(*,*)'xi = ',x_pri(1:3,j)
            go to 151
            nd = PRISMS(i)%Domain
            n_domain(nd) = n_domain(nd) + 1  
            neg_pri = neg_pri + 1
            exit
          end if
        end do     
      end do
!
      pyr_min_jac = 1.d30
      neg_pyr = 0
      do i = 1, NRPYRAM
        if (PYRAMIDS(i)%Type .eq. 'Linear') cycle
        do j = 1, n_pyr
          call pyram(i,x_pyr(1:3,j), x,dx)
          jac = dx(1,1)*dx(2,2)*dx(3,3) + &
                dx(2,1)*dx(3,2)*dx(1,3) + &
                dx(1,2)*dx(2,3)*dx(3,1) - &
                dx(3,1)*dx(2,2)*dx(1,3) - &
                dx(2,1)*dx(1,2)*dx(3,3) - &
                dx(1,1)*dx(3,2)*dx(2,3)
          pyr_min_jac = min(pyr_min_jac,jac)
          if (jac .lt. 0.d0) then
            nd = PYRAMIDS(i)%Domain
            n_domain(nd) = n_domain(nd) + 1  
            neg_pyr = neg_pyr + 1
            exit
!            write(*,*)'negative jacobian for pyramid = ',i
          end if
        end do
      end do
!
      write(*,*)'points per tet = ',n_tet
      write(*,*)'points per pri = ',n_pri
      write(*,*)'points per pyr = ',n_pyr
      write(*,*)'tets with neg jac  = ', neg_tet
      write(*,*)'pris with neg jac = ', neg_pri
      write(*,*)'pyrs with neg jac = ', neg_pyr
      write(*,*)'min jac for tet = ',tet_min_jac
      write(*,*)'min jac for pri = ',pri_min_jac
      write(*,*)'min jac for pyr = ',pyr_min_jac
      write(*,*)'elements with neg jac per domain'
      do i = 1, NRDOMAIN
        write(*,*)i,n_domain(i)
      end do
      write(*,*)'--------------------------------------'
      write(*,*)'points per tri = ',n_tri
      write(*,*)'tris with neg jac  = ', neg_tri 
      write(*,*)'triangles with neg jac per surface'
      do i = 1, NRSURFS
        write(*,*)i,n_surface(i)
      end do
!
      call pause
!      call print_GMP
!
      end subroutine check_jacobian_gmp
