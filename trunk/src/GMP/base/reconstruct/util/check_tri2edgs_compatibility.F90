!--------------------------------------------------------------------------
!
      subroutine check_tri2edgs_compatibility(Nt)
!
!--------------------------------------------------------------------------
      use GMP
      use control
      use element_data
!--------------------------------------------------------------------------
      implicit none
!--------------------------------------------------------------------------
!     DUMMY ARGUMENTS
      integer,intent(in)      :: Nt
!--------------------------------------------------------------------------
!     LOCAL VARIBLES
!  ...master edge coordinate      
      real*8                :: zeta
!  ...master triangle coordiante      
      real*8,dimension(2)   :: eta
!  ...physical space coordinates      
      real*8,dimension(3)   :: x1,x2
!  ...derivatives      
      real*8,dimension(2)   :: deta_dzeta      
      real*8,dimension(3,2) :: dx2_deta
      real*8,dimension(3)   :: dx1_dzeta,dx2_dzeta 
!  ...miscellanea
      real*8                :: s
      real*8,dimension(3)   :: temp
      integer                 :: nc,ie,nsub,i,j,iflag
!--------------------------------------------------------------------------
!
!  ...set number of subdivisions      
      nsub = 10
!      
!  ...loop over edges
      do ie = 1,3      
        nc = TRIANGLES(Nt)%EdgeNo(ie)
!  ...loop over subdivisions
        do i = 0,nsub
          zeta = i*1.d0/nsub
!
!  .......1st TERM OF COMPARISON
          if (nc.gt.0) then
            call curve(nc,     zeta,      x1,dx1_dzeta)
          else
            call curve(abs(nc),1.d0-zeta, x1,dx1_dzeta)
            dx1_dzeta = -dx1_dzeta
          endif
!          write(*,*)'x1 = ',x1
!          write(*,*)'dx1 = ',dx1_dzeta
!
!  .......2ND TERM OF COMPARISON
          call edge_param('trian',ie,zeta, eta,deta_dzeta)
          call trian(Nt,eta, x2,dx2_deta)
          dx2_dzeta = dx2_deta(1:3,1)*deta_dzeta(1) +  &
                      dx2_deta(1:3,2)*deta_dzeta(2)     
!
!  .......COMPARE
          iflag = 0
          call norm(x1-x2, s)
          if (s.gt.GEOM_TOL) then
            iflag = 1
            write(*,1000)
 1000       format(' check_tri2edgs_compatibility: inconsistency b/w function values')
            write(*,1001) Nt,ie,nc,zeta
 1001       format('   Nt = ',i8,'; ie = ',i1,'; nc = ',i5,'; zeta = ',e12.5)
            write(*,1002) x1
 1002       format('   curve:    x        = ',3(e12.5,2x))
            write(*,1003) x2
 1003       format('   triangle: x        = ',3(e12.5,2x))
          endif
!          
          call norm(dx1_dzeta-dx2_dzeta, s)
          if (s.gt.GEOM_TOL) then
            iflag = 1
            write(*,1004)
 1004       format(' check_tri2edgs_compatibility: inconsistency b/w derivatives')
            write(*,1005) Nt,ie,nc,zeta
 1005       format('   Nt = ',i8,'; ie = ',i1,'; nc = ',i5,'; zeta = ',e12.5)
            write(*,1006) dx1_dzeta
 1006       format('   curve:    dx_dzeta = ',3(e12.5,2x))
            write(*,1007) dx2_dzeta
 1007       format('   triangle: dx_dzeta = ',3(e12.5,2x))
            call cross_product(dx1_dzeta,dx2_dzeta, temp)
            call norm(temp, s)
            write(*,1008) s
 1008       format('           |vec prod| = ',e12.5)
            do j = 1,3
              if (abs(dx2_dzeta(j)).lt.GEOM_TOL) cycle
              dx2_dzeta(j) = dx1_dzeta(j)/dx2_dzeta(j)
            enddo
            write(*,1009) dx2_dzeta
 1009       format('                ratio = ',3(e12.5,2x))            
          endif
!          
!  .....end of loop over subdivisions
        enddo
        if (iflag.eq.1) then
          call pause
        endif
!  ...end of loop over edges
      enddo
!
      end
