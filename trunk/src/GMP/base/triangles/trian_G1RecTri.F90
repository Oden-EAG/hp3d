!----------------------------------------------------------------------
!      
      subroutine trian_G1RecTri(Nt,Eta, X,dXdEta)
!
!----------------------------------------------------------------------
!      
!   latest revision    - Feb 10
!
!   purpose            - routine supports G1 reconstruction for a
!                        triangle
!   arguments :
!     in:
!               No     - a GMP triangle number
!               Eta    - reference coordinates of a point
!                        in the triangle
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        to the parameters
!
!----------------------------------------------------------------------
!     MODULES      
      use GMP
!----------------------------------------------------------------------
      implicit none
!----------------------------------------------------------------------
!     DUMMY ARGUMENTS      
      integer,                 intent(in)  :: Nt
      real(8), dimension(2),   intent(in)  :: Eta
      real(8), dimension(3),   intent(out) :: X
      real(8), dimension(3,2), intent(out) :: dXdEta
!----------------------------------------------------------------------
!     LOCAL VARIABLES
      integer               :: i,j,k
      real(8)               :: poly
      real(8), dimension(2) :: dpoly
      real(8), dimension(3) :: nor
      integer               :: iprint
      integer,parameter     :: deg = 7
!----------------------------------------------------------------------
!     EXTERNAL FUNCTIONS      
      integer, external     :: bijec
!----------------------------------------------------------------------
!
      iprint=0
!
      if (iprint.eq.1) then
        write(*,1000)Nt,Eta
 1000   format(' trian_G1RecTri: Nt = ',i7,'; Eta = ',2(e12.5,2x))
      endif
!
      X = 0.d0
!  ...accumulate
      k = 0
      do j = 0,deg
        do i = 0,(deg-j) 
          call biv_bernstein_poly(i,j,deg,Eta(1),Eta(2), poly,dpoly)
          X = X + poly*TRIANGLES(Nt)%Rdata(bijec(i,j):bijec(i,j)+2)
          if (iprint.ge.2) then
            write(*,1001)i,j,X
 1001       format('  i = ',i1,', j = ',i1,'; X = ',3(e12.5,2x))           
          endif
        enddo
      enddo
!
!  
      dXdEta = 0.d0
!  ...accumulate      
      do j = 0,(deg-1)
        do i = 0,(deg-1-j)
          call biv_bernstein_poly(i,j,(deg-1),Eta(1),Eta(2), poly,dpoly)
          dXdEta(:,1) =  dXdEta(:,1) + deg*poly*                             &
                        (TRIANGLES(Nt)%Rdata(bijec(i+1,j):bijec(i+1,j)+2) -  &
                         TRIANGLES(Nt)%Rdata(bijec(i,j):bijec(i,j)+2))      
          dXdEta(:,2) =  dxdEta(:,2) + deg*poly*                             &
                        (TRIANGLES(Nt)%Rdata(bijec(i,j+1):bijec(i,j+1)+2) -  &
                         TRIANGLES(Nt)%Rdata(bijec(i,j):bijec(i,j)+2))      
        enddo
      enddo
!
      if (iprint.eq.1) then
        write(*,1002)X
 1002   format(' trian_G1RecTri: X           = ',3(e12.5,2x))
        write(*,1003)dXdEta(1:3,1)
 1003   format('                 dXdEta(:,1) = ',3(e12.5,2x))
        write(*,1004)dXdEta(1:3,2)
 1004   format('                 dXdEta(:,2) = ',3(e12.5,2x))
        call cross_product(dXdEta(1:3,1),dXdEta(1:3,2), nor)
        call normalize(nor)
        write(*,1005)nor
 1005   format('                 normal      = ',3(e12.5,2x))
      endif
!      
!
      end
