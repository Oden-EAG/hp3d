c---------------------------------------------------------------------
c      
      subroutine curve_cBezier(Nc,Eta, X,Dxdeta)
c
c----------------------------------------------------------------------
      use GMP
c----------------------------------------------------------------------
      implicit none
c----------------------------------------------------------------------
      integer, intent(in) :: Nc
      real(8), intent(in) :: Eta
      real(8), dimension(3), intent(out) :: X
      real(8), dimension(3), intent(out) :: Dxdeta
c----------------------------------------------------------------------
      integer :: i,np1,np2
      real(8), dimension(0:3)     :: B
      real(8), dimension(0:3,1:3) :: V
c----------------------------------------------------------------------
c
c  ...get curve's endpoints      
      np1 = CURVES(abs(Nc))%EndPoNo(1)
      np2 = CURVES(abs(Nc))%EndPoNo(2)
c
c  ...compute control points      
      V(0,1:3) = POINTS(np1)%Rdata(1:3)
      V(3,1:3) = POINTS(np2)%Rdata(1:3)
      call BezCurve_control_points(abs(Nc), V(1,1:3),V(2,1:3))
c
c  ...compute Bernstein polymonials      
      call bernstein_poly(Eta, B)
c
c  ...compute curve's parameterization
      X = 0.d0
      do i=0,3
        X = X + V(i,1:3)*B(i)
      enddo
c
c  ...dirty trick
      Dxdeta = 0.d0
c
c      
      end
