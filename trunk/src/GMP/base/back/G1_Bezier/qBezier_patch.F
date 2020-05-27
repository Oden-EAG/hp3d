c-----------------------------------------------------------------------
c
      subroutine qBezier_patch(Eta, X)
c
c-----------------------------------------------------------------------
c   latest revision    - Feb 10
c
c   purpose            - routine supports parametric transfinite
c                        interpolation for a triangle
c   arguments :
c     in:
c               No     - a GMP triangle number
c               Eta    - reference coordinates of a point
c                        in the triangle
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the parameters
c----------------------------------------------------------------------
c     MODULES
      use kinds
      use GMP
c----------------------------------------------------------------------
      implicit none
c----------------------------------------------------------------------
c     DUMMY ARGUMENTS
      real(DP), dimension(2)   :: Eta
      real(DP), dimension(3)   :: X
c----------------------------------------------------------------------
c     VARIABLES
      real(DP) :: u,v,w
      integer :: nv1,nv2,nv3,ne1,ne2,ne3
      real(DP), dimension(3) :: T1,T2,T3
      real(DP), dimension(3) :: B1,B2
      real(DP), dimension(3) :: B1q,B2q,B3
      real(DP), dimension(3) :: C1,C2
      real(DP), dimension(3) :: C1q,C2q,C3
      real(DP), dimension(3) :: D1,D2
      real(DP), dimension(3) :: D1q,D2q,D3
      real(DP), dimension(3) :: S1,S2,S3,S1q,S2q,S3q
      real(DP), dimension(3) :: P1,P2,P3,P1q,P2q,P3q
      real(DP), dimension(3) :: N12,N13,N23
      real(DP), dimension(3) :: L13,M13,L23,M23,L12,M12
      real(DP), dimension(3) :: K12,K13,K23
      real(DP), dimension(3) :: Z
      real(DP), dimension(3) :: R1,R2,R3
      real(DP), dimension(3) :: a1,a2
      real(DP), dimension(0:4,0:4,0:4,1:3) :: VV
      real(DP), dimension(0:4,0:4,0:4) :: BB
      integer :: deg,i,j,k
      integer :: iprint = 0
c
c----------------------------------------------------------------------
c  STEP 1: determine all control points (BRUTE FORCE APPROACH!)
c
c      if (iprint.eq.1) then
c        write(*,7005) No,Eta
c 7005   format(' trian_G1RecTri: No = ',i4,'; Eta = '2(f7.3x))
c      endif
c
      u = 1.d0 - Eta(1) - Eta(2);  v = Eta(1);  w = Eta(2)
c  .....store control points
        VV(4,0,0,1:3) = (/0, 0, 1/)
        VV(3,1,0,1:3) = (/1, 0, 0/)
        VV(2,2,0,1:3) = (/2, 0, 1/)
        VV(1,3,0,1:3) = (/3, 0, 0/)
        VV(0,4,0,1:3) = (/4, 0, 1/)
c  .............................
        VV(3,0,1,1:3) = (/0, 1, 0/)
        VV(2,1,1,1:3) = (/1, 1, 1/)
        VV(1,2,1,1:3) = (/2, 1, 0/)
        VV(0,3,1,1:3) = (/3, 1, 1/)
c  .............................
        VV(2,0,2,1:3) = (/0, 2, 1/)
        VV(1,1,2,1:3) = (/1, 2, 0/)
        VV(0,2,2,1:3) = (/2, 2, 1/)
c  .............................
        VV(1,0,3,1:3) = (/0, 3, 0/)
        VV(0,1,3,1:3) = (/1, 3, 1/)
c  .............................
        VV(0,0,4,1:3) = (/0, 4, 1/)
c  .............................
c
c
      if (iprint.eq.1) then
        write(*,*)'trian_G1RecTri: subpatch control points'
cc        write(*,7021)ipatch
 7021   format(' ipatch = ',i)
        write(*,7006)VV(4,0,0,1:3)
 7006   format(' V400 = ',3(f7.3x))
        write(*,7007)VV(3,1,0,1:3)
 7007   format(' V310 = ',3(f7.3x))
        write(*,7008)VV(2,2,0,1:3)
 7008   format(' V220 = ',3(f7.3x))
        write(*,7009)VV(1,3,0,1:3)
 7009   format(' V130 = ',3(f7.3x))
        write(*,7010)VV(0,4,0,1:3)
 7010   format(' V040 = ',3(f7.3x))
        write(*,*)'..............................'
        write(*,7011)VV(3,0,1,1:3)
 7011   format(' V301 = ',3(f7.3x))
        write(*,7012)VV(2,1,1,1:3)
 7012   format(' V211 = ',3(f7.3x))
        write(*,7013)VV(1,2,1,1:3)
 7013   format(' V121 = ',3(f7.3x))
        write(*,7014)VV(0,3,1,1:3)
 7014   format(' V031 = ',3(f7.3x))
        write(*,*)'..............................'
        write(*,7015)VV(2,0,2,1:3)
 7015   format(' V202 = ',3(f7.3x))
        write(*,7016)VV(1,1,2,1:3)
 7016   format(' V112 = ',3(f7.3x))
        write(*,7017)VV(0,2,2,1:3)
 7017   format(' V022 = ',3(f7.3x))
        write(*,*)'..............................'
        write(*,7018)VV(1,0,3,1:3)
 7018   format(' V103 = ',3(f7.3x))
        write(*,7019)VV(0,1,3,1:3)
 7019   format(' V013 = ',3(f7.3x))
        write(*,*)'..............................'
        write(*,7020)VV(0,0,4,1:3)
 7020   format(' V004 = ',3(f7.3x))
        write(*,*)'..............................'
      endif
c
c----------------------------------------------------------------------
c  STEP 3: evaluate quartic Bezier patch
c
c  ...determine bivariate Bernstein polynomials
      call biv_Bernstein_poly(u,v,w, BB)
      k = 0
      X = 0.d0
      deg = 4
      do i = 0,deg
        do j = 0,deg-i
          k = k+1
          X = X + BB(i,j,deg-i-j)*VV(i,j,deg-i-j,1:3)
        enddo
      enddo
c
c     dirty fix:
c
c
      end
