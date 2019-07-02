c-----------------------------------------------------------------------
c      
      subroutine trian_G1RecTri(No,Eta, X,Dxdeta)
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
      integer                  :: No
      real(DP), dimension(2)   :: Eta
      real(DP), dimension(3)   :: X
      real(DP), dimension(3,2) :: Dxdeta
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
      integer :: deg,i,j,ipatch
      integer :: iprint = 0
c
c----------------------------------------------------------------------      
c  STEP 1: determine all control points (BRUTE FORCE APPROACH!)
c
      if (iprint.eq.1) then
        write(*,7005) No,Eta
 7005   format(' trian_G1RecTri: No = ',i4,'; Eta = '2(f7.3x))
      endif
c      
c  ...get triangle's vertices
      nv1 = TRIANGLES(No)%VertNo(1);  T1 = POINTS(nv1)%Rdata(1:3)
      nv2 = TRIANGLES(No)%VertNo(2);  T3 = POINTS(nv2)%Rdata(1:3)
      nv3 = TRIANGLES(No)%VertNo(3);  T2 = POINTS(nv3)%Rdata(1:3)
      if (iprint.eq.1) then
        write(*,*)'trian_G1RecTri: vertices have been determined'
        write(*,7000)T1,T3,T2
 7000   format(' T1 = ',3(f7.3x)'; T3 = ',3(f7.3x),'; T2 = ',3(f7.3x))
      endif
c
c  ...get control points of edges: B1,B2; C1,C2; D1,D2;
      ne1 = TRIANGLES(No)%EdgeNo(1)
      if (ne1.gt.0) then
        call BezCurve_control_points(ne1, B1,B2)
      else
        call BezCurve_control_points(-ne1, B2,B1)
      endif
      ne2 = TRIANGLES(No)%EdgeNo(2)
      if (ne2.gt.0) then
        call BezCurve_control_points(ne2,D1,D2) 
      else
        call BezCurve_control_points(-ne2, D2,D1)
      endif
      ne3 = TRIANGLES(No)%EdgeNo(3)
      if (ne3.gt.0) then
        call BezCurve_control_points(ne3, C1,C2)
      else
        call BezCurve_control_points(-ne3, C2,C1)
      endif
      if (iprint.eq.1) then
        write(*,*)'trian_G1RecTri: edges control points determined'
        write(*,7001)B1,B2
 7001   format(' B1 = ',3(f7.3x),'; B2 = ',3(f7.3x))
        write(*,7002)D1,D2
 7002   format(' D1 = ',3(f7.3x),'; D2 = ',3(f7.3x))
        write(*,7003)C1,C2
 7003   format(' C1 = ',3(f7.3x),'; C2 = ',3(f7.3x))
      endif
c      
c  ...compute S1,S2,S3 as center of masses of suitable triangles
      S1 = (T1 + B1 + C1)/3.d0
      S2 = (T2 + C2 + D2)/3.d0
      S3 = (T3 + B2 + D1)/3.d0
      if (iprint.eq.1) then
        write(*,*)'trian_G1RecTri: have determined S1,S2,S3'
        write(*,7004)S1,S2,S3
 7004 format(' S1 = ',3(f7.3x),'; S2 =  ',3(f7.3x),'; S3 = ',3(f7.3x))
      endif
c
c  ...elevate degree to obtain corresponding control points for a 
c     quartic: B1q,B2q,B3; C1q,C2q,C3; D1q,D2q,D3;
      call cub2quar_BezCurve(T1,B1,B2,T3, B1q,B3,B2q)
      call cub2quar_BezCurve(T3,D1,D2,T2, D1q,D3,D2q)
      call cub2quar_BezCurve(T1,C1,C2,T2, C1q,C3,C2q)
c
c  ...compute 1st row of inner control points using Chiyokura's 
c     method: L13,M13; L23,M23; L12,M12
      call Chiyokura(No,1,S1-T1,S3-T3, a1,a2)
      L13 = B1 + a1;  M13 = B2 + a2
      call Chiyokura(No,2,S3-T3,S2-T2, a1,a2)
      M23 = D1 + a1;  K23 = D2 + a2
      call Chiyokura(No,3,S1-T1,S2-T2, a1,a2)
      L12 = C1 + a1;  K12 = C2 + a2
c
c  ...compute P1,P2,P3 using eq's (7.12)
      P1 = (T1 - 3.d0*S1 + 4.d0*L13 + 4.d0*L12)/6.d0
      P2 = (T2 - 3.d0*S2 + 4.d0*K23 + 4.d0*K12)/6.d0
      P3 = (T3 - 3.d0*S3 + 4.d0*M23 + 4.d0*M13)/6.d0
c
c  ...compute Z as the center of mass of P1,P2,P3
      Z = (P1 + P2 + P3)/3.d0
c
c  ...compute remaining inner control points N13,N12,N23 using eq's (7.13)
      N12 = (-S1 - S2 + S3 + 4.d0*P1 + 4.d0*P2 - 3.d0*P3)/4.d0
      N13 = (-S1 + S2 - S3 + 4.d0*P1 - 3.d0*P2 + 4.d0*P3)/4.d0
      N23 = ( S1 - S2 - S3 - 3.d0*P1 + 4.d0*P2 + 4.d0*P3)/4.d0
c
c  ...upgrade control points P's and S's of inner boundary to quartic 
c     control points
      call cub2quar_BezCurve(T1,S1,P1,Z, S1q,R1,P1q)
      call cub2quar_BezCurve(T3,S3,P3,Z, S3q,R3,P3q)
      call cub2quar_BezCurve(T2,S2,P2,Z, S2q,R2,P2q)
c
c----------------------------------------------------------------------      
c  STEP 2: determine subpatch (i-th subpatch is adjacent to i-th local
c           edge) and store appropriate control points
c
c  ...baricentrical coordinates
      u = 1.d0 - Eta(1) - Eta(2);  v = Eta(1);  w = Eta(2)
c
c  ...1st subpatch.....................................................
      if ((Eta(1)-Eta(2).ge.0.d0).and.
     .                    (1.d0-Eta(1)-2.d0*Eta(2).ge.0.d0)) then
        ipatch = 1
        if (iprint.eq.1) then
          write(*,*)'trian_G1RecTri: point belongs to 1st subpatch'
        endif
c  .....compute baricentrical coordinates wrt subpatch
        u = 1.d0 - Eta(1) - 2.d0*Eta(2)
        v = Eta(1) - Eta(2)
        w = 3.d0*Eta(2)
        if (iprint.eq.1) then
          write(*,*)'trian_G1RecTri: baricentrical coord wrt subpatch'
          write(*,7022)u,v,w
 7022     format(' u = ',f7.3x,'; v = ',f7.3x,'; w = ',f7.3)
        endif
c  .....store control points
        VV(4,0,0,1:3) = T1
        VV(3,1,0,1:3) = B1q
        VV(2,2,0,1:3) = B3
        VV(1,3,0,1:3) = B2q
        VV(0,4,0,1:3) = T3
c  ............................. 
        VV(3,0,1,1:3) = S1q
        VV(2,1,1,1:3) = L13
        VV(1,2,1,1:3) = M13
        VV(0,3,1,1:3) = S3q
c  ............................. 
        VV(2,0,2,1:3) = R1
        VV(1,1,2,1:3) = N13
        VV(0,2,2,1:3) = R3
c  ............................. 
        VV(1,0,3,1:3) = P1q
        VV(0,1,3,1:3) = P3q
c  ............................. 
        VV(0,0,4,1:3) = Z
c  ............................. 
c        
c  ...2nd subpatch.....................................................
      elseif (2.d0*Eta(1)+Eta(2)-1.d0.ge.0.d0) then
        ipatch = 2
        if (iprint.eq.1) then
          write(*,*)'trian_G1RecTri: point belongs to 2nd subpatch'
        endif
c  .....compute baricentrical coordinates wrt subpatch
        u = 2.d0*Eta(1) + Eta(2) - 1.d0
        v = Eta(1) + 2.d0*Eta(2) - 1.d0
        w = 3.d0 - 3.d0*Eta(1) - 3.d0*Eta(2)
        if (iprint.eq.1) then
          write(*,*)'trian_G1RecTri: baricentrical coord wrt subpatch'
          write(*,7022)u,v,w
        endif
c  .....store control points
        VV(4,0,0,1:3) = T3
        VV(3,1,0,1:3) = D1q
        VV(2,2,0,1:3) = D3
        VV(1,3,0,1:3) = D2q
        VV(0,4,0,1:3) = T2
c  ............................. 
        VV(3,0,1,1:3) = S3q
        VV(2,1,1,1:3) = M23
        VV(1,2,1,1:3) = K23
        VV(0,3,1,1:3) = S2q
c  ............................. 
        VV(2,0,2,1:3) = R3
        VV(1,1,2,1:3) = N23
        VV(0,2,2,1:3) = R2
c  ............................. 
        VV(1,0,3,1:3) = P3q
        VV(0,1,3,1:3) = P2q
c  ............................. 
        VV(0,0,4,1:3) = Z
c  ............................. 
c
c  ...3rd subpatch.....................................................
      else
        ipatch = 3
        if (iprint.eq.1) then
          write(*,*)'trian_G1RecTri: point belongs to 3rd subpatch'
        endif
c  .....compute baricentrical coordinates wrt subpatch
        u = Eta(2) - Eta(1)
        v = 1.d0 - 2.d0*Eta(1) - Eta(2)
        w = 3.d0*Eta(1)
        if (iprint.eq.1) then
          write(*,*)'trian_G1RecTri: baricentrical coord wrt subpatch'
          write(*,7022)u,v,w
        endif
c  .....store control points
        VV(4,0,0,1:3) = T2
        VV(3,1,0,1:3) = C2q
        VV(2,2,0,1:3) = C3
        VV(1,3,0,1:3) = C1q
        VV(0,4,0,1:3) = T1
c  ............................. 
        VV(3,0,1,1:3) = S2q
        VV(2,1,1,1:3) = K12
        VV(1,2,1,1:3) = L12
        VV(0,3,1,1:3) = S1q
c  ............................. 
        VV(2,0,2,1:3) = R2
        VV(1,1,2,1:3) = N12
        VV(0,2,2,1:3) = R1
c  ............................. 
        VV(1,0,3,1:3) = P2q
        VV(0,1,3,1:3) = P1q
c  ............................. 
        VV(0,0,4,1:3) = Z
c  ............................. 
c
      endif
c
      if (iprint.eq.1) then
        write(*,*)'trian_G1RecTri: subpatch control points'      
        write(*,7021)ipatch
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

      X = 0.d0
      deg = 4
      do i = 0,deg
        do j = 0,deg-i
          X = X + BB(i,j,deg-i-j)*VV(i,j,deg-i-j,1:3)
        enddo
      enddo
c
c     dirty fix:
      Dxdeta = 0.d0
c      
ccc      call pause
c
      end
