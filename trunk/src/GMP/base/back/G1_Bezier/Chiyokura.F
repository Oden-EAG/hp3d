c----------------------------------------------------------------------
c
      subroutine Chiyokura(No,Ie,a0,a3, a1,a2)
c
c----------------------------------------------------------------------
      use kinds
      use GMP
      use control
c----------------------------------------------------------------------
      implicit none
c----------------------------------------------------------------------
c     DUMMY ARGUMENTS
      integer, intent(in)                :: No         ! triangle #
      integer, intent(in)                :: Ie         ! local edge #
      real(DP), dimension(3),intent(in)  :: a0,a3      ! first 2 control
c                             points on other two edges of the triangle
      real(DP), dimension(3),intent(out) :: a1,a2      ! inner control
c                                                       points
c----------------------------------------------------------------------
c     LOCAL VARIABLES
      integer :: nc,np1,np2,nc1,nc2,iflag,ntf,nick,nt,lab,nnt,iie,i
      integer, dimension(3) :: nnc
      real(DP),dimension(3) :: V0,V1,V2,V3,s0,s1,s2,b0,b1,b3,void
      real(DP) :: h0,h1,k0,k1,det,prod
      real(DP), dimension(3)   :: temp
      integer :: iprint = 0
c----------------------------------------------------------------------
      real(DP), dimension(3,2) :: AA
      real(DP), dimension(3)   :: BB
      real(DP), dimension(4)   :: WORK
      integer                  :: INFO
c----------------------------------------------------------------------
c  STEP 1: get curve control points and compute tangent vectors
c
      if (iprint.eq.1) then
        write(*,7003)No,Ie
 7003   format(' Chiyokura: No = ',i5,'; Ie = ',i2)
        write(*,7011)a0,a3
 7011   format(' a0 = ',3(f7.3x),'; a3 = ',3(f7.3x))
      endif
c
c  ...get curve number
      nc = TRIANGLES(No)%EdgeNo(Ie)
      if (iprint.eq.1) then
        write(*,7007)nc
 7007   format(' Chiyokura: nc = ',i5)
      endif
c
c  ...get control points
      if (nc.gt.0) then
c  .....if positive orientation
        call BezCurve_control_points(nc, V1,V2)
        np1 = CURVES(nc)%EndPoNo(1)
        V0  = POINTS(np1)%Rdata(1:3)
        np2 = CURVES(nc)%EndPoNo(2)
        V3  = POINTS(np2)%Rdata(1:3)
      elseif (nc.lt.0) then
c  .....if negative orientation
        call BezCurve_control_points(-nc, V2,V1)
        np1 = CURVES(nc)%EndPoNo(2)
        V0  = POINTS(np1)%Rdata(1:3)
        np2 = CURVES(nc)%EndPoNo(1)
        V3  = POINTS(np2)%Rdata(1:3)
cccc        void = a0;  a0 = a3;  a3 = void
      else
        write(*,*)'Chiyokura: invalid curve number, nc = ',nc
        stop
      endif
      if (iprint.eq.1) then
        write(*,*)'Chiyokura: curve control points determined.'
        write(*,7000) V0,V1
 7000   format(' V0 = ',3(f7.3,2x),'; V1 = ',3(f7.3,2x))
        write(*,7002) V2,V3
 7002   format(' V2 = ',3(f7.3,2x),'; V3 = ',3(f7.3,2x))
      endif
c  ...disregard orientation
      nc = abs(nc)
c
c  ...compute tangent vectors along curve
      s0 = V1 - V0;  s1 = V2 - V1; s2 = V3 - V2
      if (iprint.eq.1) then
        write(*,*)'Chiyokura: tangent vectors determined.'
        write(*,7001)s0,s1,s2
 7001   format(' s0 = ',3(f7.3,2x),'; s1 = ',3(f7.3,2x),';
     .s3 = ',3(f7.3,2x))
      endif
c
c----------------------------------------------------------------------
c  STEP 2: determine two adjacent curves
c
c  ...get number of figures connected to the curve
      ntf = CURVES(nc)%NrFig
c  ...loop over connected triangles
      iflag = 0
      do i = 1, ntf
        nick = CURVES(nc)%FigNo(i)
        call decode(nick, nt,lab)
c  .....cycle if original triangle
        if (nt.eq.No) cycle
c  .....if on the G1 surface store number and increment counter
        if (TRIANGLES(nt)%Type.eq.'G1RecTri') then
          nnt = nt
          iflag = iflag + 1
        endif
c  ...end loop over connected triangles
      enddo
c
c  ...check found triangle
      if (iflag.eq.0) then
        write(*,*)'Chiyokura: no adjacent triangle has been found!'
        stop
      endif
      if (iflag.gt.1) then
        write(*,*)'Chiyokura: more then 1 adjacent triangle was found!'
        stop
      endif
      if (iprint.eq.1) then
        write(*,*)'Chiyokura: neighbooring triangle = ',nnt
      endif
c
c  ...loop over neighboor's edges
      iflag = 0
      do iie=1,3
c  .....if neighboor's edge coincides with original edge, cycle
        if (abs(TRIANGLES(nnt)%EdgeNo(iie)).eq.nc) cycle
        iflag = iflag + 1
        nnc(iflag) = abs(TRIANGLES(nnt)%EdgeNo(iie))
      enddo
      if (iprint.eq.1) then
        write(*,7008)nnc(1),nnc(2)
 7008   format(' Chihokura: nnc(1) = ',i4,'; nnc(2) = ',i4)
      endif
c
c  ...check
      if (iflag.eq.3) then
        write(*,*)'Chiyokura: more than 2 adjacent edges found!'
        stop
      endif
c
c  ...loop over curves attached to np1
      if (iprint.eq.1) then
        write(*,7009)POINTS(np1)%NrCurv
 7009   format(' Chiyokura: # curves for np1 = ',i3)
      endif
      iflag = 0
      do i=1,POINTS(np1)%NrCurv
        if (iprint.eq.1) then
          write(*,7004)i,POINTS(np1)%CurvNo(i)
 7004     format(' i = ',i2x,'; nc = ',i4)
        endif
c  .....check against 1st neighboring curve
        if (POINTS(np1)%CurvNo(i).eq.nnc(1)) then
          nc1 = nnc(1)
          iflag = 1
          exit
        endif
c  ...end of loop over curves attached to np1
      enddo
c
      if (iflag.eq.1) then
        nc2 = nnc(2)
      else
        nc1 = nnc(2)
        nc2 = nnc(1)
      endif
      if (iprint.eq.1) then
         write(*,7010)nc1,nc2
 7010    format(' Chiyokura: nc1 = ',i4'; nc2 = ',i4)
      endif
c
c----------------------------------------------------------------------
c  STEP 3: determine control points of adjacent curves and tangent
c          vectors
c
c  ...curve adjacent to np1
      if (CURVES(nc1)%EndPoNo(1).eq.np1) then
        call BezCurve_control_points(nc1, b0,void)
      else
        call BezCurve_control_points(nc1, void,b0)
      endif
c
c  ...curve adjacent to np2
      if (CURVES(nc2)%EndPoNo(1).eq.np2) then
        call BezCurve_control_points(nc2, b3,void)
      else
        call BezCurve_control_points(nc2, void,b3)
      endif
c
c  ...compute b0, b3
      b0 = b0 - V0;  b3 = b3 - V3
      if (iprint.eq.1) then
        write(*,7005)b0
 7005   format(' Chiyokura: b0 = '3(f7.3x))
        write(*,7006)b3
 7006   format(' Chiyokura: b3 = '3(f7.3x))
      endif
c
      b0 = b0 - a0;  call normalize(b0)
      b3 = b3 - a3;  call normalize(b3)
      b1 = (b0 + b3)/2.d0
c
c  ...check linear dependecy of a0,b0,s0
      call mixed_product(a0,b0,s0, prod)
      if (abs(prod).gt.GEOM_TOL) then
        write(*,*)'Chiyokura: a0,b0,s0 are L.I., prod = ',prod
        stop
      endif
c  ...solve least square system
      AA(:,1) = s0;  AA(:,2) = b0;  BB = a0
      call DGELS('N',3,2,1,AA,3,BB,3,WORK,4,INFO)
      if (INFO.ne.0) then
        write(*,*)'Chiyokura: linear system not solved (1)'
        stop
      endif
      h0 = BB(1);  k0 = BB(2)
c  ...check solution
      temp = h0*s0 + k0*b0 - a0
      call norm(temp, temp(1))
      if (temp(1).gt.GEOM_TOL) then
        write(*,*)'Chiyokura: geom toll exceeded, temp = ',temp(1)
        stop
      endif
      if (iprint.eq.1) then
        write(*,7013)h0,k0
 7013   format(' Chiyokura: h0 = ',f7.3x,'; k0 = ',f7.3x)
      endif
c


ccc  ...use Cramer rule to solve eq's (6.10)
cc      det = s0(1)*b0(2) - s0(2)*b0(1)
cc      if (abs(det).lt.GEOM_TOL) then
cc        write(*,*)'Chiyokura: det = 0 (1)!'
cc        stop
cc      endif
cc      h0 = (a0(1)*b0(2) - b0(1)*a0(2))/det
cc      k0 = (a0(2)*s0(1) - a0(1)*s0(2))/det
c  ...check compatibility
ccc      if (abs(s0(3)*h0+b0(3)*k0-a0(3)).gt.GEOM_TOL) then
ccc        write(*,*)'Chiyokura: incompatible system (1)!'
ccc        write(*,*)'3rd eq = ',abs(s0(3)*h0+b0(3)*k0-a0(3))
ccc        stop
ccc      endif


c
c  ...check linear dependecy of a3,b3,s2
      call mixed_product(a3,b3,s2, prod)
      if (abs(prod).gt.GEOM_TOL) then
        write(*,*)'Chiyokura: a3,b3,s2 are L.I., prod = ',prod
        stop
      endif
c  ...solve least square system
      AA(:,1) = s2;  AA(:,2) = b3;  BB = a3
      call DGELS('N',3,2,1,AA,3,BB,3,WORK,4,INFO)
      if (INFO.ne.0) then
        write(*,*)'Chiyokura: linear system not solved (2)'
        stop
      endif
      h1 = BB(1);  k1 = BB(2)
c  ...check solution
      temp = h1*s2 + k1*b3 - a3
      call norm(temp, temp(1))
      if (temp(1).gt.GEOM_TOL) then
        write(*,*)'Chiyokura: geom toll exceeded, temp = ',temp(1)
        stop
      endif
      if (iprint.eq.1) then
        write(*,7014)h1,k1
 7014   format(' Chiyokura: h1 = ',f7.3x,'; k1 = ',f7.3x)
      endif


c
ccc  ...use Cramer rule to solve eq's (6.10)
cc      det = s2(1)*b3(2) - s2(2)*b3(1)
cc      if (abs(det).lt.GEOM_TOL) then
cc        write(*,*)'Chiyokura: det = 0 (2)!'
cc        stop
cc      endif
cc      h1 = (a3(1)*b3(2) - b3(1)*a3(2))/det
cc      k1 = (a3(2)*s2(1) - a3(1)*s2(2))/det
c  ...check compatibility
cc      if (abs(s2(3)*h1+b3(3)*k1-a3(3)).gt.GEOM_TOL) then
cc        write(*,*)'Chiyokura: incompatible system (2)!'
cc        stop
cc      endif
c




c  ...use eq's (6.11) to compute a1, a2
      a1 = (2.d0*k0*b1 + k1*b0 + 2.d0*h0*s1 + h1*s0)/3.d0
      a2 = (k0*b3 + 2.d0*k1*b1 + h0*s2 + 2.d0*h1*s1)/3.d0
      if (iprint.eq.1) then
        write(*,*)'Chiyokura: have determined control points'
        write(*,7012)a1,a2
 7012   format(' a1 = ',3(f7.3x),'; a2 = ',3(f7.3x))
        call pause
      endif
c
c
      end
