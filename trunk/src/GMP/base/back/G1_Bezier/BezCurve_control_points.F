c----------------------------------------------------------------------       
c      
      subroutine BezCurve_control_points(Nc, V1,V2)
c      
c----------------------------------------------------------------------       
      use kinds
      use GMP
c----------------------------------------------------------------------       
      implicit none
c----------------------------------------------------------------------       
      integer,                intent(in)  :: Nc
      real(DP), dimension(3), intent(out) :: V1,V2
c----------------------------------------------------------------------
      real(DP), dimension(3) :: vel1,vel2
      integer :: is,np1,np2       
      integer :: iprint = 0
c----------------------------------------------------------------------       
c
c  ...check curve type
      if (CURVES(abs(Nc))%Type.ne.'cBezier') then
        write(*,*)'BezCurve_control_points: wrong curve type!'
        write(*,*)'        type = ',CURVES(Nc)%Type
        stop
      endif
c      
      if (iprint.eq.1) then
        write(*,7001)Nc
 7001   format(' BezCurve_control_points: Nc = ',i4)
      endif

c  ...allocate Rdata
ccc      allocate(CURVES(Nc)%Rdata(6), STAT=is)
ccc      if (is.ne.0) then
ccc        write(*,*)'BezCurve_control_points: Rdata not allocated'
ccc        write(*,*)'                         for NC = ',Nc
cccc       stop
ccc      endif
c       
c  ...get velocity vectors
      call get_velocities(abs(Nc), vel1,vel2)
c      
      if (iprint.eq.1) then
        write(*,*)'BezCurve_control_points: velocities'
        write(*,7000)vel1,vel2
 7000   format(' vel1 = ',3(f7.3x)'; vel2 = ',3(f7.3x))
      endif
c      
c  ...get curve endpoints
      np1 = CURVES(abs(Nc))%EndPoNo(1)
      np2 = CURVES(abs(Nc))%EndPoNo(2)
      if (iprint.eq.1) then
        write(*,*)'BezCurve_control_points: np1 = ',np1,'np2 = ',np2
      endif 
c  ...compute control points
ccc      CURVES(Nc)%Rdata(1:3) = POINTS(np1)%Rdata(1:3) + vel1/3.d0 
      V1 = POINTS(np1)%Rdata(1:3) + vel1/3.d0 
ccc      CURVES(Nc)%Rdata(4:6) = POINTS(np2)%Rdata(1:3) + vel2/3.d0 
      V2 = POINTS(np2)%Rdata(1:3) + vel2/3.d0 
ccc      V1 = CURVES(Nc)%Rdata(1:3);  V2 = CURVES(Nc)%Rdata(4:6)
c
c
      end
