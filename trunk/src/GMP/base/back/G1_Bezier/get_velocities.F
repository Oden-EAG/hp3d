c----------------------------------------------------------------------
c      
      subroutine get_velocities(Nc, Vel1,Vel2)
c      
c----------------------------------------------------------------------
      use kinds
      use GMP
c----------------------------------------------------------------------
      implicit none
c----------------------------------------------------------------------
      integer,                intent(in)  :: Nc
      real(DP), dimension(3), intent(out) :: Vel1,Vel2
c----------------------------------------------------------------------
      real(DP), dimension(3) :: FF = (/ 1.d0, 0.d0, 0.d0/)
      real(DP), dimension(3) :: BB = (/-1.d0, 0.d0, 0.d0/)
      real(DP), dimension(3) :: RR = (/ 0.d0, 1.d0, 0.d0/)
      real(DP), dimension(3) :: LL = (/ 0.d0,-1.d0, 0.d0/)
      real(DP), dimension(3) :: UU = (/ 0.d0, 0.d0, 1.d0/)
      real(DP), dimension(3) :: DD = (/ 0.d0, 0.d0,-1.d0/)
      real(DP) :: alpha = 1.3d0
      integer :: iprint = 0
c----------------------------------------------------------------------      
c      
      call compute_tangents(Nc, Vel1,Vel2)
      Vel1 =  alpha*Vel1
c  ...flip sign of 2nd tangent       
      Vel2 = -alpha*Vel2

c      select case(Nc)
c        case(2)
c          Vel1 = alpha*RR;  Vel2 = alpha*FF
c        case(7)
c          Vel1 = alpha*BB;  Vel2 = alpha*RR
c        case(10)
c          Vel1 = alpha*LL;  Vel2 = alpha*BB
c        case(13)
c          Vel1 = alpha*FF;  Vel2 = alpha*LL
c        case(5)
c          Vel1 = alpha*UU;  Vel2 = alpha*FF
c        case(6)
c         Vel1 = alpha*UU;  Vel2 = alpha*RR
c        case(9)
c          Vel1 = alpha*UU;  Vel2 = alpha*BB
c        case(12)
c          Vel1 = alpha*UU;  Vel2 = alpha*LL
c        case(15)
c          Vel1 = alpha*DD;  Vel2 = alpha*FF
c        case(16)
c          Vel1 = alpha*DD;  Vel2 = alpha*RR
c        case(17)
c          Vel1 = alpha*DD;  Vel2 = alpha*BB
c        case(18)
c          Vel1 = alpha*DD;  Vel2 = alpha*LL
c        case default
c          write(*,*)'get_velocities: not a Bezier curve!'
c          stop
c      endselect
c
      if (iprint.eq.1) then
        write(*,7000)Nc
 7000   format(' get_velocities: Nc = ',i4)       
        write(*,7001)Vel1,Vel2
 7001   format(' Vel1 = ',3(f7.3x),'; Vel2 = ',3(f7.3x))
      endif      
c
c
      end




c----------------------------------------------------------------------      
c      
      subroutine compute_tangents(Nc, T1, T2)
c      
c----------------------------------------------------------------------      
      use kinds
      use GMP
c----------------------------------------------------------------------      
      implicit none
c----------------------------------------------------------------------      
      integer, intent(in) :: Nc
      real(DP), dimension(3),intent(out) :: T1,T2
c----------------------------------------------------------------------      
      real(DP), dimension(3) :: N
      real(DP) :: prod
      integer :: np1,np2
c----------------------------------------------------------------------      
c
c
c  ...get endpoints       
      np1 = CURVES(abs(Nc))%EndPoNo(1)
      np2 = CURVES(abs(Nc))%EndPoNo(2)
c  ...check endpoints type      
      if (POINTS(np1)%Type.ne.'CoorNrm') then
        write(*,*)'invalid point type, np1 = ',np1
        call print_GMP
        stop
      endif
      if (POINTS(np2)%Type.ne.'CoorNrm') then
        write(*,*)'invalid point type, np2 = ',np2
        call print_GMP
        stop
      endif
c  ...store 1st endpoint normal      
      N = POINTS(np1)%Rdata(4:6)      
      T1 = POINTS(np2)%Rdata(1:3) - POINTS(np1)%Rdata(1:3); T2 = T1
      call scalar_product(T1,N, prod)
c  ...compute tanget to 1st endpoint      
      T1 = T1 - prod*N
c  ...store 2nd endpoint normal     
      N = POINTS(np2)%Rdata(4:6)      
      call scalar_product(T2,N, prod)
      T2 = T2 - prod*N
c
c
      end
