!----------------------------------------------------------------------------
      subroutine get_curve_velocity(Nc,Np, V)
!----------------------------------------------------------------------------
      use GMP
!----------------------------------------------------------------------------
      implicit none
!----------------------------------------------------------------------------
      integer, intent(in)                :: Nc,Np
      real(8), dimension(3), intent(out) :: V
!----------------------------------------------------------------------------
!
      if (CURVES(Nc)%Type.ne.'5Bezier') then
        write(*,7000)CURVES(Nc)%Type
 7000   format(' get_curve_velocity: invalid curve type ',a10)
        stop
      endif
!
!  ...consistent orientation
      if (CURVES(Nc)%EndPoNo(1).eq.Np) then
        V = 5.d0*(CURVES(Nc)%Rdata(3:5) - CURVES(Nc)%Rdata(0:2))   
!  ...inconsistent orientation        
      elseif (CURVES(Nc)%EndPoNo(2).eq.Np) then
        V = 5.d0*(CURVES(Nc)%Rdata(12:14) - CURVES(Nc)%Rdata(15:17))   
!  ...incompatible data        
      else
        write(*,7001)Nc
 7001   format(' get_curve_velocity: inconsistency, Nc = ',i4)
        stop
      endif        
!
!
      end
