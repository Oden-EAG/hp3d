subroutine add_void_curve(Iv1, Iv2, Icurve)
  use GMP
  implicit none
  integer, intent(in)  :: Iv1, Iv2
  integer, intent(out) :: Icurve
  integer :: iflag, ic
  iflag = 0
  do ic=1, NRCURVE
     if ((Iv1.eq.CURVES(ic)%EndPoNo(1)).and.(Iv2.eq.CURVES(ic)%EndPoNo(2))) then
        iflag = 1
        exit
     end if
     if ((Iv2.eq.CURVES(ic)%EndPoNo(1)).and.(Iv1.eq.CURVES(ic)%EndPoNo(2))) then
        iflag = 1
        exit
     end if
  enddo
  if (iflag.eq.0) then
     NRCURVE = NRCURVE + 1 
     Icurve  = NRCURVE; 
     CURVES(Icurve)%Type = 'Seglin'
     CURVES(Icurve)%EndPoNo(1) = Iv1
     CURVES(Icurve)%EndPoNo(2) = Iv2
     CURVES(Icurve)%NrFig = 0
  else
     Icurve = ic
  end if
end subroutine add_void_curve
!
