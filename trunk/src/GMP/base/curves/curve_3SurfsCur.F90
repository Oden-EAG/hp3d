!----------------------------------------------------------------------------
!> @brief parametrization of a curve lying on intersection of 2 algebraic
!!           surfaces
!!
!> @param[in ] No     - the curve number
!> @param[in ] Eta    - reference coordinate  (between 0 and 1)
!> @param[out] X      - physical coordinates of the point
!> @param[out] Dxdeta - derivatives of the physical coordinates
!!
!> @date Nov 12
!----------------------------------------------------------------------------
subroutine curve_3SurfsCur(Nc,Eta, R,Dr)
!
      implicit none
      integer             ,intent(in ) :: Nc
      real(8)             ,intent(in ) :: Eta
      real(8),dimension(3),intent(out) :: R
      real(8),dimension(3),intent(out) :: Dr
!
      integer :: iflag
!
!----------------------------------------------------------------------------
!
!  ...try all possible choices of 2 out of 3 surfaces
      call curve_2SurfsCur(Nc,1,2,Eta, R,Dr,iflag)
!
      if (iflag.eq.1) then
      call curve_2SurfsCur(Nc,1,3,Eta, R,Dr,iflag)
      endif
!
      if (iflag.eq.1) then
      call curve_2SurfsCur(Nc,2,3,Eta, R,Dr,iflag)
      endif
!
      if (iflag.eq.1) then
        write(*,1000)Nc,Eta
 1000   format(' curve_3SurfsCur: Nc,Eta = ',i7,2x,e12.5)
        stop
      endif
!
!
end subroutine curve_3SurfsCur
