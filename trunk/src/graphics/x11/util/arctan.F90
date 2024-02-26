#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - arctan
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine computes the angle whose tangent
!                        is Yc/Xc
!
!   arguments
!           in : Xc,Yc - usually the cartesian coordinates
!                        of given point
!           out: Ang   - the angle between radius of given
!                        point and x-axis ( 0 < Ang < 2*Pi )
!
!----------------------------------------------------------------------
      subroutine arctan(Xc,Yc,Ang)
!
      implicit none
!
      real(8), intent(in)  :: Xc,Yc
      real(8), intent(out) :: Ang
!
      real(8), parameter :: small = 1.d-10
      real(8), parameter :: pi    = 3.14159265358979d0
!
      if (dabs(Xc).le.small) then
         if (dabs(Yc).le.small) then
            Ang = 0.0d0
         else
            if (Yc.ge.0.0d0) then
               Ang = 0.5d0*pi
            else
               Ang = 1.5d0*pi
            endif
         endif
      else
         if (Xc.le.0.0d0) then
            Ang=atan(Yc/Xc)+pi
         else
            Ang=atan(Yc/Xc)
         endif
      endif
!
      end subroutine arctan

#endif
