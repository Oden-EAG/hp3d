c----------------------------------------------------------------------
c
c   routine name       - arctan
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine computes the angle whose tangent
c                        is Yc/Xc
c
c   arguments
c           in : Xc,Yc - usually the cartesian coordinates
c                        of given point
c           out: Ang   - the angle between radius of given
c                        point and x-axis ( 0 < Ang < 2*Pi )
c
c----------------------------------------------------------------------
      subroutine arctan(Xc,Yc,Ang)
c
      implicit none
c
      real(8), intent(in)  :: Xc,Yc
      real(8), intent(out) :: Ang
c
      real(8), parameter :: small = 1.d-10
      real(8), parameter :: pi    = 3.14159265358979d0
c
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
c
      end subroutine arctan
