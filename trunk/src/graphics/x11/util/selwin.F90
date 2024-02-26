#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   latest revision     - Feb 2024
!
!   purpose             - routine selects specific display
!
!   argument in : Iwind - type of display
!
!   required routines   - openwind, closwind, clrwind
!
!----------------------------------------------------------------------
!
      subroutine selwin(Iwind)
!
      use graphmod
!
      implicit none
!
      integer, intent(in) :: Iwind
!
      integer :: iwindtype = 0
!
      IWINDNUM = 1
!
!  ...these parameters are hardwired in the interface with X-windows;
!     if you want to change them, you will have to change them
!     in the interface as well !!
!!!      iwindl = 750
!!!      iwindl = 1000
!!!      iwindh = 600
!!!      iwindl = 667
!!!      iwindh = 400
!      iwindl = 1334
!      iwindh = 800
      rwindl = iwindl
      rwindh = iwindh
      ISIZE(1) = 0
      ISIZE(2) = 0
      ISIZE(3) = iwindl
      ISIZE(4) = iwindh
      call openwind(iwindtype,ISIZE,IWINDNUM)
      call setpost(1)
!
!  ...margin
      rmargin = .05d0*rwindh
!
!  ...size of the window after subtracting margin...
      xlength = rwindl - 3.d0*rmargin
      ylength = rwindh - 2.d0*rmargin
!
      call clrwind(npcol(1))
!
!
      end subroutine selwin

#endif
