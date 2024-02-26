#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - initwin
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine initiates a specific display
!
!   argument
!           in : Iwind - type of display
!
!   required routines  - openwind, closwind, clrwind
!
!----------------------------------------------------------------------
      subroutine initwin(Iwind)
!
      use graphmod
!
      implicit none
!
      integer, intent(in) :: Iwind
!
      integer :: iwindtype
!
#if HP3D_DEBUG
!  ...debug print
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'initwin: Iwind = ',Iwind
        call pause
      endif
#endif
!
!  ...store the display type
      IDISPLAY_TYPE = Iwind
!
      IWINDNUM = 1
      iwindtype = 0
!
!!!      IWINDL = 750
!!!      IWINDH = 600
!!!      IWINDL = 667
!!!      IWINDH = 400
!      IWINDL = 1334
!      IWINDH = 800
      RWINDL = IWINDL
      RWINDH = IWINDH
      ISIZE(1) = 0
      ISIZE(2) = 0
      ISIZE(3) = IWINDL
      ISIZE(4) = IWINDH
      call openwind(iwindtype,ISIZE,IWINDNUM)
!
!  ...margin
      RMARGIN = .05d0*RWINDH
!
!  ...size of the window after subtracting margin...
      XLENGTH = RWINDL - 3.d0*RMARGIN
      YLENGTH = RWINDH - 2.d0*RMARGIN
!
!  ...choose colormap
      call colorset(Iwind)
!
      call clrwind(NPCOL(1))
!
      write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
      write(*,*) 'WINDOW TO CONTINUE ...'
!
      call closwind(IWINDNUM)
!
      end subroutine initwin

#endif
