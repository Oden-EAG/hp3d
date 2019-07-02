c----------------------------------------------------------------------
c
c   routine name       - initwin
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Feb 98
c
c   purpose            - routine initiates a specific display
c
c   usage              - call initwin(Iwind)
c
c   argument 
c           in : Iwind - type of display
c
c   required  routines  - openwind, closwind, clrwind
c
c----------------------------------------------------------------------
c
      subroutine initwin(Iwind)
      use graphmod
c     
#include "syscom.blk"
c
c  ...debug print
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'initwin: Iwind = ',Iwind
        call pause
      endif
c
c  ...store the display type
      IDISPLAY_TYPE = Iwind
c
      IWINDNUM = 1
      iwindtype = 0
c
ccc      IWINDL = 750
ccc      IWINDH = 600
ccc      IWINDL = 667
ccc      IWINDH = 400
c      IWINDL = 1334
c      IWINDH = 800
      RWINDL = IWINDL
      RWINDH = IWINDH
      ISIZE(1) = 0
      ISIZE(2) = 0
      ISIZE(3) = IWINDL
      ISIZE(4) = IWINDH
      call openwind(iwindtype,ISIZE,IWINDNUM)
c
c  ...margin
      RMARGIN = .05d0*RWINDH
c
c  ...size of the window after subtracting margin...
      XLENGTH = RWINDL - 3.d0*RMARGIN
      YLENGTH = RWINDH - 2.d0*RMARGIN
c
c  ...choose colormap
      call colorset(Iwind)
c
      call clrwind(NPCOL(1))
c
      write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
      write(*,*) 'WINDOW TO CONTINUE ...'
c
      call closwind(IWINDNUM)
c
c
      return
      end

