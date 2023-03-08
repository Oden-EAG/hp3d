c----------------------------------------------------------------------
c
c   latest revision     - Feb 2023
c
c   purpose             - routine selects specific display
c
c   argument in : Iwind - type of display
c
c   required routines   - openwind, closwind, clrwind
c
c----------------------------------------------------------------------
c
      subroutine selwin(Iwind)
c
      use graphmod
c
      implicit none
c
      integer, intent(in) :: Iwind
c
      integer :: iwindtype = 0
c
      IWINDNUM = 1
c
c  ...these parameters are hardwired in the interface with X-windows;
c     if you want to change them, you will have to change them
c     in the interface as well !!
ccc      iwindl = 750
ccc      iwindl = 1000
ccc      iwindh = 600
ccc      iwindl = 667
ccc      iwindh = 400
c      iwindl = 1334
c      iwindh = 800
      rwindl = iwindl
      rwindh = iwindh
      ISIZE(1) = 0
      ISIZE(2) = 0
      ISIZE(3) = iwindl
      ISIZE(4) = iwindh
      call openwind(iwindtype,ISIZE,IWINDNUM)
      call setpost(1)
c
c  ...margin
      rmargin = .05d0*rwindh
c
c  ...size of the window after subtracting margin...
      xlength = rwindl - 3.d0*rmargin
      ylength = rwindh - 2.d0*rmargin
c
      call clrwind(npcol(1))
c
c
      end subroutine selwin
