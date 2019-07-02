c----------------------------------------------------------------------
c
c   routine name       - graphg
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec,01
c
c   purpose            - routine controls interactively calls to
c                        various graphics routines for the GMP
c                        package
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine graphg
c
      use GMP
      use graphmod
c
#include "syscom.blk"
c
c  ...initial lengths of the arrays
      MXIGTR   = 100000000
      MXRGTRZ  = 9*MXIGTR
c
c  ...select window
      write(*,*) 'PLEASE SELECT KIND OF WINDOW FOR GRAPHICS'
      write(*,*) '1-B&W, 2-Gray16, 3-Gray256, 4-Color16, 5-Color256'
      read(*,*) iwin
      call initwin(iwin)
c
c  ...allocate ALL arrays used in graphics package
      allocate(IGTRCU(MXIGTR))
      allocate(IGTR(MXIGTR))
      allocate(RGTRZ(MXIGTR))
      allocate(RGTR(MXRGTRZ))
c
c
   10 continue
      write(*,*)
      write(*,*)'graphg: SELECT:'
      write(*,*)'        RETURN.............................0'
      write(*,*)'        GMP GRAPHICS.......................1'
      write(*,*)'        GMP WIREFRAME......................2'
      write(*,*)'        print_GMP..........................3'

      read(*,*) idec
c
      select case(idec)
      case(0)
c
        deallocate(IGTRCU,IGTR,RGTRZ,RGTR)
        return
c
      case(1)
        if (NDIM.eq.2) call object2(iwin)
        if (NDIM.eq.3) call object3(iwin)
      case(2)
        call GMPwireframe(iwin)
      case(3)
        call print_GMP
      end select
c
      go to 10
c
      end

