!----------------------------------------------------------------------
!
!   routine name       - graphg
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine controls interactively calls to
!                        various graphics routines for the GMP
!                        package
!
!----------------------------------------------------------------------

#if HP3D_USE_X11

!
   subroutine graphg
!
      use GMP
      use graphmod
!
      implicit none
!
      integer :: idec,iwin
!
!  ...initial lengths of the arrays
      MXIGTR   = 100000000
      MXRGTRZ  = 9*MXIGTR
!
!  ...select window
      write(*,*) 'PLEASE SELECT KIND OF WINDOW FOR GRAPHICS'
      write(*,*) '1-B&W, 2-Gray16, 3-Gray256, 4-Color16, 5-Color256'
      read(*,*) iwin
      call initwin(iwin)
!
!  ...allocate ALL arrays used in graphics package
      allocate(IGTRCU(MXIGTR))
      allocate(IGTR(MXIGTR))
      allocate(RGTRZ(MXIGTR))
      allocate(RGTR(MXRGTRZ))
!
!
   10 continue
      write(*,*)
      write(*,*)'graphg: SELECT:'
      write(*,*)'        RETURN.............................0'
      write(*,*)'        GMP GRAPHICS.......................1'
      write(*,*)'        GMP WIREFRAME......................2'
      write(*,*)'        print_GMP..........................3'

      read(*,*) idec
!
      select case(idec)
      case(0)
!
        deallocate(IGTRCU,IGTR,RGTRZ,RGTR)
        return
!
      case(1)
        if (NDIM.eq.2) call object2(iwin)
        if (NDIM.eq.3) call object3(iwin)
      case(2)
        call GMPwireframe(iwin)
      case(3)
        call print_GMP
      end select
!
      go to 10
!
   end subroutine graphg

#else

   subroutine graphg
      
      implicit none
      
      write(*,*) 'graphg is only supported with HP3D_USE_X11 enabled.'
      
   end subroutine graphg

#endif
