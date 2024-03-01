#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - sortz
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine sorts all visible triangles in order
!                        back-to-front according to z-coordinate
!                        of their mid-points (array RGTRZ) and stores
!                        sorted numbers in the array IGTR
!
!----------------------------------------------------------------------
!
   subroutine sortz
!
      use graphmod
!
      implicit none
!
      integer :: i,index,ir,j,l
      real(8) :: q
!
#if HP3D_DEBUG
      integer :: ii
      integer :: iprint
      iprint=0
#endif
!
      if (NRVISTR.eq.0) return
!
!  ...initialize index array
      do 10 i=1,NRVISTR
        IGTR(i) = i
 10   continue
!
!  ...start sorting (heapsort algorithm)
      l=NRVISTR/2+1
      ir=NRVISTR
   20 continue
      if (l.gt.1) then
        l=l-1
        index=IGTR(l)
        q = RGTRZ(index)
      else
        index=IGTR(ir)
        q = RGTRZ(index)
        IGTR(ir) = IGTR(1)
        ir = ir - 1
        if (ir.eq.1) then
          IGTR(1) = index
#if HP3D_DEBUG
          if (iprint.eq.1) then
            write(*,7001) (IGTR(ii),ii=1,NRVISTR)
 7001       format(20i6)
            call pause
          endif
#endif
          return
        endif
      endif
      i=l
      j=l*2
   30 continue
      if (j.gt.ir) go to 40
      if (j.lt.ir) then
        if (RGTRZ(IGTR(j)).lt.RGTRZ(IGTR(j+1))) j=j+1
      endif
      if (q.lt.RGTRZ(IGTR(j))) then
        IGTR(i) = IGTR(j)
        i=j
        j=j+j
      else
        j=ir+1
      endif
      go to 30
   40 continue
      IGTR(i) = index
      go to 20
!
   end subroutine sortz

#endif
