c----------------------------------------------------------------------
c
c   routine name       - sort
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine sorts all visible triangles in order
c                        back-to-front according to z-coordinate
c                        of their mid-points (array Rgtrz) and stores
c                        sorted numbers in the array Igtr
c
c----------------------------------------------------------------------
      subroutine sortm(Nrvistr,Igtr,Rgtrz)
c
      implicit none
c
      integer, intent(in)   :: Nrvistr
c
      integer, dimension(*) :: Igtr
      real(8), dimension(*) :: Rgtrz
c
      integer :: i,j,l,index,ir
      real(8) :: q
c
c  ...initialize index array
      do i=1,Nrvistr
        Igtr(i) = i
      enddo
c
      if (Nrvistr.eq.1) return
c
c  ...start sorting (heapsort algorithm)
      l=Nrvistr/2+1
      ir=Nrvistr
   20 continue
      if (l.gt.1) then
        l=l-1
        index=Igtr(l)
        q=Rgtrz(index)
      else
        index=Igtr(ir)
        q=Rgtrz(index)
        Igtr(ir)=Igtr(1)
        ir = ir - 1
        if (ir.eq.1) then
          Igtr(1) = index

ccccc     return
          go to 700

        endif
      endif
      i=l
      j=l*2
   30 continue
      if (j.gt.ir) go to 40


cwr05.16.00
ccccc if ((j.lt.ir).and.(Rgtrz(Igtr(j)).lt.Rgtrz(Igtr(j+1)))) j=j+1
      if(j.lt.ir) then
        if(Rgtrz(Igtr(j)).lt.Rgtrz(Igtr(j+1))) j=j+1
      endif


      if (q.lt.Rgtrz(Igtr(j))) then
        Igtr(i) = Igtr(j)
        i=j
        j=j+j
      else
        j=ir+1
      endif
      go to 30
   40 continue
      Igtr(i) = index
      go to 20
c

  700 continue
c      write(*,*)'SORT: Igtr=',(Igtr(i),i=1,Nrvistr)
c      pause

      end subroutine sortm
