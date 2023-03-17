!----------------------------------------------------------------------
!
!   routine name       - sortm
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine sorts all visible triangles in order
!                        back-to-front according to z-coordinate
!                        of their mid-points (array Rgtrz) and stores
!                        sorted numbers in the array Igtr
!
!----------------------------------------------------------------------
subroutine sortm(Nrvistr,Igtr,Rgtrz)
!
   implicit none
!
   integer, intent(in)   :: Nrvistr
!
   integer, dimension(*) :: Igtr
   real(8), dimension(*) :: Rgtrz
!
   integer :: i,j,l,index,ir
   real(8) :: q
!
!..initialize index array
   do i=1,Nrvistr
      Igtr(i) = i
   enddo
!
   if (Nrvistr.eq.1) return
!
!..start sorting (heapsort algorithm)
   l=Nrvistr/2+1
   ir=Nrvistr
!
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
         goto 700
      endif
   endif
   i=l
   j=l*2
!
30 continue
   if (j.gt.ir) goto 40
!
   if(j.lt.ir) then
      if(Rgtrz(Igtr(j)).lt.Rgtrz(Igtr(j+1))) j=j+1
   endif
!
   if (q.lt.Rgtrz(Igtr(j))) then
      Igtr(i) = Igtr(j)
      i=j
      j=j+j
   else
      j=ir+1
   endif
   goto 30
!
40 continue
   Igtr(i) = index
   goto 20
!
  700 continue
!
end subroutine sortm
