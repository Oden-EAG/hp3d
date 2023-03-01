!-----------------------------------------------------------------------
!
!   routine name       - sort
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine sorts elements in List in decreasing
!                        order based on values (sort key) in Val
!
!   remark             - do not use this routine for very large lists,
!                        it is a simple O(N^2) sorting algorithm
!
!   arguments :
!     in:
!             N        - number of elements to be sorted
!     inout:
!             List     - list of elements to be sorted
!             Val      - list of values (sort key)
!
!-----------------------------------------------------------------------
subroutine sort(List,Val,N)
!
   implicit none
!
   integer, intent(in   ) :: N
   integer, intent(inout) :: List(N)
   real(8), intent(inout) :: Val(N)
!
   real(8) :: valmax,valaux
   integer :: i,imax,j,laux
!
#if DEBUG_MODE
   integer :: iprint = 0
!
   if (iprint.eq.1) then
      write(*,*) 'ORIGINAL LISTS, N=', N
      do i=1,N
         write(*,7001)i,List(i),Val(i)
      enddo
   endif
#endif
!
!..loop through elements
   do i=1,N-1
      valmax= Val(i)
      imax  = i
!
!.....loop through remaining elements
      do j=i,N
         if (Val(j).gt.valmax) then
            valmax=Val(j) ; imax=j
         endif
      enddo
!
!....switch places between the i-th and the maximum value elements
     valaux   =Val(i)    ; laux      =List(i)
     Val(i)   =Val(imax) ; List(i)   =List(imax)
     Val(imax)=valaux    ; List(imax)=laux
!
!..end loop through elements
   enddo
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,*) 'SORTED OUT, N=', N
      do i=1,N
         write(*,7001) i,List(i),Val(i)
 7001    format(1x,'i,List,Val = ',i4,i6,e12.5)
      enddo
   endif
#endif
!
end subroutine sort
