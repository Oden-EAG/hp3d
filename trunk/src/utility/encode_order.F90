!-----------------------------------------------------------------------
!
!     routine name      - encode_orderb
!
!-----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine encodes order of approximation for
!                         a brick middle node
!
!     arguments:
!       in:
!             Nord1,2,3 - orders in directions 1,2,3
!       out:
!             Norderb   - the encoded order
!
!-----------------------------------------------------------------------
subroutine encode_orderb(Nord1,Nord2,Nord3, Norderb)
!
   use parameters , only : MAXP
!
   implicit none
!
   integer, intent(in)  :: Nord1,Nord2,Nord3
   integer, intent(out) :: Norderb
!
   Norderb = (Nord1*(MAXP+1)+Nord2)*(MAXP+1)+Nord3
!
#if HP3D_DEBUG
   if (Nord1<1.or.Nord1>MAXP) then
      write(*,*)'encode_orderb: OUT OF RANGE Nord1 = ',Nord1
      stop
   endif
   if (Nord2<1.or.Nord2>MAXP) then
      write(*,*)'encode_orderb: OUT OF RANGE Nord2 = ',Nord2
      stop
   endif
   if (Nord3<1.or.Nord3>MAXP) then
      write(*,*)'encode_orderb: OUT OF RANGE Nord3 = ',Nord3
      stop
   endif
#endif
!
end subroutine encode_orderb
!
!-----------------------------------------------------------------------
!
!     routine name      - encode_orderq
!
!-----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine encodes order of approximation for
!                         a quad middle node
!
!     arguments:
!       in:
!             Nord1,2   - orders in directions 1,2
!       out:
!             Norderq   - the encoded order
!
!-----------------------------------------------------------------------
subroutine encode_orderq(Nord1,Nord2, Norderq)
!
   use parameters , only : MAXP
!
   implicit none
!
   integer, intent(in)  :: Nord1,Nord2
   integer, intent(out) :: Norderq
!
   Norderq = Nord1*(MAXP+1)+Nord2
!
#if HP3D_DEBUG
   if (Nord1<1.or.Nord1>MAXP) then
      write(*,*)'encode_orderq: OUT OF RANGE Nord1 = ',Nord1
      stop
   endif
   if (Nord2<1.or.Nord2>MAXP) then
      write(*,*)'encode_orderq: OUT OF RANGE Nord2 = ',Nord2
      stop
   endif
#endif
!
end subroutine encode_orderq
