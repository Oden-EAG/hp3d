!-----------------------------------------------------------------------
!
!     routine name      - decode_orderb
!
!-----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine decodes order of approximation for
!                         a brick middle node
!
!     arguments:
!       in:
!             Norderb   - the encoded order
!       out:
!             Nord1,2,3 - orders in directions 1,2,3
!
!-----------------------------------------------------------------------
subroutine decode_orderb(Norderb, Nord1,Nord2,Nord3)
!
   use parameters , only : MAXP
   implicit none
!
   integer, intent(in)  :: Norderb
   integer, intent(out) :: Nord1,Nord2,Nord3
   integer :: naux
!
!  RECALL: Norderb = (Nord1*(MAXP+1)+Nord2)*(MAXP+1)+Nord3
!
   naux  = Norderb/(MAXP+1)
   Nord3 = Norderb-naux*(MAXP+1)
   Nord1 = naux/(MAXP+1)
   Nord2 = naux-Nord1*(MAXP+1)
!
#if HP3D_DEBUG
   if (Nord1<1.or.Nord1>MAXP) then
      write(*,*)'decode_orderb: OUT OF RANGE Nord1 = ',Nord1
      stop
   endif
   if (Nord2<1.or.Nord2>MAXP) then
      write(*,*)'decode_orderb: OUT OF RANGE Nord2 = ',Nord2
      stop
   endif
   if (Nord3<1.or.Nord3>MAXP) then
      write(*,*)'decode_orderb: OUT OF RANGE Nord3 = ',Nord3
      stop
   endif
#endif
!
end subroutine decode_orderb
!
!
!-----------------------------------------------------------------------
!
!     routine name      - decode_orderq
!
!-----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine decodes order of approximation for
!                         a quad middle node
!
!     arguments:
!
!     in:
!             Norderq   - the encoded order
!     out:
!             Nord1,2   - orders in directions 1,2
!
!-----------------------------------------------------------------------
subroutine decode_orderq(Norderq, Nord1,Nord2)
!
   use parameters , only : MAXP
   implicit none
!
   integer, intent(in)  :: Norderq
   integer, intent(out) :: Nord1,Nord2
!
!  RECALL: Norderq = Nord1*(MAXP+1)+Nord2
!
   Nord1 = Norderq/(MAXP+1)
   Nord2 = Norderq-Nord1*(MAXP+1)
!
#if HP3D_DEBUG
   if (Nord1<1.or.Nord1>MAXP) then
      write(*,*)'decode_orderq: OUT OF RANGE Nord1 = ',Nord1
      stop
   endif
   if (Nord2<1.or.Nord2>MAXP) then
      write(*,*)'decode_orderq: OUT OF RANGE Nord2 = ',Nord2
      stop
   endif
#endif
!
end subroutine decode_orderq
