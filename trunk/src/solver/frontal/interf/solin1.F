c------------------------------------------------------------------
c
c   routine name       - solin1
c
c------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - a frontal solver interface routine
c
c   arguments :
c     in:
c              Iel     - an element number
c     out:
c              Numdes  - number of destination vectors
c              Zaldest - destination vectors
c
c----------------------------------------------------------------------
#include "typedefs.h"
      subroutine solin1(Iel, Numdes,Zaldest)
c
      use frsolmod , ONLY: IDESVE, NDESVE
      implicit none
c
      integer, intent(in)  :: Iel
      integer, intent(out) :: Numdes
c
      VTYPE  , intent(out) :: Zaldest(*)
c
      integer :: i,iaux,naux,np
c
#if DEBUG_MODE
      integer :: iprint
      iprint=0
#endif
c
      Numdes = NDESVE(1,Iel)
      np = NDESVE(2,Iel)
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,*) 'solin1: Numdes,np = ',Numdes,np
        write(*,*) 'solin1: NDESVE(1,2403),NDESVE(2,2403) = ',
     .                      NDESVE(1,2403),NDESVE(2,2403)
      endif
#endif
c
      do i=1,Numdes
        iaux = IDESVE(np+i-1)
#if C_MODE
        Zaldest(i) = dcmplx(dfloat(IDESVE(np+i-1)),0.d0)
#else
        Zaldest(i) = dfloat(IDESVE(np+i-1))
#endif
c
c  .....double check
        naux = INT(Zaldest(i))
        if (naux.ne.IDESVE(np+i-1)) then
          write(*,*) 'solin1: IDESVE(np+i-1),Zaldest(i) = ',
     .                        IDESVE(np+i-1),Zaldest(i)
          stop 1
        endif
      enddo
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7001) Iel
 7001   format('solin1: DESTINATION VECTORS FOR Iel = ',i8)
        write(*,7003) (IDESVE(np+i-1),i=1,Numdes)
 7003   format(5(6x,i10))
        write(*,7002) (Zaldest(i),i=1,Numdes)
#if C_MODE
 7002   format(5('(',f10.1,f4.1,')',2x))
#else
 7002   format(5(f10.1,2x))
#endif
        call pause
      endif
#endif
c
      end subroutine solin1
