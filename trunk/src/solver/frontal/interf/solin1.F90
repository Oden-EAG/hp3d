!------------------------------------------------------------------
!
!   routine name       - solin1
!
!------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - a frontal solver interface routine
!
!   arguments :
!     in:
!              Iel     - an element number
!     out:
!              Numdes  - number of destination vectors
!              Zaldest - destination vectors
!
!----------------------------------------------------------------------
#include "typedefs.h"
   subroutine solin1(Iel, Numdes,Zaldest)
!
      use frsolmod , ONLY: IDESVE, NDESVE
      implicit none
!
      integer, intent(in)  :: Iel
      integer, intent(out) :: Numdes
!
      VTYPE  , intent(out) :: Zaldest(*)
!
      integer :: i,iaux,naux,np
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      Numdes = NDESVE(1,Iel)
      np = NDESVE(2,Iel)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'solin1: Numdes,np = ',Numdes,np
        write(*,*) 'solin1: NDESVE(1,2403),NDESVE(2,2403) = ', &
                            NDESVE(1,2403),NDESVE(2,2403)
      endif
#endif
!
      do i=1,Numdes
        iaux = IDESVE(np+i-1)
#if HP3D_COMPLEX
        Zaldest(i) = dcmplx(dfloat(IDESVE(np+i-1)),0.d0)
#else
        Zaldest(i) = dfloat(IDESVE(np+i-1))
#endif
!
!  .....double check
        naux = INT(Zaldest(i))
        if (naux.ne.IDESVE(np+i-1)) then
          write(*,*) 'solin1: IDESVE(np+i-1),Zaldest(i) = ', &
                              IDESVE(np+i-1),Zaldest(i)
          stop 1
        endif
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Iel
 7001   format('solin1: DESTINATION VECTORS FOR Iel = ',i8)
        write(*,7003) (IDESVE(np+i-1),i=1,Numdes)
 7003   format(5(6x,i10))
        write(*,7002) (Zaldest(i),i=1,Numdes)
#if HP3D_COMPLEX
 7002   format(5('(',f10.1,f4.1,')',2x))
#else
 7002   format(5(f10.1,2x))
#endif
        call pause
      endif
#endif
!
   end subroutine solin1
