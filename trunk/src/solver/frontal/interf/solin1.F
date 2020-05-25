c------------------------------------------------------------------
c
c   routine name       - solin1
c
c------------------------------------------------------------------
c
c   latest revision    - Jan 08
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
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine solin1(Iel, Numdes,Zaldest)
c
      use frsolmod , ONLY: IDESVE, NDESVE
#include "syscom.blk"
c
      dimension Zaldest(*)
#if C_MODE
      complex(8) dcmplx
#endif
c
      if (Iel.ge.1) then
        iprint=0
      else
        iprint=0
      endif
c
      Numdes = NDESVE(1,Iel)
      np = NDESVE(2,Iel)
      if (iprint.eq.1) then
        write(*,*) 'solin1: Numdes,np = ',Numdes,np
        write(*,*) 'solin1: NDESVE(1,2403),NDESVE(2,2403) = ',
     .                      NDESVE(1,2403),NDESVE(2,2403)
      endif
c
      do i=1,Numdes
        iaux = IDESVE(np+i-1)
#if C_MODE
        zaux  = dcmplx(dble(float( IDESVE(np+i-1))),0.d0)
c
cldem, jan08, this costed me three days of work before I caught it
ccc        ZAldest(i)=dcmplx(dble(float(IDESVE(np+i-1))),0.d0)
        ZAldest(i)=dcmplx(dfloat(IDESVE(np+i-1)),0.d0)
ccc        write(*,*) 'IDESVE(np+i-1) = ',IDESVE(np+i-1)
ccc        write(*,*) 'dfloat(IDESVE(np+i-1)) = ',dfloat(IDESVE(np+i-1))
ccc        write(*,*) 'dcmplx(dfloat(IDESVE(np+i-1)),0.d0) = ',
ccc     .              dcmplx(dfloat(IDESVE(np+i-1)),0.d0)
#else
        aux  = float( IDESVE(np+i-1) )
ccc        Zaldest(i)= float( IDESVE(np+i-1) )
        Zaldest(i)= dfloat(IDESVE(np+i-1))
#endif
c
c  .....double check
        naux = Zaldest(i)
        if (naux.ne.IDESVE(np+i-1)) then
          write(*,*) 'solin1: IDESVE(np+i-1),Zaldest(i) = ',
     .                        IDESVE(np+i-1),Zaldest(i)
          stop 1
        endif
      enddo
c
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
c
c
      end
