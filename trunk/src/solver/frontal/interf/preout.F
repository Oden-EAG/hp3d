c----------------------------------------------------------------------
c
c   routine name       - preout
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - an interface routine with frontal solver,
c                        for an element "nel" routine stores the
c                        element destination vectors
c
c   arguments :
c     in:
c             Iel      - number of an element
c             N        - number of the element nicknames
c             Ia       - the element nicknames
c             Ib       - the element destination vectors
c
c----------------------------------------------------------------------
c
      subroutine preout(Iel,N,Ia,Ib)
c
      use frsolmod , ONLY: IDESVE, NDESVE, NPDESV
c
      implicit none
c
      integer, intent(in) :: Iel,N
      integer, intent(in) :: Ib(*),Ia(*)
c
      integer, save :: iel_last
c
      integer :: i,np
c
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
c
      if (iprint.eq.1) then
        write(*,*) 'preout: Iel, N, iel_last',
     .                      Iel, N, iel_last
      endif
#endif
c
c  ...check if no element has been skipped (may occur with weird
c     geometry and/or BC's)
      if ((Iel.gt.1).and.(Iel.ne.iel_last+1)) then
        write(*,7002) Iel
 7002   format('preout: AN ELEMENT BEFORE Iel = ',i8,
     .         ' HAS BEEN SKIPPED')
        stop 1
      endif
c
      NDESVE(1,Iel) = N
      NDESVE(2,Iel) = NPDESV
c
      do i=1,N
        np = NPDESV
        IDESVE(NPDESV)=Ib(i)
        NPDESV = NPDESV + 1
      enddo
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Iel, (Ib(i),i=1,N)
 7001   format('preout: DEST VECTORS FOR iel = ',i8,
     .         10(/,10i10))
        call pause
      endif
#endif
c
      iel_last = Iel
c
c
      end subroutine preout
