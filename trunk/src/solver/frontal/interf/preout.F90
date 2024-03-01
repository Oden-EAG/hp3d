!----------------------------------------------------------------------
!
!   routine name       - preout
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - an interface routine with frontal solver,
!                        for an element "nel" routine stores the
!                        element destination vectors
!
!   arguments :
!     in:
!             Iel      - number of an element
!             N        - number of the element nicknames
!             Ia       - the element nicknames
!             Ib       - the element destination vectors
!
!----------------------------------------------------------------------
!
   subroutine preout(Iel,N,Ia,Ib)
!
      use frsolmod , ONLY: IDESVE, NDESVE, NPDESV
!
      implicit none
!
      integer, intent(in) :: Iel,N
      integer, intent(in) :: Ib(*),Ia(*)
!
      integer, save :: iel_last
!
      integer :: i,np
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
!
      if (iprint.eq.1) then
        write(*,*) 'preout: Iel, N, iel_last', &
                            Iel, N, iel_last
      endif
#endif
!
!  ...check if no element has been skipped (may occur with weird
!     geometry and/or BC's)
      if ((Iel.gt.1).and.(Iel.ne.iel_last+1)) then
        write(*,7002) Iel
 7002   format('preout: AN ELEMENT BEFORE Iel = ',i8, &
               ' HAS BEEN SKIPPED')
        stop 1
      endif
!
      NDESVE(1,Iel) = N
      NDESVE(2,Iel) = NPDESV
!
      do i=1,N
        np = NPDESV
        IDESVE(NPDESV)=Ib(i)
        NPDESV = NPDESV + 1
      enddo
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Iel, (Ib(i),i=1,N)
 7001   format('preout: DEST VECTORS FOR iel = ',i8, &
               10(/,10i10))
        call pause
      endif
#endif
!
      iel_last = Iel
!
!
   end subroutine preout
