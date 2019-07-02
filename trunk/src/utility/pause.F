c----------------------------------------------------------------------
c
c   routine name       - pause
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 02
c
c   purpose            - routine emulates fortran 'pause'
c                        avoind stupid realizations in some
c                        compilers
c
c   arguments          - none
c
c   required  routines - 
c
c----------------------------------------------------------------------
c
      subroutine pause
c
      use parallelism
c
      write(*,*) 'PAUSE - type anything'
      if(MYRANK.eq.0)then
      read(*,*)
      endif
c
      return
      end




       

