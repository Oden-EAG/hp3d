c----------------------------------------------------------------------
c
c   routine name       - shape1H_0
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate 1D H1 shape functions
c
c   arguments :
c     in:
c        Xi            - master element coordinates
c        Nord          - order of approximation
c     out:
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shape1H_0(Xi,Nord, NrdofH,ShapH,DshapH)
c
      use parameters
      use element_data
      implicit none
c
      integer :: Nord,NrdofH
      real(8) :: Xi,ShapH(MAXP+1),DshapH(MAXP+1)
c
c  ...kernel shape functions
      real(8) :: shapHk(MAXP-1),dshapHk(MAXP-1)
c
      real(8) :: t,dtdxi
      integer :: j,k,norient,nrdofHk
c
c----------------------------------------------------------------------
c
      integer :: iprint=0
c
      NrdofH = Nord+1
c
c  ...calculate the linear shape functions
      ShapH(1) = 1.d0 - Xi
      DshapH(1) = -1.d0
c
      ShapH(2) = Xi
      DshapH(2) = 1.d0
c
      t = Xi; dtdxi = 1.d0
c
c  ...initiate shape function counter
      k=2
c
c  ...middle node shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  ...evaluate 1D kernel shape functions
      call shapHke(t,Nord,0, nrdofHk,shapHk,dshapHk)
c
c  ...evaluate middle node shape functions
      do j=1,nrdofHk
        k=k+1
        ShapH(k)  =  ShapH(1)* ShapH(2)* shapHk(j)
        DshapH(k) = DshapH(1)* ShapH(2)* shapHk(j)
     .            +  ShapH(1)*DshapH(2)* shapHk(j)
     .            +  ShapH(1)* ShapH(2)*dshapHk(j)*dtdxi
      enddo
c
      if (iprint.eq.1) then
        write(*,7001) Xi,Nord,Norient
 7001   format('shape1H: Xi,Nord,Norient = ',e12.5,3x,2i3)
        do k=1,Nord+1
          write(*,7002) k,ShapH(k),DshapH(k)
 7002     format('k = ',i3,2x,2e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape1H_0
