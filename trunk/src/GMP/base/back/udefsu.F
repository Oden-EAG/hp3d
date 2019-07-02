c----------------------------------------------------------------------
c
c   routine name       - udefsu
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - June,1999
c
c   purpose            - provides an oppurtunity for a user
c                        to define his surface
c
c   usage              - udefsu(Label,X, Fval,Dfdx)
c
c   arguments :
c     in:              - Label,X
c     out:             - Fval,Dfdx
c
c   required  routines - 
c
c----------------------------------------------------------------------
      subroutine udefsu(Label,X, Fval,Dfdx)
c
#include "syscom.blk"
      dimension X(3), Dfdx(3)
      save irout1, irout2
      dimension rs(3),ys(3)

      rs(1) = .7d0
      rs(2) = .83d0
      rs(3) = 1.3d0
      r = .06d0
      ys(1) = .19d0
      ys(2) = .34d0
      ys(3) = .42d0
c-----------------------------
c------------------------------
      rsph2 = rs(2)**2 
      reli2 = rs(1)**2
c  ...half-thickness of wings:
      a = .025
c  ...half-width of wings:
      d = (ys(2)-ys(1))/2
c  ...half-span of big wings:
      hb = rs(2)
c  ...half-span of small wings:
      hs = rs(1)
c  ...shift of wing:
      shift = ys(1) + d
c  ...elipsoid with radious = .5 and the center at (0,0,0)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if( Label==101 )then
         if(irout1.ne.100)then
           irout1 = 100
           write(*,*)'UDEFSU: ellipsoid, radius=.5'
           call pause
         endif
         Fval = 4.d0*X(1)*X(1) + X(2)*X(2) + 4.d0*X(3)*X(3) - reli2
c
         Dfdx(1) = 2*X(1) *4.d0
         Dfdx(2) = 2*X(2)
         Dfdx(3) = 2*X(3) *4.d0
c
c  ...sphere with radious = .83d0 and the center at (0,0,0)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      elseif( Label==102 )then
         if(irout2.ne.100)then
           irout2 = 100
           write(*,*)'UDEFSU: sphere radius=.83'
           call pause
         endif
         Fval = X(1)*X(1) + X(2)*X(2) + X(3)*X(3) - rsph2
c
         Dfdx(1) = 2*X(1) 
         Dfdx(2) = 2*X(2)
         Dfdx(3) = 2*X(3) 
c  ...upper surface of wing:
c     ~~~~~~~~~~~~~~~~~~~~~
      elseif( Label==103 )then
         sq2 = sqrt(2.d0)/2.d0
         h = hb
         x(2) = x(2)+shift
         Fval = sq2*(x(1)+x(3)) 
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * (1-.5*(x(1)-x(3))**2 / (h**2-x(2)**2))
c

         Dfdx(1) = sq2 + a * (1-x(2)**2/d**2) * (x(1)-x(3))
     .                                              / (h**2-x(2)**2)

         Dfdx(2) = -a *(
     .
     .             -2*x(2)/d**2 * (1-.5*(x(1)-x(3))**2 / (h**2-x(2)**2))
     .
     .            - (1-x(2)**2/d**2) * (x(1)-x(3))**2
     .                                / (h**2-x(2)**2)**2 *x(2)
     .                 )

         Dfdx(3) = 
     .          sq2
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * ((x(1)-x(3)) / (h**2-x(2)**2))

c     ...shift y:
         x(2) = x(2)-shift
c  ...lower surface of wing:
c~~~~~~~~~~~~~~~~~~~~~~~~~~
      elseif( Label==104 )then
         sq2 = sqrt(2.d0)/2.d0
          h = hb
c     ...shift y:
         x(2) = x(2)+shift
         Fval = -sq2*(x(1)+x(3)) 
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * (1-.5*(x(1)-x(3))**2 / (h**2-x(2)**2))
c

         Dfdx(1) = -sq2 + a * (1-x(2)**2/d**2) * (x(1)-x(3))
     .                                              / (h**2-x(2)**2)

         Dfdx(2) = -a *(
     .
     .             -2*x(2)/d**2 * (1-.5*(x(1)-x(3))**2 / (h**2-x(2)**2))
     .
     .            - (1-x(2)**2/d**2) * (x(1)-x(3))**2
     .                                / (h**2-x(2)**2)**2 *x(2)
     .                 )

         Dfdx(3) = 
     .          -sq2
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * ((x(1)-x(3)) / (h**2-x(2)**2))

c     ...shift y:
         x(2) = x(2)-shift
c  ...upper surface of small wing:
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      elseif( Label==105 )then
         sq2 = sqrt(2.d0)/2.d0
          h = hs
c     ...shift y:
         x(2) = x(2)-shift



         Fval = sq2*(x(1)+x(3)) 
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * (1-2*(x(1)-x(3))**2 / (h**2-x(2)**2))
c

         Dfdx(1) = sq2 + a * (1-x(2)**2/d**2) * 4*(x(1)-x(3))
     .                                              / (h**2-x(2)**2)

         Dfdx(2) = -a *(
     .
     .             -2*x(2)/d**2 * (1-2*(x(1)-x(3))**2 / (h**2-x(2)**2))
     .
     .            - (1-x(2)**2/d**2) * 4*(x(1)-x(3))**2
     .                                / (h**2-x(2)**2)**2 *x(2)
     .                 )

         Dfdx(3) = 
     .          sq2
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * 4*((x(1)-x(3)) / (h**2-x(2)**2))

c     ...shift y:
ccccc    x(2) = x(2)+.145    
         x(2) = x(2)+shift



c  ...lower surface of small wing:
      elseif( Label==106 )then


         sq2 = sqrt(2.d0)/2.d0

ccccc    a = .025
ccccc    d = .105
ccccc    h = .5
         h = hs

c     ...shift y:
ccccc    x(2) = x(2)-.145
         x(2) = x(2)-shift



         Fval = -sq2*(x(1)+x(3)) 
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * (1-2*(x(1)-x(3))**2 / (h**2-x(2)**2))
c

         Dfdx(1) = -sq2 + a * (1-x(2)**2/d**2) * 4*(x(1)-x(3))
     .                                              / (h**2-x(2)**2)

         Dfdx(2) = -a *(
     .
ccccc    x(2) = x(2)-.145
     .             -2*x(2)/d**2 * (1-2*(x(1)-x(3))**2 / (h**2-x(2)**2))
     .
     .            - (1-x(2)**2/d**2) * 4*(x(1)-x(3))**2
     .                                / (h**2-x(2)**2)**2 *x(2)
     .                 )

         Dfdx(3) = 
     .          -sq2
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * 4*((x(1)-x(3)) / (h**2-x(2)**2))

c     ...shift y:
ccccc    x(2) = x(2)+.145    
         x(2) = x(2)+shift


c  ...right surface of vertical small wing:
      elseif( Label==107 )then

         sq2 = sqrt(2.d0)/2.d0

ccccc    a = .025
ccccc    d = .105
ccccc    h = .5
         h = hs

c     ...shift y:
ccccc    x(2) = x(2)-.145
         x(2) = x(2)-shift



         Fval = sq2*(-x(1)+x(3)) 
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * (1-2*(x(1)+x(3))**2 / (h**2-x(2)**2))
c

         Dfdx(1) = -sq2 + a * (1-x(2)**2/d**2) * 4*(x(1)-x(3))
     .                                              / (h**2-x(2)**2)

         Dfdx(2) = -a *(
     .
     .             -2*x(2)/d**2 * (1-2*(x(1)+x(3))**2 / (h**2-x(2)**2))
     .
     .            - (1-x(2)**2/d**2) * 4*(x(1)+x(3))**2
     .                                / (h**2-x(2)**2)**2 *x(2)
     .                 )

         Dfdx(3) = 
     .          sq2
     .
     .        + a * (1-x(2)**2/d**2)
     .
     .          * 4*((x(1)-x(3)) / (h**2-x(2)**2))

c     ...shift y:
ccccc    x(2) = x(2)+.145    
         x(2) = x(2)+shift



c  ...left surface of vertical small wing:
      elseif( Label==108 )then

         sq2 = sqrt(2.d0)/2.d0

ccccc    a = .025
ccccc    d = .105
ccccc    h = .5
         h = hs

c     ...shift y:
ccccc    x(2) = x(2)-.145
         x(2) = x(2)-shift



         Fval = -sq2*(-x(1)+x(3)) 
     .
     .        - a * (1-x(2)**2/d**2)
     .
     .          * (1-2*(x(1)+x(3))**2 / (h**2-x(2)**2))
c

         Dfdx(1) = sq2 + a * (1-x(2)**2/d**2) * 4*(x(1)-x(3))
     .                                              / (h**2-x(2)**2)

         Dfdx(2) = -a *(
     .
     .             -2*x(2)/d**2 * (1-2*(x(1)+x(3))**2 / (h**2-x(2)**2))
     .
     .            - (1-x(2)**2/d**2) * 4*(x(1)+x(3))**2
     .                                / (h**2-x(2)**2)**2 *x(2)
     .                 )

         Dfdx(3) = 
     .          -sq2
     .
     .        + a * (1-x(2)**2/d**2)
     .
     .          * 4*((x(1)-x(3)) / (h**2-x(2)**2))

c     ...shift y:
ccccc    x(2) = x(2)+.145    
         x(2) = x(2)+shift


c  ...nose:
      elseif( Label==109 )then

         hnose = 1.5*d

         Fval = -x(2) -hnose*(1-x(1)**2/r**2) * (1-x(3)**2/r**2) -ys(3)
c

         Dfdx(1) =  hnose * 2*x(1)/r**2 * (1-x(3)**2/r**2) 

         Dfdx(2) = -1

         Dfdx(3) =  hnose * (1-x(1)**2/r**2) * 2*x(3)/r**2

      else
         print*,'UDEFSU: UNDEFINED SURFACE Label = ',Label
         stop
      endif

      return
c---------------------------------------------------------------------






c  ...sphere with radious = .5 and the center at (0,0,0)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if( Label==101 )then

         if(irout1.ne.100)then
           irout1 = 100
           write(*,*)'UDEFSU: sphere radius=.5'
           call pause
         endif

         Fval = X(1)*X(1) + X(2)*X(2) + X(3)*X(3) - .25d0
c
         Dfdx(1) = 2*X(1)
         Dfdx(2) = 2*X(2)
         Dfdx(3) = 2*X(3)
c
c  ...sphere with radious = .83d0 and the center at (0,0,0)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      elseif( Label==102 )then

         if(irout2.ne.100)then
           irout2 = 100
           write(*,*)'UDEFSU: sphere radius=.83'
           call pause
         endif

         Fval = X(1)*X(1) + X(2)*X(2) + X(3)*X(3) - .6889d0
c
         Dfdx(1) = 2*X(1)
         Dfdx(2) = 2*X(2)
         Dfdx(3) = 2*X(3)

      else
         print*,'UDEFSU: UNDEFINED SURFACE Label = ',Label
         stop
      endif

      return
c---------------------------------------------------------------------


c  ...sphere with radious = sqrt(3) and the center at (0,0,0)
      if( Label==101 )then
ccccc    Fval = X(1)*X(1) + X(2)*X(2) + X(3)*X(3) - 3.d0
         Fval = X(1)*X(1) + X(2)*X(2) + X(3)*X(3) - 1.d0
c
         Dfdx(1) = 2*X(1)
         Dfdx(2) = 2*X(2)
         Dfdx(3) = 2*X(3)
c
c  ...sphere with radious = sqrt(12) and the center at (0,0,0)
      elseif( Label==102 )then
ccccc    Fval = X(1)*X(1) + X(2)*X(2) + X(3)*X(3) - 12.d0
         Fval = X(1)*X(1) + X(2)*X(2) + X(3)*X(3) - .5625
c
         Dfdx(1) = 2*X(1)
         Dfdx(2) = 2*X(2)
         Dfdx(3) = 2*X(3)

      elseif(Label.eq.103)then
         Fval = X(1)*X(1) + X(2)*X(2) + X(3) - 2.d0

         Dfdx(1) = 2*X(1)
         Dfdx(2) = 2*X(2)
         Dfdx(3) = 1.d0
 

      elseif(Label.eq.104)then
         Fval = X(1)*X(1) + X(2)*X(2) + X(3) - 1.d0

         Dfdx(1) = 2*X(1)
         Dfdx(2) = 2*X(2)
         Dfdx(3) = 1.d0


      else
         print*,'UDEFSU: UNDEFINED SURFACE Label = ',Label
         stop
      endif
        
      return
      end






