c--------------------------------------------------------------------
c
c     routine name      - shape_tetrH
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of tetra H1 bubble
c                         shape functions and their derivatives
c
c     arguments:
c
c     in:
c             Xi        - master element coordinates
c             Nord      - polynomial order
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             DshapH    - derivatives of the shape functions wrt the
c                         master element coordinates
c
c-----------------------------------------------------------------------
c
      subroutine shape_tetrH(Xi,Nord, NrdofH,ShapH,DshapH)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(3),ShapH((MAXP-3)*(MAXP-2)*(MAXP-1)/6),
     .                DshapH(3,(MAXP-3)*(MAXP-2)*(MAXP-1)/6)
c
c  ...work space for routine ILegendre
      dimension u1(2:MAXP-1),du1dx(2:MAXP-1),du1dt(2:MAXP-1)
c
c  ...work space for routine IJacobi
      dimension u2(1:MAXP+1,1:MAXP),
     .          du2dx(1:MAXP+1,1:MAXP),du2dt(1:MAXP+1,1:MAXP),
     .          u3(1:MAXP+1,1:MAXP),du3dx(1:MAXP+1,1:MAXP)
      double precision z
c
c-----------------------------------------------------------------------
c      
      iprint=0
c
      x = Xi(1); y = Xi(2); z = Xi(3)
c
c  ...if no middle node bubbles, return      
      if (Nord.lt.4) then
        NrdofH=0 ; return
      endif
c
      call ILegendre(x,1.d0-y-z,Nord-1,1,
     .             u1(2:MAXP-1),du1dx(2:MAXP-1),du1dt(2:MAXP-1))
      call IJacobi(y,1.d0-z,Nord-1,Nord-2,1, u2,du2dx,du2dt)
      call IJacobi(z,1.d0,Nord-1,Nord-3,0, u3,du3dx,void)
c
      k=0
c
c  ...watch for the structure of the loop to enforce hierarchy
c     in polynomial order, np = total order
c
ccc   PG, Sep 14: do np=3,Nord
      do np=4,Nord
c
c  .....total order in x and y
        do nxy=3,np-1
c
c  .......order in x
ccc       PG, Sep 14: do i=2,np-1
ccccc     FF, Oct 14: do i=2,np-2 incorrect
          do i=2,nxy-1
c
c  .........order in y
            j = nxy-i
c
c
c  .........order in z
            l = np-nxy
c            
c  .........increment shape function counter            
            k=k+1
c            
            ShapH(k)    =    u1(i)  * u2(i,j)  * u3(nxy,l)
            DshapH(1,k) =   du1dx(i)* u2(i,j)  * u3(nxy,l)
            DshapH(2,k) = - du1dt(i)* u2(i,j)  * u3(nxy,l)
     .                    +  u1(i)  *du2dx(i,j)* u3(nxy,l)
            DshapH(3,k) = - du1dt(i)* u2(i,j)  * u3(nxy,l)
     .                    -  u1(i)  *du2dt(i,j)* u3(nxy,l)
     .                    +  u1(i)  * u2(i,j)  *du3dx(nxy,l)
          enddo
        enddo
      enddo
      NrdofH=k
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:2),Nord
 7001   format('shap_tetrH: Xi = ',2f8.3,' Nord = ',i2)
        do k=1,NrdofH
          write(*,7002) k,ShapH(k),DshapH(1:3,k)
 7002     format('k = ',i3,' ShapH, DshapH = ',e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape_tetrH
