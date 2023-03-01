!----------------------------------------------------------------------
!
!   routine name       - legendre_S_I
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 10
!
!   purpose            - evaluates Scaled Integrated Legendre-type
!                        polynomials and their derivatives
!
!   arguments :
!     in:
!                X,T   - variables
!                Nord  - order of polynomial
!     out:
!                Poly  - polynomials
!                Dpoly - derivatives
!
!----------------------------------------------------------------------

subroutine legendre_S_I(X,T,Nord, Poly,Dpoly)
  !
  use parameters
#include "syscom.blk"
  !
  !  ...exceed maximum order of approximation
  integer, parameter :: ndelta = 2
  !
  !  ...input output parameters
  dimension Poly(0:MAXP+ndelta),Dpoly(2,0:MAXP+ndelta)
  !
  !----------------------------------------------------------------------
  !
  ! ....test prints:
  iprint = 0
  !
  if (Nord.gt.MAXP+ndelta) then
     write(*,*)'legendre_S_I: Nord = ',Nord
     stop
  endif
  !
  !----------------------------------------------------------------------
  !  ...from LUCY'S routine..............................................
  !
  !  ...scaled legendre polynomials only with n>=2 used
  !     L^s_1 = x
  !     L^s_2 = 0.5 (x^2-t^2)
  !     L^s_{n+1} = (2n-1)/(n+1) x L^s_{n} - (n-2)/(n+1) t^2 L^s_{n-1}
  !
  !  ...gradients wrt (x,t)
  !     dL^s_1/dx = 1
  !     dL^s_2/dx = x
  !     dL^s_{n+1}/dx = (2n-1)/(n+1) ( L^s_n + dL^s_n )
  !                    -(n-2)/(n+1) t^2 dL^s_{n-1}/dx
  !     dL^s_1/ds = 0
  !     dL^s_2/ds =-t
  !     dL^s_{n+1}/ds = (2n-1)/(n+1) x dL^s_n/dt
  !                    -(n-2)/(n+1)(2t L^s_{n-1}+t^2 dL^s_{n-1}/dt)
  !----------------------------------------------------------------------
  !
  !
  !  ...initialize
  Poly = 0.d0; Dpoly = 0.d0
  !
  !  ...L^s_1
  Poly(1) = x
  Dpoly(1,1) = 1.d0
  !
  !  ...L^s_2, L^s_3
  Poly(2) = 0.5d0*(x*x - t*t)
  Poly(3) = 0.5d0*x*(x*x - t*t)
  Dpoly(1,2)= x; Dpoly(2,2) = -t
  Dpoly(1,3)= 3.d0/2.d0*x*x - 0.5d0*t*t; Dpoly(2,3)= -x*t
  !
  !  ...L^s_j, 4 =< j =< Nord
  do n=3,Nord-1
     !  PG, Dec 10: double check c1 and c2
     c1 = (2.d0*n-1.d0)/(n+1.d0)
     c2 = (n-2.d0)/(n+1.d0)
     !cc        c1 = (2.d0*n-3.d0)/n
     !cc        c2 = (n-3.d0)/n
     Poly(   n+1) = c1*x*Poly(n)-c2*t*t*Poly(n-1)
     Dpoly(1,n+1) = c1*(Poly(n)+x*Dpoly(1,n))-c2*t*t*Dpoly(1,n-1)
     Dpoly(2,n+1) = c1*x*Dpoly(2,n)          &
          -c2*( 2.d0*t*Poly(n-1) &
          +t*t*Dpoly(2,n-1) )
  enddo
  !
  !
  !ccc  ...L^s_2, L^s_3
  !cc      Poly(1) = 0.5d0*(x*x - t*t)
  !cc      Poly(2) = 0.5d0*x*(x*x - t*t)
  !cc      Dpoly(1,1)= x; Dpoly(2,1) = -t
  !cc      Dpoly(1,2)= 3.d0/2.d0*x*x - 0.5d0*t*t; Dpoly(2,2)= -x*t
  !ccc
  !ccc  ...L^s_j, 4 =< j =< Nord
  !cc      do n=2,Nord-1
  !cc        c1 = (2.d0*n-1.d0)/(n+1.d0)
  !cc        c2 = (n-2.d0)/(n+1.d0)
  !cc        Poly(n+1) = c1*x*Poly(n)-c2*t*t*Poly(n-1)
  !cc        if (iprint.eq.2) then
  !cc          write(*,*) n1,n2
  !cc          write(*,*) Poly(n),x,c1
  !cc          write(*,*) Poly(n-1),t*t,c2
  !cc          call pause
  !cc        endif
  !cc        Dpoly(1,n+1) = c1*(Poly(n)+x*Dpoly(1,n))-c2*t*t*Dpoly(1,n-1)
  !cc        Dpoly(2,n+1) = c1*x*Dpoly(2,n)
  !cc     .                    -c2*( 2.d0*t*Poly(n-1)
  !cc     .                          +t*t*Dpoly(2,n-1) )
  !cc      enddo
  !
  !
  !c ....without recursion formula - not tested
  !      n=0
  !      do ied=1,3
  !        ied1 = i3(ied)
  !        x    = shapH(ied)-shapH(ied1)
  !        t    = shapH(ied)+shapH(ied1)
  !c
  !c -----------------------------------------------------------
  !c L^s_2 = 0.5 (x^2-t^2)
  !        Ubb(n,ied)    = 0.5d0*x*x-0.5d0*t*t
  !        DUbb(1,n,ied) = x*(DshapH(1,ied)-DshapH(1,ied1))
  !     .                - t*(DshapH(1,ied)+DshapH(1,ied1))
  !        DUbb(2,n,ied) = x*(DshapH(2,ied)-DshapH(2,ied1))
  !     .                - t*(DshapH(2,ied)+DshapH(2,ied1))
  !        n = n+1
  !        if(n.eq.Nord) return
  !c
  !c -----------------------------------------------------------
  !c L^s_3 = 0.5 x (x^2-t^2)
  !        Ubb(n,ied)    = 0.5d0*x*x*x-0.5d0*x*t*t
  !        DUbb(1,n,ied) =
  !     .   0.5d0*(3.d0*x*x-t*t)*(DshapH(1,ied)-DshapH(1,ied1))
  !     .                  - x*t*(DshapH(1,ied)+DshapH(1,ied1))
  !        DUbb(2,n,ied) =
  !     .   0.5d0*(3.d0*x*x-t*t)*(DshapH(2,ied)-DshapH(2,ied1))
  !     .                  - x*t*(DshapH(2,ied)+DshapH(2,ied1))
  !        n = n+1
  !        if(n.eq.Nord) return
  !c
  !c -----------------------------------------------------------
  !c L^s_4 = 1/8 (5 x^4 - 6 x^2 t^2 + t^4)
  !        Ubb(n,ied)    = 5.d0/8.d0*x*x*x*x-3.d0/4.d0*x*x*t*t
  !     .                 +1.d0/8.d0*t*t*t*t
  !        DUbb(1,n,ied) = (5.d0/4.d0*x*x*x-3.d0/2.d0*x*t*t)*
  !     .                         (DshapH(1,ied)-DshapH(1,ied1))
  !     .                 +(-3.d0/4.d0*x*x*t+4.d0*t*t*t)*
  !     .                         (DshapH(1,ied)+DshapH(1,ied1))
  !        DUbb(2,n,ied) = (5.d0/4.d0*x*x*x-3.d0/2.d0*x*t*t)*
  !     .                         (DshapH(2,ied)-DshapH(2,ied1))
  !     .                 +(-3.d0/4.d0*x*x*t+4.d0*t*t*t)*
  !     .                         (DshapH(2,ied)+DshapH(2,ied1))
  !        n = n+1
  !        if(n.eq.Nord) return
  !c
  !c -----------------------------------------------------------
  !c L^s_5 = 1/8 (7 x^5 - s10 x^3 t^2 + 3 x t^4)
  !        Ubb(n,ied)    = 1.d0/8.d0*( 7.d0*x*x*x*x*x
  !     .                            -10.d0*x*x*x*t*t
  !     .                            + 3.d0*x*t*t*t*t)
  !        DUbb(1,n,ied) = (35.d0/8.d0*x*x*x*x
  !     .                  -30.d0/8.d0*x*x*t*t
  !     .                  + 3.d0/8.d0*t*t*t*t)*
  !     .                         (DshapH(1,ied)-DshapH(1,ied1))
  !     .                 +(5.d0/4.d0*x*x*x*t+3.d0/4.d0*x*t*t*t)*
  !     .                         (DshapH(1,ied)+DshapH(1,ied1))
  !        DUbb(2,n,ied) = (35.d0/8.d0*x*x*x*x
  !     .                  -30.d0/8.d0*x*x*t*t
  !     .                  + 3.d0/8.d0*t*t*t*t)*
  !     .                         (DshapH(2,ied)-DshapH(2,ied1))
  !     .                 +(5.d0/4.d0*x*x*x*t+3.d0/4.d0*x*t*t*t)*
  !     .                         (DshapH(2,ied)+DshapH(2,ied1))
  !        n = n+1
  !        if(n.eq.Nord) return
  !c
  !      enddo
  !
end subroutine legendre_S_I
!
!
!----------------------------------------------------------------------
!
!   routine name       - legendre_S
!
!----------------------------------------------------------------------
!
!   latest revision    - Oct 10
!
!   purpose            - evaluates scaled Legendre polynomials and
!                        their derivatives up to degree Nord
!
!   arguments :
!     in:
!                  X,T - coordinates
!                 Nord - degree
!     out:
!                 Poly - polynomials up to degree Nord
!                Dpoly - derivatives
!
!----------------------------------------------------------------------
!
subroutine legendre_S(X,T,Nord, Poly,Dpoly)
  !
  use parameters
#include "syscom.blk"
  !
  !  ...exceed maximum order of approximation
  integer, parameter :: ndelta = 2
  !
  !  ...output
  dimension Poly(0:MAXP+ndelta),Dpoly(2,0:MAXP+ndelta)
  !
  !  ...local variables
  !dimension aux(0:MAXP+ndelta),daux(0:MAXP+ndelta)
  !
  !----------------------------------------------------------------------
  !
  !  ...definition:
  !     l^{S}_n(x,t) := t^n*l_n(x/t)
  !
  !  ...2-terms recurrence relation from Sabine's dissertation
  !           l^{S}_0     = 1
  !           l^{S}_1     = x
  !     (n+1)*l^{S}_{n+1} = (2n+1)*x*l^{S}_n - n*t^2*l^{S}_{n-1}
  !
  !----------------------------------------------------------------------
  !
  if (Nord.gt.MAXP+ndelta) then
     write(*,*)'legendre_S: maximum order of approx. exceeded!'
     write(*,*)'Nord = ',Nord
     stop
  endif
  !
  !  ...l^{S}_0, \nabla l^{S}_0
  Poly(   0) = 1.d0
  Dpoly(1,0) = 0.d0
  Dpoly(2,0) = 0.d0
  !  ...l^{S}_1, \nabla l^{S}_1
  Poly(   1) = X
  Dpoly(1,1) = 1.d0
  Dpoly(2,1) = 0.d0
  !
  !  ...2-terms recurrence relation
  do n=1,Nord-1
     c1 = (2.d0*n+1.d0)/(n+1.d0)
     c2 = n/(n+1.d0)
     Poly(   n+1) = c1*X*Poly(n) - c2*T**2*Poly(n-1)
     Dpoly(1,n+1) = c1*(Poly(n) + X*Dpoly(1,n)) - &
          c2*T**2*Dpoly(1,n-1)
     Dpoly(2,n+1) = c1*X*Dpoly(2,n) -             &
          c2*(2.d0*T*Poly(n-1) + T**2*Dpoly(2,n-1))
  enddo
  !
  !
  !ccc  ...scaled variable
  !cc      s = X/T
  !ccc
  !ccc  ...compute Legendre polynomials
  !cc      call Legendre(s,Nord, aux,daux)
  !ccc  ...initialize
  !cc      Poly = 0.d0; Dpoly = 0.d0; Poly(0) = 1.d0; Dpoly(1,1) = 1.d0
  !ccc
  !cc      do i=1,Nord
  !cc        Poly(i) = T**i*aux(i)
  !cc
  !ccc  .....gradient
  !ccc       dl^S_{n}/dx = t^(n-1) l_{n}'
  !ccc       dl^S_{n}/dt = n t^(n-1) l_{n} - x t^(n-2) l_{n}'
  !ccc
  !cc        if (i.ge.2) then
  !cc          Dpoly(1,i) = T**(i-1)*daux(i)
  !cc          Dpoly(2,i) = i*T**(i-1)*aux(i) - X*T**(i-2)*daux(i)
  !cc        endif
  !cc      enddo
  !
  !
end subroutine legendre_S

!----------------------------------------------------------------------
!
!   routine name       - legendre
!
!----------------------------------------------------------------------
!
!   latest revision    - Oct 10
!
!   purpose            - evaluates Legendre polynomials and
!                        their derivatives up to degree Nord
!
!   arguments :
!     in:
!                    X - coordinate
!                 Nord - degree
!     out:
!                 Poly - polynomials up to degree Nord
!                Dpoly - derivatives
!
!----------------------------------------------------------------------
!
subroutine legendre_I(X,Nord, Poly,Dpoly)
  !
  use parameters
#include "syscom.blk"
  !
  !  ...exceed maximum order of approximtion
  integer, parameter :: ndelta = 2
  !
  !  ...output
  dimension Poly(0:MAXP+ndelta),Dpoly(0:MAXP+ndelta)
  !
  if (Nord.gt.MAXP+ndelta) then
     write(*,*)'legendre_I: maximum order of approx. exceeded!'
     write(*,*)'Nord = ',Nord
     stop
  endif
  !
  !  ...l_j, j=0,1
  Poly(0) = 0.d0;               Dpoly(0) = 0.d0
  Poly(1) = X;                  Dpoly(1) = 1.d0
  Poly(2) = 0.5d0*(X*X-1.d0);   Dpoly(2) = X
  !
  !  ...l_j, j>=2
  do n=2,Nord-1
     c1 = (2.d0*n-1.d0)/(n+1.d0)
     c2 = (n-2.d0)/(n+1.d0)
     Poly (n+1) = c1*x*Poly(n)-c2*Poly(n-1)
     Dpoly(n+1) = c1*(Poly(n) + X*Dpoly(n)) - c2*Dpoly(n-1)
  enddo
end subroutine legendre_I

!----------------------------------------------------------------------
!
!   routine name       - legendre
!
!----------------------------------------------------------------------
!
!   latest revision    - Oct 10
!
!   purpose            - evaluates Legendre polynomials and
!                        their derivatives up to degree Nord
!
!   arguments :
!     in:
!                    X - coordinate
!                 Nord - degree
!     out:
!                 Poly - polynomials up to degree Nord
!                Dpoly - derivatives
!
!----------------------------------------------------------------------
!
subroutine legendre(X,Nord, Poly,Dpoly)
  !
  use parameters
#include "syscom.blk"
  !
  !  ...exceed maximum order of approximtion
  integer, parameter :: ndelta = 2
  !
  !  ...output
  dimension Poly(0:MAXP+ndelta),Dpoly(0:MAXP+ndelta)
  !
  !  ...local variables
  !dimension sleg(MAXP+ndelta+1),dsleg(MAXP+ndelta+1)
  !
  !----------------------------------------------------------------------
  !  ...from LUCY'S routine..............................................
  !  ...Legendre polynomials
  !     l_0 = 1                                         ;   n == 0
  !     l_1 = x                                         ;   n == 1
  !     l_{n+1} = (2n+1)/(n+1) x l_{n} - n/(n+1) l_{n-1};   n >= 1
  !
  !  ...recurrence formula for computing derivatives
  !     dl_0/dx = 0
  !     dl_1/dx = 1
  !     dl_{n+1}/dx = (2n+1)/(n+1)( l_n + dl_n/dx )- n/(n+1) dl_{n-1}
  !----------------------------------------------------------------------
  !
  if (Nord.gt.MAXP+ndelta) then
     write(*,*)'legendre: maximum order of approx. exceeded!'
     write(*,*)'Nord = ',Nord
     stop
  endif
  !
  !  ...l_j, j=0,1
  Poly( 0) = 1.d0
  Dpoly(0) = 0.d0
  Poly( 1) = X
  Dpoly(1) = 1.d0
  !
  !cc      sleg(1) = 1.d0
  !cc      sleg(2) = x
  !cc      dsleg(1)= 0.d0
  !cc      dsleg(2)= 1.d0
  !
  !  ...l_j, j>=2
  do n=1,Nord-1
     c1 = (2.d0*n+1.d0)/(n+1.d0)
     c2 = n/(n+1.d0)
     Poly (n+1) = c1*x*Poly(n)-c2*Poly(n-1)
     Dpoly(n+1) = c1*(Poly(n)+x*Dpoly(n))-c2*Dpoly(n-1)
  enddo

  !cc      do n=2,Nord
  !cc        c1 = (2.d0*n+1.d0)/(n+1.d0)
  !cc        c2 = n/(n+1.d0)
  !cc        sleg (n+1) = c1*x*sleg(n)-c2*sleg(n-1)
  !cc        dsleg(n+1) = c1*(sleg(n)+x*dsleg(n))-c2*dsleg(n-1)
  !cc      enddo
  !
  !  ...copy to output
  !cc      Poly(0:MAXP+ndelta)  = sleg(1:MAXP+ndelta+1)
  !cc      Dpoly(0:MAXP+ndelta) = dsleg(1:MAXP+ndelta+1)
  !
  !
  !c
  !c -----------------------------------------------------------
  !c ....without recursion formula - not tested
  !      n=0
  !      do ino=1,3
  !        y=ShapH
  !c -----------------------------------------------------------
  !c l_0(x) = 1  => v_1 = shapH
  !        Vbb(n,ino) = y
  !        DVbb(1,n) = DshapH(1)
  !        DVbb(2,n) = DshapH(2)
  !        n = n+1
  !        if(n.eq.Nord) return
  !c -----------------------------------------------------------
  !c l_1(x) = x  => v_1 = shapH*(2 shapH -1)
  !        Vbb(n,ino) = y*(2.d0*y-1.d0)
  !        DVbb(1,n) = 4.d0*DshapH(1)
  !        DVbb(2,n) = 4.d0*DshapH(2)
  !        n = n+1
  !        if(n.eq.Nord) return
  !c -----------------------------------------------------------
  !c l_2(x) = 1/2 (3x^2-1)  => v_2 =shapH  1/2 (3 (2 shapH -1)^2 -1)
  !        Vbb(n) = y*(8.d0*y*y-2.d0*y+1.d0)
  !        DVbb(1,n) = (24.d0*y*y-6.d0*y+1.d0)*DshapH(1)
  !        DVbb(2,n) = (24.d0*y*y-6.d0*y+1.d0)*DshapH(2)
  !c -----------------------------------------------------------
  !c l_3(x) = 1/2 x (5x^2-3) => v_4 =shapH 1/2 (2 shapH -1) (5(2 shapH -1)^2-3)
  !        Vbb(n) = 20.d0*y*y*y*y-30.d0*y*y*y+12.d0*y*y-y
  !        DVbb(1,n) = (80.d0*y*y*y-90.d0*y*y+24.d0+*y-1.d0)
  !     .                  *DshapH(1)
  !        DVbb(2,n) = (80.d0*y*y*y-90.d0*y*y+24.d0+*y-1.d0)
  !     .                  *DshapH(2)
  !      enddo
  !c
end subroutine legendre


!--------------------------------------------------------------
!  Apr 2009
!
!  ROUTINES FOR COMPUTING LEGENDRE POLYNOMIALS DOWNLOADED FROM:
!
!  http://jin.ece.uiuc.edu/routines/routines.html
!--------------------------------------------------------------
!
!
!======================================================
SUBROUTINE LPN(N,X,PN,PD)
  !======================================================
  !       Purpose: Compute Legendre polynomials Pn(x)
  !                and their derivatives Pn'(x)
  !       Input :  x --- Argument of Pn(x)
  !                n --- Degree of Pn(x) ( n = 0,1,...)
  !       Output:  PN(n) --- Pn(x)
  !                PD(n) --- Pn'(x)
  !       ===============================================
  !
  IMPLICIT DOUBLE PRECISION (P,X)
  DIMENSION PN(0:N+1),PD(0:N+1)
  PN(0)=1.0D0
  PN(1)=X
  PD(0)=0.0D0
  PD(1)=1.0D0
  P0=1.0D0
  P1=X
  DO K=2,N
     PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
     PN(K)=PF
     IF (DABS(X).EQ.1.0D0) THEN
        PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
     ELSE
        PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
     ENDIF
     P0=P1
     P1=PF
  END DO
END SUBROUTINE LPN
