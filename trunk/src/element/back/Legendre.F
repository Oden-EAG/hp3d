c----------------------------------------------------------------------
c
c     routine name      - Legendre
c
c----------------------------------------------------------------------
c
c     latest revision:  - Sep 13
c
c     purpose:          - routine returns values of shifted scaled
c                         Legendre polynomials and their derivatives
c                         wrt the scaling parameter
c
c     arguments:
c
c     in:
c             X         - coordinate from [0,1]
c             T         - scaling parameter
c             Nord      - polynomial order
c
c     out:
c             P         - polynomial values
c
c-----------------------------------------------------------------------
c
      subroutine LegendreNEW(X,T,Nord, P)
c
      use parametersDPG
      implicit none
      integer          :: Nord,Idec, i,iprint
      double precision :: X,T, tt,y
      double precision, dimension(0:MAXP) :: P
c
      iprint=0
c
c  ...shift argument
      y = 2.d0*X - T
c
c  ...use recurrence formula
      P(0) = 1.d0
      P(1) = y
c
c  ...i stands for the order of the polynomial, stored in P(i)
      tt = T**2
      do i=2,Nord
        P(i) = (2*i-1)*y*P(i-1) - (i-1)*tt*P(i-2)
        P(i) = P(i)/i
      enddo
c
      if (iprint.eq.1) then
        write(*,7001) Idec,Nord, X,T
 7001   format('Legendre: Idec = ',i1,' Nord = ',i2,' X,T = ',2F8.3)
        do i=0,Nord
          write(*,7002) i,P(i)
 7002     format('i = ',i2,' P = ',e25.15)
        enddo
        call pause
      endif
c
c
      end subroutine LegendreNEW
c
c--------------------------------------------------------------------
c
c     routine name      - ILegendre
c
c--------------------------------------------------------------------
c
c     latest revision:  - Sep 13
c
c     purpose:          - routine returns values of shifted scaled
c                         integrated polynomials and their derivatives
c                         starting with p=2
c
c     arguments:
c
c     in:
c             X         - coordinate from [0,1]
c             T         - scaling parameter
c             Nord      - polynomial order
c             Idec      - decision flag to compute:
c                       = 0  integrated Legendre polynomials with x derivatives
c                       = 1  integrated polynomials with t derivatives
c
c     out:
c             RL        - polynomial values
c             DLdx      - derivatives in x
c             DLdt      - derivatives in t
c
c-----------------------------------------------------------------------
c
      subroutine ILegendre(X,T,Nord,Idec, RL,DLdx,DLdt)
c
      use parameters
      implicit none
      integer          :: Nord,Idec, iprint,i,ifact
      double precision :: X,T, tt,y
      double precision, dimension(2:MAXP)   :: RL,DLdx,DLdt
      double precision, dimension(0:MAXP)   :: p
c  ...work space for Legendre
c
      iprint=0
      if (Nord.lt.2) return
c
      call LegendreNEW(X,T,Nord, p)
c
c  ...Integrated Legendre polynomial of order i is stored in RL(i)
      tt = T**2
      do i=2,Nord
        ifact = 4*i-2
        RL(i) = (p(i) - tt*p(i-2))/ifact
        DLdx(i) = p(i-1)
        if (Idec.ne.0) then
          DLdt(i) = -(p(i-1)+T*p(i-2))/2
        endif
      enddo
c
c
      if (iprint.eq.1) then
        write(*,7001) Idec,Nord, X,T
 7001   format('ILegendre: Idec = ',i1,' Nord = ',i2,' X,T = ',2F8.3)
        do i=2,Nord
          select case(Idec)
          case(0)
            write(*,7002) i,RL(i),DLdx(i)
 7002       format('i = ',i2,' L,DLdx = ',2e25.15)
          case default
            write(*,7003) i,RL(i),DLdx(i),DLdt(i)
 7003       format('i = ',i2,' L,DLdx,DLdt = ',3e25.15)
          end select
        enddo
        call pause
      endif
c
      end subroutine ILegendre

c
c----------------------------------------------------------------------
c
c     routine name      - Jacobi
c
c----------------------------------------------------------------------
c
c     latest revision:  - Sep 14
c
c     purpose:          - routine returns values of shifted scaled
c                         Jacobi polynomials P^\alpha_i and their
c                         derivatives wrt the scaling parameter.
c                         P^\alpha_i With alpha = 2*a+ralpha of
c                         order i is stored in P(a,i).
c
c     arguments:
c
c     in:
c             X         - coordinate from [0,1]
c             T         - scaling parameter
c             Maxalpha  - max value of alpha (integer)
c             Nord      - max polynomial order
c
c     out:
c             P         - polynomial values
c
c-----------------------------------------------------------------------
c
      subroutine Jacobi(X,T,Maxalpha,Nord, P)
c
      use parameters
      implicit none
      integer :: Maxalpha,Nord, i,a,alpha,aa,ai,bi,ci,di,iprint,
     .           nalpha,ralpha
      double precision :: X,T, y,tt
      double precision, dimension(0:MAXP+1,0:MAXP) :: P

      iprint=0
      if (Maxalpha.lt.1) then
        write(*,7001) Maxalpha
 7001   format('Jacobi: Maxalpha = ',i3)
        stop 1
      endif
      if (Nord.gt.MAXP) then
        write(*,7002) Nord, MAXP
 7002   format('Jacobi: Nord = ',i3,' MAXP = ',i3)
        stop 1
      endif
c
c  ...Jacobi polynomial of weight 2*a+ralpha of order i is stored in P(a,i)
      ralpha = mod(Maxalpha,2)
      nalpha = (Maxalpha-ralpha)/2
c
      if ((nalpha.lt.1).or.(nalpha.gt.MAXP+1)) then
        write(*,*) 'Jacobi: nalpha,MAXP+1 = ',nalpha,MAXP+1
        stop 1
      endif
c
c  ...initiate zero order polynomials for all weights
      P(0:nalpha,0) = 1.d0
      tt = T**2
      y  = 2.d0*X - T
      do a=0,nalpha
        alpha = 2*a+ralpha
        aa    = alpha**2
c
        P(a,1) = y+alpha*X
c
c    ...use recursion in order, i, to compute P^alpha_i for i>=2
        if (alpha.ne.0) then
          do i=2,Nord
            ai = 2*i*(i+alpha)*(2*i+alpha-2)
            bi = 2*i+alpha-1
            ci = (2*i+alpha)*(2*i+alpha-2)
            di = 2*(i+alpha-1)*(i-1)*(2*i+alpha)
c
            P(a,i) = bi*(ci*y +aa*T)*P(a,i-1)-di*tt*P(a,i-2)
            P(a,i) = P(a,i)/ai
          enddo
c
c    ...simplified computation
        else
          do i=2,Nord
            P(0,i) = (2*i-1)*y*P(0,i-1)-(i-1)*tt*P(0,i-2)
            P(0,i) = P(0,i)/i
          enddo
        endif
c
      enddo
c
c
 999  continue
      if (iprint.eq.1) then
        write(*,7003) Nord, X,T
 7003   format('Jacobi: Nord = ',i2,' X,T = ',2F8.3)
        do alpha=Maxalpha,1,-1
          write(*,7004) alpha,P(alpha,0:Nord)
 7004     format(' P(',i2,',0:Nord) = ',10e12.5)
        enddo
        call pause
      endif
c
      end subroutine Jacobi
c
c----------------------------------------------------------------------
c
c     routine name      - Jacobia_temp
c
c----------------------------------------------------------------------
c
c     latest revision:  - Mar 14
c
c     purpose:          - routine returns values of shifted scaled
c                         Jacobi polynomials P^\alpha_n and their
c                         derivatives wrt the scaling parameter
c
c     arguments:
c
c     in:
c             Xi        - coordinate from [0,1]
c             T         - scaling parameter
c             Maxalpha  - max weight
c             Nord      - max polynomial order
c             Idec      - decision flag to compute:
c                       = 0  Jacobi polynomials
c                       = 1  Jacobi polynomials with derivatives in t
c
c     out:
c             P         - polynomial values
c             DPdt      - derivatives in t
c
c-----------------------------------------------------------------------
c
      subroutine Jacobia_temp(Xi,T,Maxalpha,Nord,Idec, P,DPdt)
c
      use parameters
#include "syscom.blk"
c
      dimension P(1:3*MAXP,0:MAXP),DPdt(1:3*MAXP,0:MAXP)
      integer alpha
c
      iprint=0
      if ((Maxalpha.lt.1).or.(Maxalpha.gt.3*MAXP)) then
        write(*,7001) Maxalpha,3*MAXP
 7001   format('Jacobia_temp: Maxalpha,3*MAXP = ',2i3)
        stop 1
      endif
c
c  ...shift argument
      x = 2.d0*Xi - T
c
c  ...initiate zero order polynomials for all weights
      P(1:Maxalpha,0) = 1.d0
      if (Idec.eq.1) DPdt(1:Maxalpha,0) =  0.d0
c
c  ...use recursion in degree n to compute P^1_n
      alpha=1
      P(1,1) = ((x-T) + (alpha+1)*(x+T))/2.d0
      if (Idec.eq.1) DPdt(1,1) = -1.d0
      naux = alpha**2
      do n=2,Nord
        na = 2*n*(n+alpha)*(2*n+alpha-2)
        nb = 2*n+alpha-1
        nc = (2*n+alpha)*(2*n+alpha-2)
        nd = 2*(n+alpha-1)*(n-1)*(2*n+alpha)
        P(1,n) = nb*(nc*x +naux*T)*P(1,n-1) - nd*T**2*P(1,n-2)
        P(1,n) = P(1,n)/na
        if (Idec.eq.1) then
          DPdt(1,n) = nb*(-nc+naux)*P(1,n-1)
     .            + nb*(nc*x +naux*T)*DPdt(1,n-1)
     .            - 2*nd*T*P(1,n-2)
     .            - nd*T**2*DPdt(1,n-2)
          DPdt(1,n) = DPdt(1,n)/na
        endif
      enddo
c
c  ...use recursion in weight and degree to compute the remaining
c     polynomials
      do n=1,Nord
        do alpha=2,Maxalpha
          P(alpha,n) = (alpha+2*n)*P(alpha-1,n) + n*T*P(alpha,n-1)
          P(alpha,n) = P(alpha,n)/(alpha+n)
          if (Idec.eq.1) then
            DPdt(alpha,n) = (alpha+2*n)*DPdt(alpha-1,n)
     .                    + n*P(alpha,n-1) + n*T*DPdt(alpha,n-1)
            DPdt(alpha,n) = DPdt(alpha,n)/(alpha+n)
          endif
        enddo
      enddo
c
c
 999  continue
      if (iprint.eq.1) then
        write(*,7002) Idec,Nord, X,T
 7002   format('Jacobia: Idec = ',i1,
     .         ' Nord = ',i2,' X,T = ',2F8.3)
        do alpha=Maxalpha,1,-1
          write(*,7003) alpha,P(alpha,0:Nord)
 7003     format(' P(',i2,',0:Nord) = ',10e12.5)
        enddo
        if (Idec.eq.1) then
          write(*,*) '  '
          do alpha=Maxalpha,1,-1
            write(*,7004) alpha,DPdt(alpha,0:n)
 7004       format(' DPdt(',i2,',0:Nord) = ',10e12.5)
          enddo
        endif
        call pause
      endif
c
      end subroutine Jacobia_temp
c
c--------------------------------------------------------------------
c
c     routine name      - IJacobi
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of integrated shifted
c                         scaled Jacobi polynomials and their derivatives
c                         starting with p=1
c
c     arguments:
c
c     in:
c             X         - coordinate from [0,1]
c             T         - scaling parameter
c             Imax      - range for parameter i specifying the weight
c                         alpha = 2i-1
c             Nord      - max polynomial order in x
c             Idec      - decision flag to compute:
c                       = 0  polynomials with x derivatives
c                       = 1  polynomials with t derivatives
c
c     out:
c             RL        - polynomial values
c             DLdx      - derivatives in x
c             DLdt      - derivatives in t
c
c-----------------------------------------------------------------------
c
      subroutine IJacobi(X,T,Imax,Nord,Idec, RL,DLdx,DLdt)
c
      use parameters
#include "syscom.blk"
c
      dimension RL(1:MAXP+1,1:MAXP),DLdx(1:MAXP+1,1:MAXP),
     .          DLdt(1:MAXP+1,1:MAXP)
ccc  P. Gatto, Nov 24, 2014
ccc      dimension RL(1:MAXP+1,1:MAXP-1),DLdx(1:MAXP+1,1:MAXP-1),
ccc     .          DLdt(1:MAXP+1,1:MAXP-1)
      integer alpha
c
c  ...work space for Jacobi
      dimension p(1:3*MAXP,0:MAXP),dPdt(1:3*MAXP,0:MAXP)
c
c
      iprint=0
      if ((Nord.lt.1).or.(Nord.gt.MAXP)) then
        write(*,*) 'IJacobi: Nord,MAXP = ', Nord,MAXP
        stop 1
      endif

      if ((Imax.lt.1).or.(Imax.gt.MAXP+1)) then
        write(*,*) 'IJacobi: Imax,MAXP+1 = ',Imax,MAXP+1
        stop 1
      endif
c
      call Jacobia_temp(X,T,2*Imax-1,Nord,Idec, p,dPdt)
c
      RL(1:Imax,1) = X; DLdx(1:Imax,1) = 1.d0
      if (Idec.ne.0) DLdt(1:Imax,1) = 0.d0
c
      do n=2,Nord
        do i=1,Imax
c
c  .......set the weight in proportion to order i
          alpha = 2*i-1
          talpha = alpha+alpha
          tn = n+n
          tna = tn+alpha
          tnam1 = tna-1.d0
          tnam2 = tna-2.d0
          an = (tn+talpha)/(tnam1*tna)
          bn = talpha/(tnam2*tna)
          cn = (tn-2.d0)/(tnam2*tnam1)
          RL(i,n) = an*P(alpha,n)+bn*T*P(alpha,n-1)-cn*T**2*P(alpha,n-2)
          RL(i,n) = RL(i,n)/2.d0
          DLdx(i,n) =  P(alpha,n-1)
          if (Idec.ne.0) then
            DLdt(i,n) = an*dPdt(alpha,n)
     .                + bn*(P(alpha,n-1)+T*dPdt(alpha,n-1))
     .                - cn*T*(2.d0*P(alpha,n-2)+T*dPdt(alpha,n-2))
            DLdt(i,n) = DLdt(i,n)/2.d0
          endif
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,7001) Idec,Imax,Nord, X,T
 7001   format('IJacobi: Idec = ',i1,' Imax = ',i2,
     .         ' Nord = ',i2,' X,T = ',2F8.3)
        do i=1,Imax
          alpha = 2*i-1
          write(*,7002) i,alpha,RL(i,1:Nord)
 7002     format('i = ',i1,' alpha = ',i2,' RL(i,1:Nord)   = ',10e12.5)
        enddo
        write(*,*) '  '
        do i=1,Imax
          alpha = 2*i-1
          write(*,7003) i,alpha,DLdx(i,1:Nord)
 7003     format('i = ',i1,' alpha = ',i2,' DLdx(1:Nord) = ',10e12.5)
        enddo
        if (Idec.ne.0) then
          write(*,*) '  '
          do i=1,Imax
            alpha = 2*i-1
            write(*,7004) i,alpha,DLdt(i,1:Nord)
 7004       format('i = ',i1,' alpha = ',i2,
     .             ' DLdt(1:Nord) = ',10e12.5)
          enddo
        endif
        call pause
      endif
c
      end subroutine IJacobi


      subroutine my_tests
c
      use parameters
      use parametersDPG
c
 10   write(*,*) 'SELECT'
      write(*,*) 'EXIT........................................0'
c

      end subroutine my_tests
