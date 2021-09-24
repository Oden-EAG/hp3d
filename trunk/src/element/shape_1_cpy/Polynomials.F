c Routines:
c  - PolyLegendre
c  - PolyILegendre
c  - PolyJacobi
c  - PolyIJacobi
c  - HomLegendre
c  - HomILegendre
c  - HomJacobi
c  - HomIJacobi
c
c----------------------------------------------------------------------
c
c     routine name      - PolyLegendre
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of shifted scaled
c                         Legendre polynomials
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
c----------------------------------------------------------------------
c
      subroutine PolyLegendre(X,T,Nord, P)
c
      implicit none
      integer,          intent(in)  ::                             Nord
      double precision, intent(in)  ::                              X,T
      double precision, intent(out) ::                        P(0:Nord)
      integer ::                                               i,iprint
      double precision ::                                          tt,y
c
      iprint=0
c
c  ...i stands for the order of the polynomial, stored in P(i)
c  ...lowest order case (order 0)
      P(0) = 1.d0
c  ...first order case (order 1) if necessary
      if (Nord.ge.1) then
        y = 2.d0*X - T
        P(1) = y
      endif
c  ...higher order if necessary - use recurrence formula
      if (Nord.ge.2) then
        tt = T**2
        do i=2,Nord
          P(i) = (2*i-1)*y*P(i-1) - (i-1)*tt*P(i-2)
          P(i) = P(i)/i
        enddo
      endif
c
c  ...catching problems (debugging)
      if (iprint.eq.1) then
        write(*,7002) Nord, X,T
 7002   format('PolyLegendre: Nord = ',i2,' X,T = ',2F8.3)
        do i=0,Nord
          write(*,7003) i,P(i)
 7003     format('i = ',i2,' P = ',e25.15)
        enddo
        call pause
      endif
c
      end subroutine PolyLegendre
c
c----------------------------------------------------------------------
c
c     routine name      - PolyILegendre
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of shifted scaled
c                         integrated Legendre polynomials and their
c                         derivatives starting with p=2
c
c     arguments:
c
c     in:
c             X         - coordinate from [0,1]
c             T         - scaling parameter
c             Nord      - polynomial order
c             Idec      - decision flag to compute:
c                       = FALSE polynomials with x and t derivatives
c                       = TRUE  polynomials with x derivatives only
c
c     out:
c             L         - polynomial values
c             P         - derivatives in x
c             R         - derivatives in t
c
c----------------------------------------------------------------------
c
      subroutine PolyILegendre(X,T,Nord,Idec, L,P,R)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                             Nord
      double precision, intent(in)  ::                              X,T
      double precision, intent(out) ::    L(2:Nord),P(Nord-1),R(Nord-1)
      integer ::                                         i,ifact,iprint
      double precision ::                              ptemp(0:Nord),tt
c
      iprint=0
c
c  ...calling Legendre for required information
      call PolyLegendre(X,T,Nord, ptemp)
      P = ptemp(1:Nord-1)
c
c  ...Integrated polynomial of order i is stored in L(i)
      tt = T**2
c
c  ...simplified case: no need to compute R
      if (Idec) then
        do i=2,Nord
          ifact = 4*i-2
          L(i) = (ptemp(i) - tt*ptemp(i-2))/ifact
        enddo
c
c  ...general case: compute R
      else
        do i=2,Nord
          ifact = 4*i-2
          L(i) = (ptemp(i) - tt*ptemp(i-2))/ifact
          R(i-1) = -(ptemp(i-1)+T*ptemp(i-2))/2
        enddo
      endif
c
c  ...catching problems (debugging)
      if (iprint.eq.1) then
        write(*,7002) Idec,Nord, X,T
 7002   format('PolyILegendre: Idec = ',i1,' Nord = ',i2,
     .         ' X,T = ',2F8.3)
        do i=2,Nord
          select case(Idec)
          case(.TRUE.)
            write(*,7003) i,L(i),P(i)
 7003       format('i = ',i2,' L,P = ',2e25.15)
          case default
            write(*,7004) i,L(i),P(i),R(i)
 7004       format('i = ',i2,' L,P,R, = ',3e25.15)
          end select
        enddo
        call pause
      endif
c
      end subroutine PolyILegendre
c
c
c----------------------------------------------------------------------
c
c     routine name      - PolyJacobi
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of shifted scaled
c                         Jacobi polynomials P^\alpha_i. Result is a
c                         'half' of a  matrix with each row
c                         associated to a fixed alpha. Alpha grows
c                         by 2 in each row.
c
c     arguments:
c
c     in:
c             X         - coordinate from [0,1]
c             T         - scaling parameter
c             Nord      - max polynomial order
c             Minalpha  - first row value of alpha (integer)
c
c     out:
c             P         - polynomial values
c
c----------------------------------------------------------------------
c
      subroutine PolyJacobi(X,T,Nord,Minalpha, P)
c
      implicit none
      integer,          intent(in)  ::                    Nord,Minalpha
      double precision, intent(in)  ::                              X,T
      double precision, intent(out) ::                 P(0:Nord,0:Nord)
      integer :: minI,maxI,i,ni,a,aa,al,ai,bi,ci,di,iprint,alpha(0:Nord)
      double precision ::                                          y,tt
c
c  ...clearly (minI,maxI)=(0,Nord), but the syntax is written as it is
c     because it reflects how the indexing is called from outside
      minI = 0; maxI = minI+Nord
c
c  ...in our work Minalpha>=1
      iprint = 0
      if (Minalpha.lt.1) then
        write(*,7001) Minalpha
 7001   format('PolyJacobi: Minalpha = ',i3)
        stop 1
      endif
c
c  ...create vector alpha first
      do a=minI,maxI
         alpha(a) = Minalpha+2*(a-minI)
      enddo
c
c  ...initiate first column (order 0)
      P(minI:maxI,0) = 1.d0
c  ...initiate second column (order 1) if necessary
      if (Nord.ge.1) then
        y = 2.d0*X - T
        P(minI:maxI-1,1) = y+alpha(minI:maxI-1)*X
      endif
c  ...fill the last columns if necessary
      if (Nord.ge.2) then
        tt = T**2
        ni = -1
        do a=minI,maxI-2
          al=alpha(a)
          aa = al**2
          ni=ni+1
c      ...use recursion in order, i, to compute P^alpha_i for i>=2
          do i=2,Nord-ni
            ai = 2*i*(i+al)*(2*i+al-2)
            bi = 2*i+al-1
            ci = (2*i+al)*(2*i+al-2)
            di = 2*(i+al-1)*(i-1)*(2*i+al)
c
            P(a,i) = bi*(ci*y+aa*T)*P(a,i-1)-di*tt*P(a,i-2)
            P(a,i) = P(a,i)/ai
          enddo
        enddo
      endif
c
c  ...catching problems (debugging)
      if (iprint.eq.1) then
        write(*,7003) Nord, X,T
 7003   format('PolyJacobi: Nord = ',i2,' X,T = ',2F8.3)
        do a=minI,maxI
          write(*,7004) a,P(a,0:Nord)
 7004     format(' P(',i2,',0:Nord) = ',10e12.5)
        enddo
        call pause
      endif
c
      end subroutine PolyJacobi
c
c----------------------------------------------------------------------
c
c     routine name      - PolyIJacobi
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of integrated
c                         shifted scaled Jacobi polynomials and
c                         their derivatives starting with p=1
c                         Result is 'half' of a  matrix
c                         with each row  associated to a fixed alpha.
c                         Alpha grows by 2 in each row.
c
c     arguments:
c
c     in:
c             X         - coordinate from [0,1]
c             T         - scaling parameter
c             Nord      - max polynomial order
c             Minalpha  - first row value of alpha (integer)
c             Idec      - decision flag to compute:
c                       = FALSE polynomials with x and t derivatives
c                       = TRUE  polynomials with x derivatives only
c
c     out:
c             L        - polynomial values
c             P        - derivatives in x (Jacobi polynomials)
c             R        - derivatives in t
c
c----------------------------------------------------------------------
c
      subroutine PolyIJacobi(X,T,Nord,Minalpha,Idec, L,P,R)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                    Nord,Minalpha
      double precision, intent(in)  ::                              X,T
      double precision, intent(out) ::                 L(1:Nord,1:Nord),
     .                            P(1:Nord,0:Nord-1),R(1:Nord,0:Nord-1)
      integer ::             minI,maxI,i,ni,a,al,tia,tiam1,tiam2,iprint,
     .                                                    alpha(1:Nord)
      double precision ::            ai,bi,ci,tt,ptemp(1:Nord+1,0:Nord)
c
      iprint=0
c
c  ...clearly (minI,maxI)=(1,Nord), but the syntax is written as it is
c     because it reflects how the indexing is called from outside
      minI = 1; maxI = minI+Nord-1
c
c  ...in our work Minalpha>=1
      if (Minalpha.lt.1) then
        write(*,7001) Minalpha
 7001   format('PolyIJacobi: Minalpha = ',i3)
        stop 1
      endif
c
c  ...calling Jacobi for required information
      call PolyJacobi(X,T,Nord,Minalpha, ptemp)
c  ...define P. Note that even though P is defined at all entries,
c     because of the way Jacobi computes ptemp, only the necessary entries,
c     and those on the first subdiagonal (which are never used later)
c     are actually accurate.
      P = ptemp(minI:maxI,0:Nord-1)
c
c  ...create vector alpha first
      do a=minI,maxI
         alpha(a) = Minalpha+2*(a-minI)
      enddo
c
c  ...initiate first column (order 1 in L)
      L(minI:maxI,1) = X
c
c  ...simplified case, do not compute R
      if (Idec) then
c  ...fill the last columns if necessary
        if (Nord.ge.2) then
          tt = T**2
          ni = -1
          do a=minI,maxI-1
            al=alpha(a)
            ni = ni+1
            do i=2,Nord-ni
              tia = i+i+al
              tiam1 = tia-1
              tiam2 = tia-2
              ai = dble(i+al)/(tiam1*tia)
              bi = dble(al)/(tiam2*tia)
              ci = (i-1.d0)/(tiam2*tiam1)
c
              L(a,i) = ai*ptemp(a,i)+bi*T*ptemp(a,i-1)
     .              -ci*tt*ptemp(a,i-2)
c              P(a,i-1) =  ptemp(a,i-1)
            enddo
          enddo
        endif
c
c  ...general case; compute R
      else
      R(minI:maxI,0) = 0.d0
c  ...fill the last columns if necessary
        if (Nord.ge.2) then
          tt = T**2
          ni = -1
          do a=minI,maxI-1
            al=alpha(a)
            ni = ni+1
            do i=2,Nord-ni
              tia = i+i+al
              tiam1 = tia-1
              tiam2 = tia-2
              ai = dble(i+al)/(tiam1*tia)
              bi = dble(al)/(tiam2*tia)
              ci = (i-1.d0)/(tiam2*tiam1)
c
              L(a,i) = ai*ptemp(a,i)+bi*T*ptemp(a,i-1)
     .              -ci*tt*ptemp(a,i-2)
c              P(a,i-1) =  ptemp(a,i-1)
              R(a,i-1) = -(i-1)*(ptemp(a,i-1)+T*ptemp(a,i-2))
              R(a,i-1) = R(a,i-1)/tiam2
            enddo
          enddo
        endif
      endif
c
c  ...catching problems (debugging)
      if (iprint.eq.1) then
        write(*,7003) minI,Nord, X,T
 7003   format('PolyIJacobi: minI = ',i2,
     .         ' Nord = ',i2,' X,T = ',2F8.3)
        do a=minI,maxI
          al = alpha(a)
          write(*,7004) a,al,L(a,1:Nord)
 7004     format('a = ',i1,' alpha = ',i2,
     .         ' L(a,1:Nord)   = ',10e12.5)
        enddo
        write(*,*) '  '
        do a=minI,maxI
          al = alpha(a)
          write(*,7005) a,al,P(a,0:Nord-1)
 7005     format('a = ',i1,' alpha = ',i2,
     .         ' P(a,0:Nord-1) = ',10e12.5)
        enddo
        if (.NOT.Idec) then
          write(*,*) '  '
          do a=minI,maxI
            al = alpha(a)
            write(*,7006) a,al,R(a,0:Nord-1)
 7006       format('a = ',i1,' alpha = ',i2,
     .             ' R(a,0:Nord-1) = ',10e12.5)
          enddo
        endif
        call pause
      endif
c
      end subroutine PolyIJacobi
c
c----------------------------------------------------------------------
c
c     routine name      - HomLegendre
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of homogenized
c                         Legendre polynomials
c
c     arguments:
c
c     in:
c             S         - affine(like) coordinates
c             Nord      - polynomial order
c
c     out:
c             HomP      - polynomial values
c
c----------------------------------------------------------------------
c
      subroutine HomLegendre(S,Nord, HomP)
c
      implicit none
      integer,          intent(in)  ::                             Nord
      double precision, intent(in)  ::                           S(0:1)
      double precision, intent(out) ::                     HomP(0:Nord)
c
c  ...simply the definition of homogenized polynomials
      call PolyLegendre(S(1),S(0)+S(1),Nord, HomP)
c
      end subroutine HomLegendre
c
c
c----------------------------------------------------------------------
c
c     routine name      - HomILegendre
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of homogenized
c                         integrated Legendre polynomials and their
c                         gradient (wrt to affine like coordinates)
c
c     arguments:
c
c     in:
c             S         - (s0,s1) affine(like) coordinates
c             DS        - gradients of S (in R^N)
c             Nord      - polynomial order
c             Idec      - decision flag to compute:
c                         = FALSE s0+s1 != 1 -> general case
c                         = TRUE  s0+s1  = 1 -> simple case
c             N         - number of spatial dimensions (R^N)
c
c     out:
c             HomL        - polynomial values
c             DHomL       - gradients of L in R^N
c
c----------------------------------------------------------------------
c
      subroutine HomILegendre(S,DS,Nord,Idec,N, HomL,DHomL)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:1),DS(1:N,0:1)
      double precision, intent(out) ::   HomL(2:Nord),DHomL(1:N,2:Nord)
      integer ::                                                      i
      double precision ::       homP(1:Nord-1),homR(1:Nord-1),DS01(1:N)
c
c  ...Idec is the flag to compute x AND t derivatives
c  ...If sum of S equal 1 -> Idec=.TRUE.
      if (Idec) then
        call PolyILegendre(S(1),1.d0,Nord,Idec, HomL,homP,homR)
        do i=2,Nord
          DHomL(1:N,i) = homP(i-1)*DS(1:N,1)
        enddo
c
c  ...If sum of S different from 1 -> Idec=.FALSE.
      else
        call PolyILegendre(S(1),S(0)+S(1),Nord,Idec, HomL,homP,homR)
        DS01 = DS(1:N,0)+DS(1:N,1)
        do i=2,Nord
          DHomL(1:N,i) = homP(i-1)*DS(1:N,1)+homR(i-1)*DS01
        enddo
      endif
c
      end subroutine HomILegendre
c
c
c----------------------------------------------------------------------
c
c     routine name      - HomJacobi
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of homogenized
c                         Jacobi polynomials P^\alpha_i. Result is a
c                         'half' of a  matrix with each row
c                         associated to a fixed alpha. Alpha grows
c                         by 2 in each row.
c
c     arguments:
c
c     in:
c             S         - affine(like) coordinates
c             Nord      - max polynomial order
c             Minalpha  - first row value of alpha (integer)
c
c     out:
c             HomP      - polynomial values
c
c----------------------------------------------------------------------
c
      subroutine HomJacobi(S,Nord,Minalpha, HomP)
c
      implicit none
      integer,          intent(in)  ::                    Nord,Minalpha
      double precision, intent(in)  ::                           S(0:1)
      double precision, intent(out) ::              HomP(0:Nord,0:Nord)
      integer ::                                              minI,maxI
c
c  ...clearly (minI,maxI)=(0,Nord), but the syntax is written as it is
c     because it reflects how the indexing is called from outside
      minI = 0; MaxI = MinI+Nord
c
c  ...simply the definition of homogenized polynomials
      call PolyJacobi(S(1),S(0)+S(1),Nord,Minalpha, HomP)
c
      end subroutine HomJacobi
c
c
c----------------------------------------------------------------------
c
c     routine name      - HomIJacobi
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of integrated
c                         homogenized Jacobi polynomials and
c                         their gradients.
c                         Result is 'half' of a  matrix
c                         with each row  associated to a fixed alpha.
c                         Alpha grows by 2 in each row.
c
c     arguments:
c
c     in:
c             S         - (s0,s1) affine(like) coordinates
c             DS        - gradients of S (in R^N)
c             Nord      - max polynomial order
c             Minalpha  - first row value of alpha (integer)
c             Idec      - decision flag to compute:
c                         = FALSE s0+s1 != 1 -> general case
c                         = TRUE  s0+s1  = 1 -> simple case
c             N         - number of spatial dimensions (R^N)
c
c     out:
c             HomL      - polynomial values
c             DHomL     - derivatives in x (Jacobi polynomials)
c
c----------------------------------------------------------------------
c
      subroutine HomIJacobi(S,DS,Nord,Minalpha,Idec,N, HomL,DHomL)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                  Nord,Minalpha,N
      double precision, intent(in)  ::               S(0:1),DS(1:N,0:1)
      double precision, intent(out) ::              HomL(1:Nord,1:Nord),
     .                                         DHomL(1:N,1:Nord,1:Nord)
      integer ::                                       minI,maxI,a,i,ni
      double precision ::   homP(1:Nord,0:Nord-1),homR(1:Nord,0:Nord-1),
     .                                                          DS01(N)
c
c  ...clearly (minI,maxI)=(1,Nord), but the syntax is written as it is
c     because it reflects how the indexing is called from outside
      minI = 1; maxI = minI+Nord-1
c
c  ...Idec is the flag to compute x AND t derivatives
c  ...If sum of S equal 1 -> Idec=.TRUE.
      if (Idec) then
        call PolyIJacobi(S(1),1.d0,Nord,Minalpha,Idec, HomL,homP,homR)
        ni = -1
        do a=minI,maxI
          ni = ni+1
          do i=1,Nord-ni
            DHomL(1:N,a,i) = homP(a,i-1)*DS(1:N,1)
          enddo
        enddo
c
c  ...If sum of S different from 1 -> Idec=.FALSE.
      else
        call PolyIJacobi(S(1),S(0)+S(1),Nord,Minalpha,Idec,
     .                                                 HomL,homP,homR)
        ni = -1
        DS01 = DS(1:N,0)+DS(1:N,1)
        do a=minI,maxI
          ni = ni+1
          do i=1,Nord-ni
            DHomL(1:N,a,i) = homP(a,i-1)*DS(1:N,1)+homR(a,i-1)*DS01
          enddo
        enddo
      endif
c
      end subroutine HomIJacobi
c
