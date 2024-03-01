! Routines:
!  - PolyLegendre
!  - PolyILegendre
!  - PolyJacobi
!  - PolyIJacobi
!  - HomLegendre
!  - HomILegendre
!  - HomJacobi
!  - HomIJacobi
!
!----------------------------------------------------------------------
!
!     routine name      - PolyLegendre
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!> @brief         - routine returns values of shifted scaled
!                         Legendre polynomials
!
!     arguments:
!
!     in:
!             X         - coordinate from [0,1]
!             T         - scaling parameter
!             Nord      - polynomial order
!
!     out:
!             P         - polynomial values
!
!----------------------------------------------------------------------
!
   subroutine PolyLegendre(X,T,Nord, P)
!
      implicit none
      integer,          intent(in)  ::                             Nord
      double precision, intent(in)  ::                              X,T
      double precision, intent(out) ::                        P(0:Nord)
      integer ::                                                      i
      double precision ::                                          tt,y
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!  ...i stands for the order of the polynomial, stored in P(i)
!  ...lowest order case (order 0)
      P(0) = 1.d0
!  ...first order case (order 1) if necessary
      if (Nord.ge.1) then
        y = 2.d0*X - T
        P(1) = y
      endif
!  ...higher order if necessary - use recurrence formula
      if (Nord.ge.2) then
        tt = T**2
        do i=2,Nord
          P(i) = (2*i-1)*y*P(i-1) - (i-1)*tt*P(i-2)
          P(i) = P(i)/i
        enddo
      endif
!
#if HP3D_DEBUG
!  ...catching problems (debugging)
      if (iprint.eq.1) then
        write(*,7002) Nord, X,T
 7002   format('PolyLegendre: Nord = ',i2,' X,T = ',2F8.3)
        do i=0,Nord
          write(*,7003) i,P(i)
 7003     format('i = ',i2,' P = ',e25.15)
        enddo
        call pause
      endif
#endif
!
   end subroutine PolyLegendre
!
!----------------------------------------------------------------------
!
!     routine name      - PolyILegendre
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!> @brief         - routine returns values of shifted scaled
!                         integrated Legendre polynomials and their
!                         derivatives starting with p=2
!
!     arguments:
!
!     in:
!             X         - coordinate from [0,1]
!             T         - scaling parameter
!             Nord      - polynomial order
!             Idec      - decision flag to compute:
!                       = FALSE polynomials with x and t derivatives
!                       = TRUE  polynomials with x derivatives only
!
!     out:
!             L         - polynomial values
!             P         - derivatives in x
!             R         - derivatives in t
!
!----------------------------------------------------------------------
!
   subroutine PolyILegendre(X,T,Nord,Idec, L,P,R)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                             Nord
      double precision, intent(in)  ::                              X,T
      double precision, intent(out) ::    L(2:Nord),P(Nord-1),R(Nord-1)
      integer ::                                                i,ifact
      double precision ::                              ptemp(0:Nord),tt
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!  ...calling Legendre for required information
      call PolyLegendre(X,T,Nord, ptemp)
      P = ptemp(1:Nord-1)
!
!  ...Integrated polynomial of order i is stored in L(i)
      tt = T**2
!
!  ...simplified case: no need to compute R
      if (Idec) then
        do i=2,Nord
          ifact = 4*i-2
          L(i) = (ptemp(i) - tt*ptemp(i-2))/ifact
        enddo
!
!  ...general case: compute R
      else
        do i=2,Nord
          ifact = 4*i-2
          L(i) = (ptemp(i) - tt*ptemp(i-2))/ifact
          R(i-1) = -(ptemp(i-1)+T*ptemp(i-2))/2
        enddo
      endif
!
#if HP3D_DEBUG
!  ...catching problems (debugging)
      if (iprint.eq.1) then
        write(*,7002) Idec,Nord, X,T
 7002   format('PolyILegendre: Idec = ',i1,' Nord = ',i2, &
               ' X,T = ',2F8.3)
        do i=2,Nord
          select case(Idec)
          case(.true.)
            write(*,7003) i,L(i),P(i)
 7003       format('i = ',i2,' L,P = ',2e25.15)
          case default
            write(*,7004) i,L(i),P(i),R(i)
 7004       format('i = ',i2,' L,P,R, = ',3e25.15)
          end select
        enddo
        call pause
      endif
#endif
!
   end subroutine PolyILegendre
!
!
!----------------------------------------------------------------------
!
!     routine name      - PolyJacobi
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!> @brief         - routine returns values of shifted scaled
!                         Jacobi polynomials P^\alpha_i. Result is a
!                         'half' of a  matrix with each row
!                         associated to a fixed alpha. Alpha grows
!                         by 2 in each row.
!
!     arguments:
!
!     in:
!             X         - coordinate from [0,1]
!             T         - scaling parameter
!             Nord      - max polynomial order
!             Minalpha  - first row value of alpha (integer)
!
!     out:
!             P         - polynomial values
!
!----------------------------------------------------------------------
!
   subroutine PolyJacobi(X,T,Nord,Minalpha, P)
!
      implicit none
      integer,          intent(in)  ::                    Nord,Minalpha
      double precision, intent(in)  ::                              X,T
      double precision, intent(out) ::                 P(0:Nord,0:Nord)
      integer ::       minI,maxI,i,ni,a,aa,al,ai,bi,ci,di,alpha(0:Nord)
      double precision ::                                          y,tt
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!  ...clearly (minI,maxI)=(0,Nord), but the syntax is written as it is
!     because it reflects how the indexing is called from outside
      minI = 0; maxI = minI+Nord
!
!  ...in our work Minalpha>=1
      if (Minalpha.lt.1) then
        write(*,7001) Minalpha
 7001   format('PolyJacobi: Minalpha = ',i3)
        stop 1
      endif
!
!  ...create vector alpha first
      do a=minI,maxI
         alpha(a) = Minalpha+2*(a-minI)
      enddo
!
!  ...initiate first column (order 0)
      P(minI:maxI,0) = 1.d0
!  ...initiate second column (order 1) if necessary
      if (Nord.ge.1) then
        y = 2.d0*X - T
        P(minI:maxI-1,1) = y+alpha(minI:maxI-1)*X
      endif
!  ...fill the last columns if necessary
      if (Nord.ge.2) then
        tt = T**2
        ni = -1
        do a=minI,maxI-2
          al=alpha(a)
          aa = al**2
          ni=ni+1
!      ...use recursion in order, i, to compute P^alpha_i for i>=2
          do i=2,Nord-ni
            ai = 2*i*(i+al)*(2*i+al-2)
            bi = 2*i+al-1
            ci = (2*i+al)*(2*i+al-2)
            di = 2*(i+al-1)*(i-1)*(2*i+al)
!
            P(a,i) = bi*(ci*y+aa*T)*P(a,i-1)-di*tt*P(a,i-2)
            P(a,i) = P(a,i)/ai
          enddo
        enddo
      endif
!
#if HP3D_DEBUG
!  ...catching problems (debugging)
      if (iprint.eq.1) then
        write(*,7003) Nord, X,T
 7003   format('PolyJacobi: Nord = ',i2,' X,T = ',2F8.3)
        do a=minI,maxI
          write(*,7004) a,P(a,0:Nord)
 7004     format(' P(',i2,',0:Nord) = ',10e12.5)
        enddo
        call pause
      endif
#endif
!
   end subroutine PolyJacobi
!
!----------------------------------------------------------------------
!
!     routine name      - PolyIJacobi
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!> @brief         - routine returns values of integrated
!                         shifted scaled Jacobi polynomials and
!                         their derivatives starting with p=1
!                         Result is 'half' of a  matrix
!                         with each row  associated to a fixed alpha.
!                         Alpha grows by 2 in each row.
!
!     arguments:
!
!     in:
!             X         - coordinate from [0,1]
!             T         - scaling parameter
!             Nord      - max polynomial order
!             Minalpha  - first row value of alpha (integer)
!             Idec      - decision flag to compute:
!                       = FALSE polynomials with x and t derivatives
!                       = TRUE  polynomials with x derivatives only
!
!     out:
!             L        - polynomial values
!             P        - derivatives in x (Jacobi polynomials)
!             R        - derivatives in t
!
!----------------------------------------------------------------------
!
   subroutine PolyIJacobi(X,T,Nord,Minalpha,Idec, L,P,R)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                    Nord,Minalpha
      double precision, intent(in)  ::                              X,T
      double precision, intent(out) ::                 L(1:Nord,1:Nord), &
                                  P(1:Nord,0:Nord-1),R(1:Nord,0:Nord-1)
      integer ::                    minI,maxI,i,ni,a,al,tia,tiam1,tiam2, &
                                                          alpha(1:Nord)
      double precision ::            ai,bi,ci,tt,ptemp(1:Nord+1,0:Nord)
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!  ...clearly (minI,maxI)=(1,Nord), but the syntax is written as it is
!     because it reflects how the indexing is called from outside
      minI = 1; maxI = minI+Nord-1
!
!  ...in our work Minalpha>=1
      if (Minalpha.lt.1) then
        write(*,7001) Minalpha
 7001   format('PolyIJacobi: Minalpha = ',i3)
        stop 1
      endif
!
!  ...calling Jacobi for required information
      call PolyJacobi(X,T,Nord,Minalpha, ptemp)
!  ...define P. Note that even though P is defined at all entries,
!     because of the way Jacobi computes ptemp, only the necessary entries,
!     and those on the first subdiagonal (which are never used later)
!     are actually accurate.
      P = ptemp(minI:maxI,0:Nord-1)
!
!  ...create vector alpha first
      do a=minI,maxI
         alpha(a) = Minalpha+2*(a-minI)
      enddo
!
!  ...initiate first column (order 1 in L)
      L(minI:maxI,1) = X
!
!  ...simplified case, do not compute R
      if (Idec) then
!  ...fill the last columns if necessary
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
!
              L(a,i) = ai*ptemp(a,i)+bi*T*ptemp(a,i-1) &
                    -ci*tt*ptemp(a,i-2)
!              P(a,i-1) =  ptemp(a,i-1)
            enddo
          enddo
        endif
!
!  ...general case; compute R
      else
      R(minI:maxI,0) = 0.d0
!  ...fill the last columns if necessary
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
!
              L(a,i) = ai*ptemp(a,i)+bi*T*ptemp(a,i-1) &
                    -ci*tt*ptemp(a,i-2)
!              P(a,i-1) =  ptemp(a,i-1)
              R(a,i-1) = -(i-1)*(ptemp(a,i-1)+T*ptemp(a,i-2))
              R(a,i-1) = R(a,i-1)/tiam2
            enddo
          enddo
        endif
      endif
!
#if HP3D_DEBUG
!  ...catching problems (debugging)
      if (iprint.eq.1) then
        write(*,7003) minI,Nord, X,T
 7003   format('PolyIJacobi: minI = ',i2, &
               ' Nord = ',i2,' X,T = ',2F8.3)
        do a=minI,maxI
          al = alpha(a)
          write(*,7004) a,al,L(a,1:Nord)
 7004     format('a = ',i1,' alpha = ',i2, &
               ' L(a,1:Nord)   = ',10e12.5)
        enddo
        write(*,*) '  '
        do a=minI,maxI
          al = alpha(a)
          write(*,7005) a,al,P(a,0:Nord-1)
 7005     format('a = ',i1,' alpha = ',i2, &
               ' P(a,0:Nord-1) = ',10e12.5)
        enddo
        if (.not.Idec) then
          write(*,*) '  '
          do a=minI,maxI
            al = alpha(a)
            write(*,7006) a,al,R(a,0:Nord-1)
 7006       format('a = ',i1,' alpha = ',i2, &
                   ' R(a,0:Nord-1) = ',10e12.5)
          enddo
        endif
        call pause
      endif
#endif
!
   end subroutine PolyIJacobi
!
!----------------------------------------------------------------------
!
!     routine name      - HomLegendre
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!> @brief         - routine returns values of homogenized
!                         Legendre polynomials
!
!     arguments:
!
!     in:
!             S         - affine(like) coordinates
!             Nord      - polynomial order
!
!     out:
!             HomP      - polynomial values
!
!----------------------------------------------------------------------
!
   subroutine HomLegendre(S,Nord, HomP)
!
      implicit none
      integer,          intent(in)  ::                             Nord
      double precision, intent(in)  ::                           S(0:1)
      double precision, intent(out) ::                     HomP(0:Nord)
!
!  ...simply the definition of homogenized polynomials
      call PolyLegendre(S(1),S(0)+S(1),Nord, HomP)
!
   end subroutine HomLegendre
!
!
!----------------------------------------------------------------------
!
!     routine name      - HomILegendre
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!> @brief         - routine returns values of homogenized
!                         integrated Legendre polynomials and their
!                         gradient (wrt to affine like coordinates)
!
!     arguments:
!
!     in:
!             S         - (s0,s1) affine(like) coordinates
!             DS        - gradients of S (in R^N)
!             Nord      - polynomial order
!             Idec      - decision flag to compute:
!                         = FALSE s0+s1 != 1 -> general case
!                         = TRUE  s0+s1  = 1 -> simple case
!             N         - number of spatial dimensions (R^N)
!
!     out:
!             HomL        - polynomial values
!             DHomL       - gradients of L in R^N
!
!----------------------------------------------------------------------
!
   subroutine HomILegendre(S,DS,Nord,Idec,N, HomL,DHomL)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:1),DS(1:N,0:1)
      double precision, intent(out) ::   HomL(2:Nord),DHomL(1:N,2:Nord)
      integer ::                                                      i
      double precision ::       homP(1:Nord-1),homR(1:Nord-1),DS01(1:N)
!
!  ...Idec is the flag to compute x AND t derivatives
!  ...If sum of S equal 1 -> Idec=.true.
      if (Idec) then
        call PolyILegendre(S(1),1.d0,Nord,Idec, HomL,homP,homR)
        do i=2,Nord
          DHomL(1:N,i) = homP(i-1)*DS(1:N,1)
        enddo
!
!  ...If sum of S different from 1 -> Idec=.false.
      else
        call PolyILegendre(S(1),S(0)+S(1),Nord,Idec, HomL,homP,homR)
        DS01 = DS(1:N,0)+DS(1:N,1)
        do i=2,Nord
          DHomL(1:N,i) = homP(i-1)*DS(1:N,1)+homR(i-1)*DS01
        enddo
      endif
!
   end subroutine HomILegendre
!
!
!----------------------------------------------------------------------
!
!     routine name      - HomJacobi
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!> @brief         - routine returns values of homogenized
!                         Jacobi polynomials P^\alpha_i. Result is a
!                         'half' of a  matrix with each row
!                         associated to a fixed alpha. Alpha grows
!                         by 2 in each row.
!
!     arguments:
!
!     in:
!             S         - affine(like) coordinates
!             Nord      - max polynomial order
!             Minalpha  - first row value of alpha (integer)
!
!     out:
!             HomP      - polynomial values
!
!----------------------------------------------------------------------
!
   subroutine HomJacobi(S,Nord,Minalpha, HomP)
!
      implicit none
      integer,          intent(in)  ::                    Nord,Minalpha
      double precision, intent(in)  ::                           S(0:1)
      double precision, intent(out) ::              HomP(0:Nord,0:Nord)
      integer ::                                              minI,maxI
!
!  ...clearly (minI,maxI)=(0,Nord), but the syntax is written as it is
!     because it reflects how the indexing is called from outside
      minI = 0; MaxI = MinI+Nord
!
!  ...simply the definition of homogenized polynomials
      call PolyJacobi(S(1),S(0)+S(1),Nord,Minalpha, HomP)
!
   end subroutine HomJacobi
!
!
!----------------------------------------------------------------------
!
!     routine name      - HomIJacobi
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 14
!
!> @brief         - routine returns values of integrated
!                         homogenized Jacobi polynomials and
!                         their gradients.
!                         Result is 'half' of a  matrix
!                         with each row  associated to a fixed alpha.
!                         Alpha grows by 2 in each row.
!
!     arguments:
!
!     in:
!             S         - (s0,s1) affine(like) coordinates
!             DS        - gradients of S (in R^N)
!             Nord      - max polynomial order
!             Minalpha  - first row value of alpha (integer)
!             Idec      - decision flag to compute:
!                         = FALSE s0+s1 != 1 -> general case
!                         = TRUE  s0+s1  = 1 -> simple case
!             N         - number of spatial dimensions (R^N)
!
!     out:
!             HomL      - polynomial values
!             DHomL     - derivatives in x (Jacobi polynomials)
!
!----------------------------------------------------------------------
!
   subroutine HomIJacobi(S,DS,Nord,Minalpha,Idec,N, HomL,DHomL)
!
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                  Nord,Minalpha,N
      double precision, intent(in)  ::               S(0:1),DS(1:N,0:1)
      double precision, intent(out) ::              HomL(1:Nord,1:Nord), &
                                               DHomL(1:N,1:Nord,1:Nord)
      integer ::                                       minI,maxI,a,i,ni
      double precision ::   homP(1:Nord,0:Nord-1),homR(1:Nord,0:Nord-1), &
                                                                DS01(N)
!
!  ...clearly (minI,maxI)=(1,Nord), but the syntax is written as it is
!     because it reflects how the indexing is called from outside
      minI = 1; maxI = minI+Nord-1
!
!  ...Idec is the flag to compute x AND t derivatives
!  ...If sum of S equal 1 -> Idec=.true.
      if (Idec) then
        call PolyIJacobi(S(1),1.d0,Nord,Minalpha,Idec, HomL,homP,homR)
        ni = -1
        do a=minI,maxI
          ni = ni+1
          do i=1,Nord-ni
            DHomL(1:N,a,i) = homP(a,i-1)*DS(1:N,1)
          enddo
        enddo
!
!  ...If sum of S different from 1 -> Idec=.false.
      else
        call PolyIJacobi(S(1),S(0)+S(1),Nord,Minalpha,Idec, &
                                                       HomL,homP,homR)
        ni = -1
        DS01 = DS(1:N,0)+DS(1:N,1)
        do a=minI,maxI
          ni = ni+1
          do i=1,Nord-ni
            DHomL(1:N,a,i) = homP(a,i-1)*DS(1:N,1)+homR(a,i-1)*DS01
          enddo
        enddo
      endif
!
   end subroutine HomIJacobi
!
