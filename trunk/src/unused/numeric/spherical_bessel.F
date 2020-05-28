
c--------------------------------------------------------------
c  Apr 2009
c
c  ROUTINES FOR COMPUTING SPHERICAL BESSEL FUCTIONS OF 1st AND
c  2nd TYPE DOWNLOADED FROM:

c  http://jin.ece.uiuc.edu/routines/routines.html
c--------------------------------------------------------------
c
c
c
c--------------------------------------------------------------
c  S P H .  B E S S E L  F U N C T I O N  O F  1st  T Y P E
c==============================================================
        SUBROUTINE SPHJ(N,X,NM,SJ,DJ)
C==============================================================
C       Purpose: Compute spherical Bessel functions jn(x) and
C                their derivatives
C       Input :  x --- Argument of jn(x)
C                n --- Order of jn(x)  ( n = 0,1,úúú )
C       Output:  SJ(n) --- jn(x)
C                DJ(n) --- jn'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SJ(0:N+1),DJ(0:N+1)
        NM=N
        IF (DABS(X).EQ.1.0D-100) THEN
           DO 10 K=0,N
              SJ(K)=0.0D0
10            DJ(K)=0.0D0
           SJ(0)=1.0D0
           DJ(1)=.3333333333333333D0
           RETURN
        ENDIF
        SJ(0)=DSIN(X)/X
        SJ(1)=(SJ(0)-DCOS(X))/X
        IF (N.GE.2) THEN
           SA=SJ(0)
           SB=SJ(1)
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) SJ(K)=F
              F0=F1
15            F1=F
           IF (DABS(SA).GT.DABS(SB)) CS=SA/F
           IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
           DO 20 K=0,NM
20            SJ(K)=CS*SJ(K)
        ENDIF
        DJ(0)=(DCOS(X)-DSIN(X)/X)/X
        DO 25 K=1,NM
25         DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
        RETURN
        END
c
c
        INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that the magnitude of
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END
c
c
        INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)+1
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
c
c
c--------------------------------------------------------------
c  S P H .  B E S S E L  F U N C T I O N  O F  2nd  T Y P E
c==============================================================
        SUBROUTINE SPHY(N,X,NM,SY,DY)
C==============================================================
C       Purpose: Compute spherical Bessel functions yn(x) and
C                their derivatives
C       Input :  x --- Argument of yn(x) ( x ò 0 )
C                n --- Order of yn(x) ( n = 0,1,úúú )
C       Output:  SY(n) --- yn(x)
C                DY(n) --- yn'(x)
C                NM --- Highest order computed
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SY(K)=-1.0D+300
10            DY(K)=1.0D+300
           RETURN
        ENDIF
        SY(0)=-DCOS(X)/X
        SY(1)=(SY(0)-DSIN(X))/X
        F0=SY(0)
        F1=SY(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X-F0
           SY(K)=F
           IF (DABS(F).GE.1.0D+300) GO TO 20
           F0=F1
15         F1=F
20      NM=K-1
           DY(0)=(DSIN(X)+DCOS(X)/X)/X
           DO 25 K=1,NM
25            DY(K)=SY(K-1)-(K+1.0D0)*SY(K)/X
        RETURN
        END



C$$$      subroutine spherical_bessel(Kind, N, Z, ZS, ZD)
C$$$#include "syscom.blk"
C$$$      dimension ZS(0:N+1), ZD(0:N+1)
C$$$c
C$$$      iprint = 0
C$$$      if (N.eq.0) return
C$$$c
C$$$      select case (Kind)
C$$$      case (1)
C$$$        ZS(0) = ( sin(Z)/Z )
C$$$        ZS(1) = ( (ZS(0) - cos(Z))/Z )
C$$$      case (2)
C$$$        ZS(0) = ( -cos(Z)/Z )
C$$$        ZS(1) = ( (ZS(0) - sin(Z))/Z )
C$$$      end select
C$$$c
C$$$      do i=1, N
C$$$        ZS(i+1) = (2.d0*i+1.d0)/Z*ZS(i) - ZS(i-1)
C$$$      enddo
C$$$c
C$$$      ZD(0) = -ZS(1);
C$$$      do i=1, N
C$$$        ZD(i) = ( i*ZS(i-1)-(1.d0+i)*ZS(i+1) )/(2.d0*i+1.d0)
C$$$      enddo
C$$$c
C$$$      if (iprint.eq.1) then
C$$$        write(*,*) 'spherical_bessel : value, derivative'
C$$$        do i=0, N
C$$$          write(*,*) i, ZS(i), ZD(i)
C$$$ 7000     format('i= ',i3, 2('(',f10.1,f4.1,')',2x))
C$$$        enddo
C$$$      endif
C$$$      end subroutine
