        SUBROUTINE CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY)
C
C       ==========================================================
C       Purpose: Compute spherical Bessel functions jn(z) & yn(z)
C                and their derivatives for a complex argument
C       Input :  z --- Complex argument
C                n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
C       Output:  CSJ(n) --- jn(z)
C                CDJ(n) --- jn'(z)
C                CSY(n) --- yn(z)
C                CDY(n) --- yn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ==========================================================
C
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION A0
        DIMENSION CSJ(0:N),CDJ(0:N),CSY(0:N),CDY(0:N)
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-60) THEN
           DO 10 K=0,N
              CSJ(K)=0.0D0
              CDJ(K)=0.0D0
              CSY(K)=-1.0D+300
10            CDY(K)=1.0D+300
           CSJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(.333333333333333D0,0.0D0)
           RETURN
        ENDIF
        CSJ(0)=CDSIN(Z)/Z
        CSJ(1)=(CSJ(0)-CDCOS(Z))/Z
        IF (N.GE.2) THEN
           CSA=CSJ(0)
           CSB=CSJ(1)
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF0=0.0D0
           CF1=1.0D0-100
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF1/Z-CF0
              IF (K.LE.NM) CSJ(K)=CF
              CF0=CF1
15            CF1=CF
           IF (CDABS(CSA).GT.CDABS(CSB)) CS=CSA/CF
           IF (CDABS(CSA).LE.CDABS(CSB)) CS=CSB/CF0
           DO 20 K=0,NM
20            CSJ(K)=CS*CSJ(K)
        ENDIF
        CDJ(0)=(CDCOS(Z)-CDSIN(Z)/Z)/Z
        DO 25 K=1,NM
25         CDJ(K)=CSJ(K-1)-(K+1.0D0)*CSJ(K)/Z
        CSY(0)=-CDCOS(Z)/Z
        CSY(1)=(CSY(0)-CDSIN(Z))/Z
        CDY(0)=(CDSIN(Z)+CDCOS(Z)/Z)/Z
        CDY(1)=(2.0D0*CDY(0)-CDCOS(Z))/Z
        DO 30 K=2,NM
           IF (CDABS(CSJ(K-1)).GT.CDABS(CSJ(K-2))) THEN
              CSY(K)=(CSJ(K)*CSY(K-1)-1.0D0/(Z*Z))/CSJ(K-1)
           ELSE
              CSY(K)=(CSJ(K)*CSY(K-2)-(2.0D0*K-1.0D0)/Z**3)/CSJ(K-2)
           ENDIF
30      CONTINUE
        DO 35 K=2,NM
35         CDY(K)=CSY(K-1)-(K+1.0D0)*CSY(K)/Z
        RETURN
        END


