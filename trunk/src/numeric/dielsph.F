      SUBROUTINE dielsph(xyz,wk,zE)

      IMPLICIT NONE
c      double precision, parameter :: pi       = 3.141592653589793, 
c     &                               c0       = 2.99792458d8,
c     &                               my0      = 4*pi*1.d-7,
c     &                               epsilon0 = 1.d0/my0/c0**2

      double precision, parameter ::  pi       = 3.141592653589793, 
     &                                c0       = 1.d0,
     &                                my0      = 1.d0,
     &                                epsilon0 = 1.d0
      complex*16          zE(3)
      double precision::  xyz(3),wk
      double precision::  r32(2), epsilon1, epsilon2, 
     &                    epsilon3, eps(3)
      double precision::  my1, my2, my3, my(3), sigma1, sigma2, 
     &                    sigma3, sigma(3)
      character*8::       Etot

      double precision:: u0, omega
      complex*16, parameter:: im=(0.d0,1.d0)

      data r32/0.25d0,0.5d0/
      data epsilon1,epsilon2,epsilon3/4.0d0,4.d0,4.0d0/
      data my1,my2,my3               /1.0d0,1.0d0      ,1.0d0/
      data sigma1,sigma2,sigma3      /0.0d0,0.0d0      ,0.0d0/
      data Etot/'no      '/


      eps(1) = epsilon0*epsilon1
      eps(2) = epsilon0*epsilon2
      eps(3) = epsilon0*epsilon3

      my(1) = my0*my1
      my(2) = my0*my2
      my(3) = my0*my3

      sigma(1) = sigma1
      sigma(2) = sigma2
      sigma(3) = sigma3

      omega=wk*c0
      u0 = im*cdsqrt(im*my(1)*omega/(sigma(1) - im*omega*eps(1)))
ccc      write(*,*) 'u0=',u0
c      u0=1.d0
      call e_field(xyz,r32,eps,my,sigma,omega,u0,Etot,zE)

      RETURN
      END SUBROUTINE dielsph





      SUBROUTINE e_field(xyz,r32,eps,my,sigma,omega,u0,Etot,zE)

      IMPLICIT NONE
      integer:: i, j, k, l, tao
      integer, parameter:: m1 = 2, lmax = 40



      double precision:: xyz(3), r32(2), eps(3), my(3), sigma(3),
     &                   u0
      double precision:: r,x,y
      double precision:: omega, theta, phi, cphi, sphi

      complex*16, parameter:: im=(0.d0,1.d0)
      complex*16:: zE(3)
      complex*16:: kr,  setkr, kr32(4) 
      complex*16:: xx, seteta, eta(3)
      complex*16, allocatable:: mainbe(:,:),  mainha(:,:)
      complex*16, allocatable:: tbe(:), tha(:), ps(:,:,:,:)
      complex*16, allocatable:: fa(:,:)
      complex*16, allocatable:: tmpps(:,:,:,:)
      character*8:: Etot

      allocate(mainbe(4,lmax+1), mainha(4,lmax+1), 
     &         tbe(lmax+1), tha(lmax+1))
      allocate(fa(8,lmax))
      allocate(ps(3,2,m1+1,lmax+1),tmpps(3,2,m1+1,lmax+1))


      k = 0
      do j=1,2
         do i=3-j,4-j
            k = k+1
            xx = setkr(my(i), eps(i), sigma(i), omega, r32(j))
            kr32(k) = xx
            call csfhan(xx, lmax+1, tbe, tha)
            mainbe(k,:) = tbe
            mainha(k,:) = tha
         enddo                  
      enddo                    

      deallocate(tbe, tha)
      do i=1,3
         eta(i)=seteta(my(i), eps(i), sigma(i), omega)
      enddo

      do l=2,lmax-1
         call coeff(fa(1,l),fa(2,l),fa(3,l),fa(4,l),fa(5,l),fa(6,l),
     &         fa(7,l),fa(8,l),l,u0,eta,kr32,
     &         mainbe(:,l-1:l+1), mainha(:,l-1:l+1))
      enddo                    


      y = xyz(2)
      x = xyz(1)
      r = dsqrt(x**2 + y**2 + xyz(3)**2)
      ! theta = datan2(dsqrt(x**2 + y**2), xyz(3))
      theta = dacos(xyz(3)/r)
      if (theta .eq. 0) theta = 1e-8
      if ((x.eq.0d0).and.(y.eq.0d0)) then
         phi = 0d0
      else
         phi = datan2(y, x)
      end if

      kr = setkr(my(1), eps(1), sigma(1), omega, r)


      

      zE(1:3) =(0.0d0,0.0d0)
      ps(:,:,:,:)=(0.d0,0.d0)
      call PSI(m1+1, lmax+1, theta, kr, tmpps)
      cphi=dcos(phi)
      sphi=dsin(phi)
      do tao=1,2
         ps(1,tao,m1,2:lmax-1)=tmpps(1,tao,m1,2:lmax-1)*cphi
         ps(2,tao,m1,2:lmax-1)=tmpps(2,tao,m1,2:lmax-1)*cphi
         ps(3,tao,m1,2:lmax-1)=tmpps(3,tao,m1,2:lmax-1)*sphi
      enddo 
      do l=2,lmax-1
         do j=1,3 
            zE(j) = zE(j)  
     &            + fa(3,l) * ps(j,1,m1,l) 
     &            + fa(4,l) * ps(j,2,m1,l)
         enddo
      enddo
c
c     Add incoming field
c
      if (Etot(1:3) .eq. 'yes') then
         zE(1)=zE(1)
     &        +u0*dsin(theta)*dcos(phi)*cdexp(im*kr*dcos(theta))
         zE(2)=zE(2)
     &        +u0*dcos(theta)*dcos(phi)*cdexp(im*kr*dcos(theta))
         zE(3)=zE(3)
     &        -u0*dsin(phi)*cdexp(im*kr*dcos(theta))
      endif
c      write(*,*) 'x,y,z,r=',x,y,xyz(3),r
c      write(*,*) 'E=',zE

      deallocate(mainbe, mainha, ps, tmpps, fa)

      RETURN
      END SUBROUTINE e_field


C===============================================================================
C
CSubroutine coeff, by Henrik Stensson 980623
C
CCalculates the coefficients for the scattered and incoming fields.
C
CInput: order > 1, u0 = incoming field amplitude, eta, r32 = r1 and r2,
Cbe and ha = a 4 x 3 matrix with bessel and hankel function of 
Corder-2, order-1, order stored in columns. Note that order represents 
Cthe index in the bessel arrays. If l is the desired order to calculate 
Cthen: l = order - 1
C
COutput: The calculated coefficients are stored in a1, a2, f1, and f2
C
CUses: sfderive and Tcalc in this file
C
C===============================================================================
      SUBROUTINE coeff
     &(f1, f2, f3, f4, fapp1, fapp2, a1, a2, 
     & order, u0, eta, kr32, be, ha)

      implicit none
      integer:: i, order
      double precision, parameter  :: pi = 3.141592653589793
      double precision:: u0
      complex*16, parameter :: im=(0.d0,1.d0) 
      complex*16:: ne(4,3), dbe(4), dha(4), dne(4), T(4), 
     &             f1, f2, a1, a2
      complex*16:: f3, f4, fapp1, fapp2, ny, kr32(4), 
     &             eta(3), be(4,3), ha(4,3)

      do i=1,4
         ne(i,1:3) = -im*(ha(i,1:3)-be(i,1:3))
      enddo 
      call sfderive(kr32, be, ha, ne, dbe, dha, dne, order-1)
      call Tcalc(T, eta, be(:,2), ha(:,2), ne(:,2), dbe, dha, dne)
      a1 = u0*im**(order-1)*dsqrt(2*pi*(2*(order-1)+1))
      a2 = u0*im**(order-2)*dsqrt(2*pi*(2*(order-1)+1))
      f1 = T(1)*a1
      f2 = T(2)*a2
      f3 = T(3)*a1
      f4 = T(4)*a2
      ny = 0    ! eta(2)/eta(1)
      fapp1 = 0 ! ((a2-ny*a1)*ny*dbe(3)*ha(3,2)-(a1+ny*a2)*dha(3)*be(3,2))/(1+ny)/dha(3)/ha(3,2)
      fapp2 = 0 ! ((ny*a1-a2)*ny*dbe(3)*ha(3,2)-(a1+ny*a2)*dha(3)*be(3,2))/(1+ny)/dha(3)/ha(3,2)

      RETURN
      END SUBROUTINE coeff

C===============================================================================
C
CSubroutine Tcalc, by Henrik Stensson 980525
C
CCalculates T with order specified by bef, haf etc.
C
CInput: Spherical bessel functions and constant eta
C
COutput: The result T1 and T2 are stored in T
C
C===============================================================================
      SUBROUTINE Tcalc(T, eta, bef, haf, nef, dbe, dha, dne)

      implicit none
      complex*16:: T(4), eta(3), bef(4), haf(4), nef(4), 
     &             dbe(4), dha(4), dne(4)
      complex*16:: m1, m2, t1, t2, n1, n2
      complex*16:: t3, n3, t4, n4 ! , m3, m4 C remember T(4)!

      m1 = (bef(2)*eta(3)*dne(1)-nef(1)*eta(2)*dbe(2)) / 
     &     (bef(1)*eta(2)*dbe(2)-bef(2)*eta(3)*dbe(1)) 
      m2 = (nef(1)*eta(3)*dbe(2)-bef(2)*eta(2)*dne(1)) / 
     &     (bef(2)*eta(2)*dbe(1)-bef(1)*eta(3)*dbe(2))

      if ((dabs(real(m1)) .lt. 1d-16) .and. 
     &(dabs(imag(m1)+1) .lt. 1d-16)) then
         t1 = -eta(2)*dbe(3)+eta(1)*bef(3)*dbe(4)/bef(4)
         n1 = eta(2)*dha(3)-eta(1)*haf(3)*dbe(4)/bef(4)
cccc        print '(a0, \)', 's1'
      else
         t1 = (bef(4)*m1+nef(4))*eta(2)*dbe(3)-
     &        (dbe(4)*m1+dne(4))*eta(1)*bef(3)
         n1 = (dbe(4)*m1+dne(4))*eta(1)*haf(3)-
     &        (bef(4)*m1+nef(4))*eta(2)*dha(3)
      endif

      if ((dabs(real(m2)) .lt. 1d-16) .and. 
     &    (dabs(imag(m2)+1) .lt. 1d-16)) then
         t2 = -eta(1)*dbe(3)+eta(2)*bef(3)*dbe(4)/bef(4)
         n2 = eta(1)*dha(3)-eta(2)*haf(3)*dbe(4)/bef(4)
C        print '(a0, \)', 's2'
      else
         t2 = (dbe(4)*m2+dne(4))*eta(2)*bef(3)-
     &        (bef(4)*m2+nef(4))*eta(1)*dbe(3)
         n2 = (bef(4)*m2+nef(4))*eta(1)*dha(3)-
     &        (dbe(4)*m2+dne(4))*eta(2)*haf(3)
      endif

      if (n1 .eq. (0,0)) then
         T(1) = 0
ccc         print '(a, \)', 'z'
      else
         T(1) = t1/n1
      endif

      if (n2 .eq. (0,0)) then
         T(2) = 0
ccc         print '(a, \)', 'z'
      else
         T(2) = t2/n2
      endif

C
CIndices 3 and 4 only needed when comparing with the results from G. Kristensson
C
      t3 = eta(2)*bef(4)*dbe(3)-eta(1)*dbe(4)*bef(3)
      n3 = eta(1)*haf(3)*dbe(4)-eta(2)*dha(3)*bef(4)
      t4 = eta(1)*bef(4)*dbe(3)-eta(2)*dbe(4)*bef(3)
      n4 = eta(2)*haf(3)*dbe(4)-eta(1)*dha(3)*bef(4)
      T(3) = t3/n3
      T(4) = t4/n4

C
CAlternative calculation of t1,n2, etc.
C
Cm3 = (eta(2)*dbe(1)*bef(2)-eta(3)*dbe(2)*bef(1))/(eta(3)*dbe(2)*nef(1)-eta(2)*dne(1)*bef(2))
Cm4 = (eta(3)*dbe(1)*bef(2)-eta(2)*dbe(2)*bef(1))/(eta(2)*dbe(2)*nef(1)-eta(3)*dne(1)*bef(2))
Ct3 = eta(1)*dbe(3)*(bef(4)+m3*nef(4))-eta(2)*bef(3)*(dbe(4)+m3*dne(4))
Cn3 = eta(2)*haf(3)*(dbe(4)+m3*dne(4))-eta(1)*dha(3)*(bef(4)+m3*nef(4))
Ct4 = eta(2)*dbe(3)*(bef(4)+m4*nef(4))-eta(1)*bef(3)*(dbe(4)+m4*dne(4))
Cn4 = eta(1)*haf(3)*(dbe(4)+m4*dne(4))-eta(2)*dha(3)*(bef(4)+m4*nef(4))
CT(3) = t3/n3
CT(4) = t4/n4

      RETURN
      END SUBROUTINE Tcalc
C===============================================================================
C
CSubroutine sfderive, by Henrik Stensson 980610
C
CCalculates the derivative of spherical bessel functions for use only with Tcalc
C
CInput: Spherical bessel functions, radius and order of the functions. 
CThe input functions needs to be 4 x 3 matrices with the order l-1 to l+1 
Cstored in the columns.
C
COutput: The results are put in dbe, dha and dne.
C
C===============================================================================
      SUBROUTINE sfderive(kr32, be, ha, ne, dbe, dha, dne, ord)

      IMPLICIT NONE
      complex*16:: be(4, 3), ha(4, 3), ne(4, 3), dbe(4), 
     &             dha(4), dne(4), bep(4)
      complex*16:: hap(4), nep(4), kr32(4)
      integer:: k, ord

      do k=1,4
         bep(k) = (ord*be(k,1)-(ord+1)*be(k,3))/(2*ord+1)
         hap(k) = (ord*ha(k,1)-(ord+1)*ha(k,3))/(2*ord+1)
         nep(k) = (ord*ne(k,1)-(ord+1)*ne(k,3))/(2*ord+1)
      enddo                

      do k=1,4
         dbe(k) = be(k,2)/kr32(k)+bep(k)
         dha(k) = ha(k,2)/kr32(k)+hap(k)
         dne(k) = ne(k,2)/kr32(k)+nep(k)
      enddo 

      RETURN
      END SUBROUTINE sfderive

C===============================================================================
C
C
CFunction setkr, by Henrik Stensson 980420
C
CCalculates the variable kr for use with the spherical bessel functions.
C
CInput: m = my, e = epsilon, s = sigma, w = omega, r = radius.
C
COutput: k * r, where k^2 = w^2 * my * epsilon * (1 + i * sigma / w / epsilon) 
C
C===============================================================================
      COMPLEX*16 FUNCTION setkr(m, e, s, w, r)

      complex*16, parameter :: im=(0.d0,1.d0)

      double precision:: m, e, s, w, r

      setkr = r*w*cdsqrt(m*(e+im*s/w))
      
      RETURN
      END FUNCTION setkr

C===============================================================================
C
CFunction seteta, by Henrik Stensson 980525
C
CCalculates the constant eta for use with emscattering problems.
C
CInput: m = my, e = epsilon, s = sigma, w = omega.
C
COutput: sqrt(my / epsilon~) / i, where epsilon~ = epsilon + i * sigma / w
C
CNote: The only subroutine that uses this function is 'Tcalc'. Therefore it aint
Cnecessaire to divide by im.
C
C===============================================================================
      COMPLEX*16 FUNCTION seteta(m, e, s, w)

      complex*16, parameter :: im=(0.d0,1.d0)

      double precision:: m, e, s, w

      seteta = cdsqrt(m/(e+im*s/w)) !/im

      RETURN
      END FUNCTION seteta

      SUBROUTINE CSFHAN(X,LDIM,B,H)
      IMPLICIT   NONE
      INTEGER    L1,LDIM,L
      COMPLEX*16 B(1),B1,B2,Y1,Y2,H(*),CI,X,COS,SIN
      COMPLEX*16, ALLOCATABLE :: Y(:)

      ALLOCATE (Y(LDIM))
      
      SIN=CDSIN(X)
      COS=CDCOS(X)
      Y1=-COS/X
      Y2=(Y1-SIN)/X
      Y(1)=Y1
      Y(2)=Y2
      IF(cdabs(X).GT.LDIM) GOTO 15
      CALL CSFBES(X,LDIM-1,B1)
      CALL CSFBES(X,LDIM-2,B2)
      B(LDIM)=B1
      B(LDIM-1)=B2
      DO 11 L1=3,LDIM
        L=LDIM-L1+1
        B(L)=(2*L+1)*B2/X-B1
        B1=B2
        B2=B(L)
        Y(L1)=(2*L1-3)*Y2/X-Y1
        Y1=Y2
        Y2=Y(L1)
11    CONTINUE
      GOTO 19
15    B1=SIN/X
      B2=(B1-COS)/X
      B(1)=B1
      B(2)=B2
      DO 18 L1=3,LDIM
        B(L1)=(2*L1-3)*B2/X-B1
        B1=B2
        B2=B(L1)
        Y(L1)=(2*L1-3)*Y2/X-Y1
        Y1=Y2
        Y2=Y(L1)
18    CONTINUE
19    CI=(0.D0,1.D0)
      DO 22 L1=1,LDIM
        H(L1)=B(L1)+CI*Y(L1)
22    CONTINUE
      DEALLOCATE(Y)
      RETURN
      END


C===============================================================================
C      
C Subroutine PSI
C
C COMPUTES PSI AT KS AND TH 
C
C===============================================================================
      SUBROUTINE PSI(MMAX1,LMAX1,TH,KS,PS)
      IMPLICIT   NONE
      INTEGER    MMAX1,LMAX1,LMAX2,L1,IVEC,J1,M,M1,M10
      DOUBLE PRECISION TH,SI,CO,R,A1,A2,A3


c      COMPLEX*16 PS(3,2,MMAX1,LMAX1)
      COMPLEX*16 PS(3,2,MMAX1,LMAX1)
      COMPLEX*16 KS

      COMPLEX*16, ALLOCATABLE :: BS(:),HS(:)      
      COMPLEX*16, ALLOCATABLE  :: H1(:),H2(:), H3(:)
      DOUBLE PRECISION, ALLOCATABLE  :: G(:,:,:), P(:)

c      ALLOCATE( H1(LMAX1), H2(LMAX1), H3(LMAX1) )
c      ALLOCATE( BS(LMAX1+1),HS(LMAX1+1) )
c      ALLOCATE( G(3,MMAX1,LMAX1) , P(LMAX1+1) )

      ALLOCATE( H1(LMAX1), H2(LMAX1), H3(LMAX1) )
      ALLOCATE( BS(LMAX1+1),HS(LMAX1+1) )
      ALLOCATE( G(3,MMAX1,LMAX1+1) , P(LMAX1+1) )
       
      LMAX2=LMAX1+1
      SI=DSIN(TH)
      CO=DCOS(TH)
      CALL CSFHAN(KS,LMAX2,BS,HS)
      call gamma(MMAX1,LMAX2,G)
      DO 10 L1=1,LMAX1
      H1(L1)=L1*HS(L1)/KS-HS(L1+1)
      H2(L1)=H1(L1)-2*HS(L1)/KS
      H3(L1)=(2*L1*(L1-2)/KS-KS)*HS(L1)+2*HS(L1+1)
   10 CONTINUE
      DO 11 IVEC=1,3
      DO 11 J1=1,2
      PS(IVEC,J1,1,1)=(0.D0,0.D0)
   11 CONTINUE
      DO 26 M1=1,MMAX1
         M=M1-1
         M10=M1
         IF(M.EQ.0) M10=2
         CALL LEG(TH,M,LMAX2,P)
         DO 26 L1=M10,LMAX1
            R=DSQRT(L1*(L1-1.D0))
            A1=M*G(3,M1,L1)*P(L1)/SI
            A2=-G(3,M1,L1)*(L1*CO*P(L1)-(L1-M)*P(L1+1))/SI
            A3=G(3,M1,L1)*P(L1)
            PS(1,1,M1,L1)=(0.D0,0.D0)
            PS(2,1,M1,L1)=HS(L1)*A1/R
            PS(3,1,M1,L1)=-HS(L1)*A2/R
            PS(1,2,M1,L1)=R*HS(L1)/KS*A3
            PS(2,2,M1,L1)=H1(L1)*A2/R
            PS(3,2,M1,L1)=-H1(L1)*A1/R
   26 CONTINUE
      DEALLOCATE(H1,H2,H3,HS,BS,G,P)
      RETURN
      END

C===============================================================================
C
C Subroutine GAMMA
C
C COMPUTES NORMALIZATION CONSTANTS FOR REAL SPHERICAL
C HARMONICS FOR L<LMAX1 AND M<MMAX1
C
C===============================================================================
      SUBROUTINE GAMMA(MDIM,LDIM,G)
      DOUBLE PRECISION G(3,MDIM,LDIM),X,PI,DSQRT

      DOUBLE PRECISION, ALLOCATABLE :: FCT(:)
      ALLOCATE (FCT(MDIM+LDIM))
 
      PI=3.141592653589793D0
      FCT(1)=1.D0
      LM=LDIM+MDIM-2
      DO 11 I=1,LM
   11 FCT(I+1)=I*FCT(I)
      DO 22 M1=1,MDIM
         MEPS=2
         IF(M1.EQ.1) MEPS=1
         DO 22 L1=M1,LDIM
            LMM=L1-M1+1
            LMP=L1+M1-1
            G(3,M1,L1)=DSQRT(MEPS*(2*L1-1)*FCT(LMM)/(4*PI*FCT(LMP)))
            X=DSQRT(L1*(L1-1.D0))
            IF(L1.NE.1) X=1/X
            G(1,M1,L1)=X*G(3,M1,L1)
 22         G(2,M1,L1)=G(1,M1,L1)

      DEALLOCATE(FCT)
      RETURN
      END
C===============================================================================
C
C Subroutine LEG
C
C COMPUTES ASSOCIATED LEGENDRE FUNCTIONS PLM(COS(TH))
C WITH M FIXED AND L BETWEEN M AND LDIM-1
C
C===============================================================================
      SUBROUTINE LEG(TH,M,LDIM,PLM)
      DOUBLE PRECISION TH,PLM(LDIM),P0,P1,SI,CO,DSIN,DCOS
      P0=0.D0
      P1=1.D0
      IF(M.EQ.0) GOTO 13
      SI=DSIN(TH)
      DO 11 I=1,M
   11 P1=P1*SI*(2*M+1-I)/2
   13 PLM(M+1)=P1
      IF(M+1.EQ.LDIM) RETURN
      CO=DCOS(TH)
      M2=M+2
      DO 17 L1=M2,LDIM
      PLM(L1)=(CO*P1*(2*L1-3)-P0*(L1+M-2))/(L1-M-1)
      P0=P1
   17 P1=PLM(L1)
      RETURN
      END
C===============================================================================
C      
C Subroutine CSFBES
C
C COMPUTES THE SPHERICAL BESSEL FUNCTION JL(X), L>0
C
C===============================================================================
      SUBROUTINE CSFBES(X,L,B)
      IMPLICIT   NONE
      INTEGER    L,J,K,I
      COMPLEX*16 X,B,A,X2
      
      B=1.D0
      A=1.D0
      X2=-X*X/2
      J=1
      K=2*L+3
      DO 15 I=1,100
      A=A*X2/J/K
      B=B+A
      IF(CDABS(A).LT.1.D-16*CDABS(B)) GOTO 25
      J=J+1
   15 K=K+2
      STOP
   25 J=3
      DO 29 I=1,L
      B=B*X/J
   29 J=J+2
      RETURN
      END
