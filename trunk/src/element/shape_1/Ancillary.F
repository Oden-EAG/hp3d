c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c                                  EDGES
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c     routine name      - AncPhiE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - compute edge H1 ancillary functions and
c                         their gradients
c
c     arguments:
c
c     in:
c             S         - (s0,s1) affine coordinates associated to edge
c             DS        - gradients of S in R^N
c             Nord      - polynomial order
c             Idec      - Binary flag:
c                         = FALSE  s0+s1 != 1
c                         = TRUE   s0+s1  = 1
c             N         - spatial dimension
c
c     out:
c             PhiE      - values of edge H1 ancillary functions
c             DPhiE     - gradients of edge H1 ancillary functions
c
c----------------------------------------------------------------------
c
      subroutine AncPhiE(S,DS,Nord,Idec,N, PhiE,DPhiE)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:1),DS(1:N,0:1)
      double precision, intent(out) ::   PhiE(2:Nord),DPhiE(1:N,2:Nord)
      integer ::                                              minI,maxI
c
c  ...local parameters
      minI = 2
      maxI = Nord
c
      if (N.lt.1) then
        write(*,7001) N
 7001   format('AncPhiE: N = ',i2)
      endif
c
c  ...these are precisely the homogenized Legendre polynomials
      call HomILegendre(S,DS,Nord,Idec,N, PhiE,DPhiE)
c
      end subroutine AncPhiE
c
c----------------------------------------------------------------------
c
c     routine name      - AncEE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - compute edge Hcurl ancillary functions and
c                         their curls
c
c     arguments:
c
c     in:
c             S         - (s0,s1) affine coordinates associated to edge
c             DS        - derivatives of S in R^N
c             Nord      - polynomial order
c             Idec      - Binary flag:
c                         = FALSE  s0+s1 != 1
c                         = TRUE   s0+s1  = 1
c             N         - spatial dimension
c
c     out:
c             EE        - edge Hcurl ancillary functions
c             CurlEE    - curls of edge Hcurl ancillary functions
c
c----------------------------------------------------------------------
c
      subroutine AncEE(S,DS,Nord,Idec,N, EE,CurlEE)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:1),DS(1:N,0:1)
      double precision, intent(out) ::                 EE(1:N,0:Nord-1),
     .                                         CurlEE(1:2*N-3,0:Nord-1)
      integer ::                                      minI,maxI,Ncurl,i
      double precision ::    homP(0:Nord-1),whiE(1:N),curlwhiE(1:2*N-3)
c
c  ...local parameters
      minI = 0
      maxI = Nord-1
      Ncurl = 2*N-3
c
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncEE: N = ',i2)
      endif
c
c  ...extract homogenized Legendre polyomials first
      call HomLegendre(S,maxI, homP)
c
c  ...simplified case
      if (Idec) then
        do i=minI,maxI
          EE(1:N,i) = homP(i)*DS(1:N,1)
        enddo
c    ...no need to compute Whitney function or curl
        CurlEE(1:Ncurl,minI:maxI) = 0.d0
c
c  ...in general
      else
c    ...lowest order Whitney function and its curl
        whiE = S(0)*DS(1:N,1)-S(1)*DS(1:N,0)
        call cross(N,DS(1:N,0),DS(1:N,1), curlwhiE)
c    ...now construct the higher order elements
        do i=minI,maxI
          EE(1:N,i) = homP(i)*whiE
          CurlEE(1:Ncurl,i) = (i+2)*homP(i)*curlwhiE
        enddo
      endif
c
      end subroutine AncEE
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c                         QUADRILATERAL FACES
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c     routine name      - AncPhiQuad
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - compute quadrilateral face H1 ancillary
c                         functions and their gradients
c
c     arguments:
c
c     in:
c             ST       - affine coordinates associated to face
c                        2x2 matrix [s0,s1;t0,t1]
c             DST      - gradients of ST
c             Nord     - (NordS,NordT) vector polynomial order
c             Idec     - (IdecS,IdecT) vector binary flag:
c                        IdecS = FALSE if s0+s1 != 1
c                        IdecS = TRUE  if s0+s1  = 1, same with IdecT
c             N        - spatial dimension
c
c     out:
c             PhiQuad  - quad H1 ancillary functions
c             DPhiQuad - gradients of quad H1 ancillary functions
c
c----------------------------------------------------------------------
c
      subroutine AncPhiQuad(ST,DST,Nord,Idec,N, PhiQuad,DPhiQuad)
c
      implicit none
      logical,          intent(in)  ::                          Idec(2)
      integer,          intent(in)  ::                        Nord(2),N
      double precision, intent(in)  ::     ST(1:2,0:1),DST(1:N,1:2,0:1)
      double precision, intent(out) ::     PhiQuad(2:Nord(1),2:Nord(2)),
     .                                DPhiQuad(1:N,2:Nord(1),2:Nord(2))
      integer ::                                minI,maxI,minJ,maxJ,i,j
      double precision ::        phiES(2:Nord(1)),DphiES(1:N,2:Nord(1)),
     .                           phiET(2:Nord(2)),DphiET(1:N,2:Nord(2))
c
c  ...local parameters
      minI = 2; maxI = Nord(1)
      minJ = 2; maxJ = Nord(2)
c
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncPhiQuad: N = ',i2)
      endif
c
c  ...get PhiE for each coordinate pair
      call AncPhiE(ST(1,0:1),DST(1:N,1,0:1),Nord(1),Idec(1),N,
     .          phiES,DphiES)
      call AncPhiE(ST(2,0:1),DST(1:N,2,0:1),Nord(2),Idec(2),N,
     .          phiET,DphiET)
c
c  ...the final result is the product of the two phiE
      do j=minJ,maxJ
        do i=minI,maxI
          PhiQuad(i,j) = phiES(i)*phiET(j)
          DphiQuad(1:N,i,j) = phiES(i)*DphiET(1:N,j)
     .                      + phiET(j)*DphiES(1:N,i)
        enddo
      enddo
c
      end subroutine AncPhiQuad
c
c----------------------------------------------------------------------
c
c     routine name      - AncEQuad
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - compute quadrilateral face Hcurl ancillary
c                         functions and their curls
c
c     arguments:
c
c     in:
c             ST       - affine coordinates associated to face
c                        2x2 matrix [s0,s1;t0,t1]
c             DST      - gradients of ST
c             Nord     - (NordS,NordT) vector polynomial order
c             Idec     - (IdecS,IdecT) vector binary flag:
c                        IdecS = FALSE if s0+s1 != 1
c                        IdecS = TRUE  if s0+s1  = 1, same with IdecT
c             N        - spatial dimension
c
c     out:
c             EQuad     - quad Hcurl ancillary functions
c             CurlEQuad - curls of quad Hcurl ancillary functions
c
c----------------------------------------------------------------------
c
      subroutine AncEQuad(ST,DST,Nord,Idec,N, EQuad,CurlEQuad)
c
      implicit none
      logical,          intent(in)  ::                          Idec(2)
      integer,          intent(in)  ::                        Nord(2),N
      double precision, intent(in)  ::     ST(1:2,0:1),DST(1:N,1:2,0:1)
      double precision, intent(out) :: EQuad(1:N,0:Nord(1)-1,2:Nord(2)),
     .                         CurlEQuad(1:2*N-3,0:Nord(1)-1,2:Nord(2))
      integer ::                          minI,maxI,minJ,maxJ,Ncurl,i,j
      double precision ::         S(0:1),DS(1:N,0:1),T(0:1),DT(1:N,0:1),
     .                EES(1:N,0:Nord(1)-1),curlEES(1:2*N-3,0:Nord(1)-1),
     .       phiET(2:Nord(2)),DphiET(1:N,2:Nord(2)),DphiETxEES(1:2*N-3)
c
c  ...local parameters
      minI = 0; maxI = Nord(1)-1
      minJ = 2; maxJ = Nord(2)
      Ncurl = 2*N-3
c
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncEQuad: N = ',i2)
      endif
c
      S = ST(1,0:1); DS = DST(1:N,1,0:1)
      T = ST(2,0:1); DT = DST(1:N,2,0:1)
c
      call AncEE(S,DS,Nord(1),Idec(1),N, EES,curlEES)
      call AncphiE(T,DT,Nord(2),Idec(2),N, phiET,DphiET)
c
      do j=minJ,maxJ
        do i=minI,maxI
          EQuad(1:N,i,j) = EES(1:N,i)*phiET(j)
c
          call cross(N,DphiET(1:N,j),EES(1:N,i), DphiETxEES)
c
          CurlEQuad(1:Ncurl,i,j) = curlEES(1:Ncurl,i)*phiET(j)
     .                           + DphiETxEES(1:Ncurl)
        enddo
      enddo

      end subroutine AncEQuad
c
c----------------------------------------------------------------------
c
c     routine name      - AncVQuad
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - compute quadrilateral face Hdiv ancillary
c                         functions and their divergences
c
c     arguments:
c
c     in:
c             ST       - affine coordinates associated to face
c                        2x2 matrix [s0,s1;t0,t1]
c             DST      - gradients of ST
c             Nord     - (NordS,NordT) vector polynomial order
c             Idec     - (IdecS,IdecT) vector binary flag:
c                        IdecS = FALSE if s0+s1 != 1
c                        IdecS = TRUE  if s0+s1  = 1, same with IdecT
c             N        - spatial dimension
c
c     out:
c             VQuad    - quad Hdiv ancillary functions
c             DivVQuad - divs of quad Hdiv ancillary functions
c
c----------------------------------------------------------------------
c
      subroutine AncVQuad(ST,DST,Nord,Idec,N, VQuad,DivVQuad)
c
      implicit none
      logical,          intent(in)  ::                          Idec(2)
      integer,          intent(in)  ::                        Nord(2),N
      double precision, intent(in)  ::     ST(1:2,0:1),DST(1:N,1:2,0:1)
      double precision, intent(out) ::
     .                               VQuad(1:N,0:Nord(1)-1,0:Nord(2)-1),
     .                                DivVQuad(0:Nord(1)-1,0:Nord(2)-1)
      integer ::                          minI,maxI,minJ,maxJ,Ncurl,i,j
      double precision ::         S(0:1),T(0:1),DS(1:N,0:1),DT(1:N,0:1),
     .          EES(1:N,0:Nord(1)-1),curlEES(1:2*N-3,0:Nord(1)-1),prod1,
     .          EET(1:N,0:Nord(2)-1),curlEET(1:2*N-3,0:Nord(2)-1),prod2
c
c  ...local parameters
      minI = 0; maxI = Nord(1)-1
      minJ = 0; maxJ = Nord(2)-1
      Ncurl = 2*N-3
c
      if (N.lt.3) then
        write(*,7001) N
 7001   format('AncVQuad: N = ',i2)
      endif
c
      S = ST(1,0:1); DS = DST(1:N,1,0:1)
      T = ST(2,0:1); DT = DST(1:N,2,0:1)

      call AncEE(S,DS,Nord(1),Idec(1),N, EES,curlEES)
      call AncEE(T,DT,Nord(2),Idec(2),N, EET,curlEET)
c
c      ...slight speedup when Idec=(.TRUE.,.TRUE.)
      if (Idec(1).and.Idec(2)) then
        do j=minJ,maxJ
          do i=minI,maxI
            call cross(N,EES(1:N,i),EET(1:N,j), VQuad(1:N,i,j))
          enddo
        enddo
        DivVQuad(minI:maxI,minJ:maxJ) = 0.d0
      else
        do j=minJ,maxJ
          do i=minI,maxI
            call cross(N,EES(1:N,i),EET(1:N,j), VQuad(1:N,i,j))
c
            call dot_product(EET(1:N,j),curlEES(1:N,i), prod1)
            call dot_product(EES(1:N,i),curlEET(1:N,j), prod2)
c
            DivVQuad(i,j) = prod1-prod2
          enddo
        enddo
      endif
c
      end subroutine AncVQuad
c
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c                           TRIANGULAR FACES
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
c     routine name      - AncPhiTri
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - compute triangle face H1 ancillary
c                         functions and their gradients
c
c     arguments:
c
c     in:
c             S        - (s0,s1,s2) affine coordinates associated to
c                        triangle face
c             DS       - derivatives of S0,S1,S2
c             Nord     - polynomial order
c             Idec     - Binary flag:
c                        = FALSE s0+s1+s2 != 1
c                        = TRUE  s0+s1+s2  = 1
c             N        - spatial dimension
c
c     out:
c             PhiTri   - triangle H1 ancillary functions
c             DPhiTri  - grads of triangle H1 ancillary functions
c
c----------------------------------------------------------------------
c
      subroutine AncPhiTri(S,DS,Nord,Idec,N, PhiTri,DPhiTri)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:2),DS(1:N,0:2)
      double precision, intent(out) ::        PhiTri(2:Nord-1,1:Nord-2),
     .                                   DPhiTri(1:N,2:Nord-1,1:Nord-2)
      logical ::                                                  IdecE
      integer ::       minI,maxI,minJ,maxJ,minIJ,maxIJ,minalpha,i,j,nij
      double precision ::                          sL(0:1),DsL(1:N,0:1),
     .     phiE(2:Nord-1),DphiE(1:N,2:Nord-1),homLal(2:Nord-1,1:Nord-2),
     .                                   DhomLal(1:N,2:Nord-1,1:Nord-2)
c
c  ...local parameters
      minI = 2; maxI = Nord-1
      minJ = 1; maxJ = Nord-2
      minIJ = minI+minJ; maxIJ = Nord
      minalpha = 2*minI
      IdecE = .FALSE.
c
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncPhiTri: N = ',i2)
      endif
c
c  ...get PhiE - this is never a simplified case (IdecE=0)
      call AncPhiE(S(0:1),DS(1:N,0:1),Nord-minJ,IdecE,N, phiE,DphiE)
c
c  ...get homogenized Jacobi polynomials, homLal, and gradients
      sL(0) = S(0)+S(1); sL(1) = S(2)
      DsL(1:N,0) = DS(1:N,0)+DS(1:N,1)
      DsL(1:N,1) = DS(1:N,2)
      call HomIJacobi(sL,DsL,maxJ,minalpha,Idec,N, homLal,DhomLal)
c
c  ...simply complete the required information
      do nij=minIJ,maxIJ
        do i=minI,nij-minJ
          j=nij-i
          PhiTri(i,j) = phiE(i)*homLal(i,j)
          DPhiTri(1:N,i,j) = homLal(i,j)*DphiE(1:N,i)
     .                        + phiE(i)*DhomLal(1:N,i,j)
        enddo
      enddo
c
      end subroutine AncPhiTri
c
c----------------------------------------------------------------------
c
c     routine name      - AncETri
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - compute triangle face Hcurl ancillary
c                         functions and their curls
c
c     arguments:
c
c     in:
c             S        - (s0,s1,s2) affine coordinates associated to
c                        triangle face
c             DS       - derivatives of S0,S1,S2
c             Nord     - polynomial order
c             Idec     - Binary flag:
c                        = FALSE s0+s1+s2 != 1
c                        = TRUE  s0+s1+s2  = 1
c             N        - spatial dimension
c
c     out:
c             ETri     - triangle Hcurl ancillary functions
c             CurlETri - curls of triangle Hcurl ancillary functions
c
c----------------------------------------------------------------------
c
      subroutine AncETri(S,DS,Nord,Idec,N, ETri,CurlETri)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:2),DS(1:N,0:2)
      double precision, intent(out) ::      ETri(1:N,0:Nord-2,1:Nord-1),
     .                              CurlETri(1:2*N-3,0:Nord-2,1:Nord-1)
      logical ::                                                  IdecE
      integer :: minI,maxI,minJ,maxJ,minIJ,maxIJ,minalpha,Ncurl,i,j,nij
      double precision ::     EE(1:N,0:Nord-2),curlEE(1:2*N-3,0:Nord-2),
     .                   sL(0:1),DsL(1:N,0:1),homLal(0:Nord-2,1:Nord-1),
     .               DhomLal(1:N,0:Nord-2,1:Nord-1),DhomLalxEE(1:2*N-3)
c
c  ...local parameters
      minI = 0; maxI = Nord-2
      minJ = 1; maxJ = Nord-1
      minIJ = minI+minJ; maxIJ = Nord-1
      minalpha = 2*minI+1
      Ncurl = 2*N-3
      IdecE = .FALSE.
c
      if (N.lt.2) then
        write(*,7001) N
 7001   format('AncETri: N = ',i2)
      endif
c
c  ...get EE - this is never a simplified case (IdecE=0)
      call AncEE(S(0:1),DS(1:N,0:1),Nord-minJ,IdecE,N, EE,curlEE)
c
c  ...get homogenized Integrated Jacobi polynomials, homLal, and gradients
      sL(0) = S(0)+S(1); sL(1) = S(2)
      DsL(1:N,0) = DS(1:N,0)+DS(1:N,1)
      DsL(1:N,1) = DS(1:N,2)
      call HomIJacobi(sL,DsL,maxJ,minalpha,Idec,N, homLal,DhomLal)
c
c  ...simply complete the required information
      do nij=minIJ,maxIJ
        do i=minI,nij-minJ
          j=nij-i
            ETri(1:N,i,j) = EE(1:N,i)*homLal(i,j)
c
            call cross(N,DhomLal(1:N,i,j),EE(1:N,i), DhomLalxEE)
c
            CurlETri(1:Ncurl,i,j) = homLal(i,j)*curlEE(1:Ncurl,i)
     .                            + DhomLalxEE
        enddo
      enddo
c
      end subroutine AncETri
c
c----------------------------------------------------------------------
c
c     routine name      - AncVTri
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - compute triangle face Hcurl ancillary
c                         functions and their curls (family II)
c
c     arguments:
c
c     in:
c             S        - (s0,s1,s2) affine coordinates associated to
c                        triangle face
c             DS       - derivatives of S0,S1,S2
c             Nord     - polynomial order
c             Idec     - Binary flag:
c                        = FALSE s0+s1+s2 != 1
c                        = TRUE  s0+s1+s2  = 1
c             N        - spatial dimension
c
c     out:
c             VTri     - triangle Hdiv ancillary functions
c             DivVTri  - divs of triangle Hdiv ancillary functions
c
c----------------------------------------------------------------------
c
      subroutine AncVTri(S,DS,Nord,Idec,N, VTri,DivVTri)
c
      implicit none
      logical,          intent(in)  ::                             Idec
      integer,          intent(in)  ::                           Nord,N
      double precision, intent(in)  ::               S(0:2),DS(1:N,0:2)
      double precision, intent(out) ::      VTri(1:N,0:Nord-1,0:Nord-1),
     .                                       DivVTri(0:Nord-1,0:Nord-1)
      integer ::       minI,maxI,minJ,maxJ,minIJ,maxIJ,minalpha,i,j,nij
      double precision ::      homP(0:Nord-1),homPal(0:Nord-1,0:Nord-1),
     .        DS0xDS1(N),DS1xDS2(N),DS2xDS0(N),V00(N),tripleprod,psiTri
c
c  ...local parameters
      minI = 0; maxI = Nord-1
      minJ = 0; maxJ = Nord-1
      minIJ = minI+minJ; maxIJ = Nord-1
      minalpha = 2*minI+1
c
      if (N.lt.3) then
        write(*,7001) N
 7001   format('AncVTri: N = ',i2)
      endif
c
c  ...get homogenized Legendre polynomials, homP
      call HomLegendre(S(0:1),Nord-1-minJ, homP)
c
c  ...get homogenized Jacobi polynomials, homPal
      call HomJacobi((/S(0)+S(1),S(2)/),maxJ,minalpha, homPal)
c
c  ...simplified case
      if (Idec) then
c    ...construct V00
        call cross(N,DS(1:N,1),DS(1:N,2), V00)
c    ...loop
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            VTri(1:N,i,j) = homP(i)*homPal(i,j)*V00(1:N)
          enddo
        enddo
c
        DivVTri = 0.d0
c
c  ...general case
      else
c    ...construct V00
        call cross(N,DS(1:N,0),DS(1:N,1), DS0xDS1)
        call cross(N,DS(1:N,1),DS(1:N,2), DS1xDS2)
        call cross(N,DS(1:N,2),DS(1:N,0), DS2xDS0)
        V00 = S(0)*DS1xDS2+S(1)*DS2xDS0+S(2)*DS0xDS1
c    ...loop
        do nij=minIJ,maxIJ
          do i=minI,nij-minJ
            j=nij-i
            psiTri = homP(i)*homPal(i,j)
c
            VTri(1:N,i,j) = psiTri*V00
            DivVTri(i,j)  = (nij+3)*psiTri
          enddo
        enddo
c
        call dot_product(DS(1:N,0),DS1xDS2, tripleprod)
c
        DivVTri = DivVTri*tripleprod
      endif
c
      end subroutine AncVTri
c
